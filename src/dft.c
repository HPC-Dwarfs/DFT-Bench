/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of DFT-Bench.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#include <fftw3.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "allocate.h"
#include "blasLapack.h"
#include "comm.h"
#include "constants.h"
#include "dft.h"
#include "eigen.h"
#include "gtoOnGrid.h"
#include "likwid-marker.h"
#include "poscarReader.h"
#include "timing.h"
#include "util.h"

enum {
  PH_GTO,
  PH_TRANSPOSE,
  PH_OVERLAP,
  PH_REDUCE,
  PH_EIGEN,
  PH_GATHER,
  PH_ROTATE,
  PH_UNTRANSPOSE,
  PH_KINETIC,
  NPHASES
};
static const char *const PHASE_NAMES[NPHASES] = { "gto",
  "transpose",
  "overlap",
  "reduce",
  "eigen",
  "gather",
  "rotate",
  "untranspose",
  "kinetic" };
/* Reported phases: all in the MPI build, the original five in the sequential build */
static const int MPI_PHASES[NPHASES] = { PH_GTO,
  PH_TRANSPOSE,
  PH_OVERLAP,
  PH_REDUCE,
  PH_EIGEN,
  PH_GATHER,
  PH_ROTATE,
  PH_UNTRANSPOSE,
  PH_KINETIC };
static const int SEQUENTIAL_PHASES[] = {
  PH_GTO, PH_OVERLAP, PH_EIGEN, PH_ROTATE, PH_KINETIC
};

/* How the transposed form is held: not at all (single group, it is the orbital form),
 * completely (stored), or one chunk of rows at a time (streamed) */
enum { MODE_NONE, MODE_STORED, MODE_STREAMED };
static const char *const MODE_NAMES[] = { "none (single group)", "stored", "streamed" };

#define GB (1024.0 * 1024.0 * 1024.0)

typedef struct {
  fftw_plan planf, planb;
  double *rbuf; /* T|psi> in real space */
  fftw_complex *cbuf; /* psi in reciprocal space */
  double *k2[3]; /* |k|^2 components along each direction */
} KineticWorkType;

typedef struct {
  const CommType *comm;
  GridType grid;
  ElectronStateType es;
  DistType dist;
  EigenType eig;
  int mode;
  int nchunks;
  size_t *choff; /* row offsets of the group's ranks in the current transposed chunk */
  double *sendbuf, *recvbuf; /* pack buffers of the transposes */
  size_t *scounts, *sdispls, *rcounts, *rdispls;
  double *overlap, *xmat, *eval; /* group-shared S, X and eigenvalues */
  double *wa; /* rotation work, wrows x norb */
  size_t wrows;
  double *segbuf; /* packed S segment of the lane reduction */
  CommWindowType winPsi, winPsit, winOverlap, winXmat, winEval;
  bool timed; /* record phase times, byte counts and LIKWID regions */
  double t[NPHASES]; /* phase times of the current repetition */
  double bytes[NPHASES]; /* bytes received from other groups per repetition */
} DftType;

static void *allocateArray(size_t bytes)
{
  return allocate(ARRAY_ALIGNMENT, bytes > 0 ? bytes : ARRAY_ALIGNMENT);
}

/* Starts timing a phase without a barrier (streamed mode accumulates chunks) */
static double phaseResume(const DftType *d, int phase)
{
  (void)phase;
  if (d->timed) {
    LIKWID_MARKER_START(PHASE_NAMES[phase]);
  }
  return getTimeStamp();
}

/* Starts timing a phase after a barrier over all ranks */
static double phaseStart(const DftType *d, int phase)
{
  commBarrier(COMM_WORLD);
  return phaseResume(d, phase);
}

static void phaseStop(DftType *d, int phase, double t0)
{
  double t = getTimeStamp() - t0;
  if (d->timed) {
    LIKWID_MARKER_STOP(PHASE_NAMES[phase]);
    d->t[phase] += t;
  }
}

static void addBytes(DftType *d, int phase, double bytes)
{
  if (d->timed) {
    d->bytes[phase] += bytes;
  }
}

/* Smallest even n' >= n whose only prime factors are 2, 3 and 5 */
static int fftGoodSize(int n)
{
  for (int m = MAX(n, 2);; m++) {
    int r = m;
    if (r % 2)
      continue;
    while (r % 2 == 0)
      r /= 2;
    while (r % 3 == 0)
      r /= 3;
    while (r % 5 == 0)
      r /= 5;
    if (r == 1)
      return m;
  }
}

static size_t round8(size_t n)
{
  return (n + 7) / 8 * 8;
}

static void setupGrid(const PoscarFileType *pf, const ConfigType *cfg, GridType *grid)
{
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i != j && fabs(pf->cell[i][j]) > 1.0E-10 * fabs(pf->cell[i][i])) {
        dftAbort("ERROR: only orthogonal cells are supported\n");
      }
    }
    if (pf->cell[i][i] <= 0.0) {
      dftAbort("ERROR: cell vector %d must point in positive direction\n", i + 1);
    }
  }

  grid->ngtot = 1;
  for (int d = 0; d < 3; d++) {
    double alat = pf->cell[d][d];
    if (cfg->ng[0] > 0) {
      grid->ng[d] = cfg->ng[d];
    } else {
      grid->ng[d] = fftGoodSize((int)ceil(alat / cfg->spacing));
    }
    grid->h[d]      = alat / grid->ng[d];
    grid->origin[d] = 0.0;
    grid->ngtot *= (size_t)grid->ng[d];
  }
  // pad the columns of the orbital block to keep every orbital aligned
  size_t pad  = ARRAY_ALIGNMENT / sizeof(double);
  grid->ld    = (grid->ngtot + pad - 1) / pad * pad;
  grid->voxel = grid->h[0] * grid->h[1] * grid->h[2];

  if (grid->ld > INT_MAX) {
    dftAbort(
        "ERROR: %zu grid points exceed the 32 bit integer BLAS interface\n", grid->ngtot);
  }
}

/* Splits orbitals and grid rows over the ranks in group order, see DistType */
static void setupDistribution(DftType *d)
{
  const CommType *c = d->comm;
  DistType *ds      = &d->dist;
  int nproc         = c->size;
  int norb          = d->es.norb;
  size_t ngtot      = d->grid.ngtot;
  size_t units      = (ngtot + 7) / 8;
  int *order        = malloc(nproc * sizeof(int));
  int *first        = calloc(c->ngroups + 1, sizeof(int));

  for (int g = 0; g < c->ngroups; g++) {
    first[g + 1] = first[g] + c->groupSize[g];
  }
  for (int p = 0; p < nproc; p++) {
    order[first[c->groupOf[p]] + c->grankOf[p]] = p;
  }

  ds->orb0   = malloc(nproc * sizeof(int));
  ds->norb   = malloc(nproc * sizeof(int));
  ds->row0   = malloc(nproc * sizeof(size_t));
  ds->nrows  = malloc(nproc * sizeof(size_t));
  ds->gorb0  = calloc(c->ngroups, sizeof(int));
  ds->gnorb  = calloc(c->ngroups, sizeof(int));
  ds->grow0  = calloc(c->ngroups, sizeof(size_t));
  ds->gnrows = calloc(c->ngroups, sizeof(size_t));
  ds->idle   = 0;

  int o      = 0;
  size_t r   = 0;
  for (int i = 0; i < nproc; i++) {
    int p        = order[i];
    int g        = c->groupOf[p];
    size_t nu    = units / nproc + ((size_t)i < units % nproc);
    ds->orb0[p]  = o;
    ds->norb[p]  = norb / nproc + (i < norb % nproc);
    ds->row0[p]  = MIN(r, ngtot);
    ds->nrows[p] = MIN(r + 8 * nu, ngtot) - ds->row0[p];
    o += ds->norb[p];
    r += 8 * nu;
    if (c->grankOf[p] == 0) {
      ds->gorb0[g] = ds->orb0[p];
      ds->grow0[g] = ds->row0[p];
    }
    ds->gnorb[g] += ds->norb[p];
    ds->gnrows[g] += ds->nrows[p];
    ds->idle += (ds->norb[p] == 0);
  }

  ds->members = malloc(c->gsize * sizeof(int));
  memcpy(ds->members, order + first[c->group], c->gsize * sizeof(int));
  free(order);
  free(first);
}

static void freeDistribution(DistType *ds)
{
  free(ds->orb0);
  free(ds->norb);
  free(ds->row0);
  free(ds->nrows);
  free(ds->gorb0);
  free(ds->gnorb);
  free(ds->grow0);
  free(ds->gnrows);
  free(ds->members);
}

/* Rows [start, start + len) of rank p's row range that belong to chunk k. Chunks
 * cover the same fraction of every rank's range, in multiples of 8 rows. */
static void chunkRows(const DftType *d, int p, int k, size_t *start, size_t *len)
{
  size_t n     = d->dist.nrows[p];
  size_t units = (n + 7) / 8;
  size_t a     = MIN(units * k / d->nchunks * 8, n);
  size_t b     = MIN(units * (k + 1) / d->nchunks * 8, n);
  *start       = a;
  *len         = b - a;
}

/* Row offsets of the group's ranks in the transposed form for chunk k. Returns the
 * number of rows of the group in the chunk. */
static size_t chunkOffsets(DftType *d, int k)
{
  const DistType *ds = &d->dist;
  int g              = d->comm->group;
  size_t rows        = 0;

  for (int r = 0; r < d->comm->gsize; r++) {
    int q = ds->members[r];
    size_t s, n;
    chunkRows(d, q, k, &s, &n);
    d->choff[r] = (d->mode == MODE_STREAMED) ? rows : ds->row0[q] - ds->grow0[g] + s;
    rows += n;
  }
  return rows;
}

static double physicalMemory(void)
{
  long pages = sysconf(_SC_PHYS_PAGES);
  long size  = sysconf(_SC_PAGESIZE);
  return (pages > 0 && size > 0) ? (double)pages * size : 1.0E30;
}

/* Pack buffer bytes of the transposes of group g, without chunking */
static double bufferBytes(const DftType *d, int g)
{
  const DistType *ds = &d->dist;
  double bytes       = 0.0;
  for (int p = 0; p < d->comm->size; p++) {
    if (d->comm->groupOf[p] == g) {
      double send = (double)ds->norb[p] * (d->grid.ngtot - ds->gnrows[g]);
      double recv = (double)ds->nrows[p] * (d->es.norb - ds->gnorb[g]);
      bytes += 2.0 * MAX(send, recv) * sizeof(double);
    }
  }
  return bytes;
}

/* Memory per group: orbital form + S, X, eigenvalues, plus transposed form and pack
 * buffers, which shrink with the chunk count in the streamed mode */
static double groupMemory(const DftType *d, int g, int mode, int nchunks)
{
  double norb = d->es.norb;
  double fixed =
      ((double)d->grid.ld * d->dist.gnorb[g] + 2.0 * norb * norb + norb) * sizeof(double);
  if (mode == MODE_NONE) {
    return fixed;
  }
  double tr  = (double)round8(d->dist.gnrows[g]) * norb * sizeof(double);
  double buf = bufferBytes(d, g) / nchunks;
  return fixed + buf + (mode == MODE_STORED ? tr : tr / nchunks);
}

/* Chooses the transpose mode and chunk count within the memory budget per group */
static void planMemory(DftType *d, const ConfigType *cfg, double *budget, double *need)
{
  const CommType *c = d->comm;
  double b          = (cfg->memBudget > 0.0)
                          ? cfg->memBudget * GB
                          : DEFAULT_MEMORY_FRACTION * physicalMemory() / c->nodeGroups;
  // nodes may differ, all ranks use the smallest budget
  commAllreduce(&b, 1, COMM_MIN, COMM_WORLD);

  size_t maxRows = 0;
  for (int p = 0; p < c->size; p++) {
    maxRows = MAX(maxRows, d->dist.nrows[p]);
  }
  int kmax   = (int)MAX((size_t)1, (maxRows + 7) / 8);

  d->mode    = (c->ngroups > 1) ? MODE_STORED : MODE_NONE;
  d->nchunks = 1;
  for (int g = 0; g < c->ngroups; g++) {
    if (d->mode == MODE_STORED && groupMemory(d, g, MODE_STORED, kmax) > b) {
      d->mode = MODE_STREAMED;
    }
  }
  double required = 0.0;
  for (int g = 0; g < c->ngroups; g++) {
    required = MAX(required, groupMemory(d, g, d->mode, kmax));
  }
  if (required > b) {
    dftAbort("ERROR: the memory budget of %.2f GB per group is too small, the orbitals "
             "need at least %.2f GB per group\n",
        b / GB,
        required / GB);
  }
  // smallest chunk count that fits every group, the memory decreases with it
  for (int g = 0; g < c->ngroups; g++) {
    int lo = d->nchunks, hi = kmax;
    while (lo < hi) {
      int mid = lo + (hi - lo) / 2;
      if (groupMemory(d, g, d->mode, mid) > b) {
        lo = mid + 1;
      } else {
        hi = mid;
      }
    }
    d->nchunks = lo;
  }
  *need = 0.0;
  for (int g = 0; g < c->ngroups; g++) {
    *need = MAX(*need, groupMemory(d, g, d->mode, d->nchunks));
  }
  *budget = b;
}

/* Allocates the group-shared windows and the buffers of the transposes */
static void allocateState(DftType *d)
{
  const CommType *c     = d->comm;
  const DistType *ds    = &d->dist;
  ElectronStateType *es = &d->es;
  int norb              = es->norb;
  int me                = c->rank;

  es->psi =
      commAllocShared((size_t)ds->norb[me] * d->grid.ld * sizeof(double), &d->winPsi);
  d->choff = malloc(c->gsize * sizeof(size_t));

  if (d->mode == MODE_NONE) {
    es->psit = es->psi;
    es->ldt  = d->grid.ld;
  } else {
    // stored: all rows of the group, streamed: the largest chunk
    size_t rows = ds->gnrows[c->group];
    if (d->mode == MODE_STREAMED) {
      rows = 0;
      for (int k = 0; k < d->nchunks; k++) {
        rows = MAX(rows, chunkOffsets(d, k));
      }
    }
    es->ldt  = MAX((size_t)8, round8(rows));
    es->psit = commAllocShared(
        (c->grank == 0) ? es->ldt * norb * sizeof(double) : 0, &d->winPsit);

    // first touch: every rank zeroes its rows (of chunk 0), the last one the padding
    size_t s, n;
    chunkOffsets(d, 0);
    chunkRows(d, me, 0, &s, &n);
    size_t r0 = d->choff[c->grank];
    size_t r1 = r0 + ((d->mode == MODE_STREAMED) ? n : ds->nrows[me]);
    if (d->mode == MODE_STORED) {
      r0 = ds->row0[me] - ds->grow0[c->group];
      r1 = r0 + ds->nrows[me];
    }
    if (c->grank == c->gsize - 1) {
      r1 = es->ldt;
    }
    for (int j = 0; j < norb; j++) {
      memset(es->psit + (size_t)j * es->ldt + r0, 0, (r1 - r0) * sizeof(double));
    }

    // pack buffers: the largest exchange of any chunk in either direction
    size_t outRows = 0, maxBuf = 0;
    int outOrbs = norb - ds->gnorb[c->group];
    for (int k = 0; k < d->nchunks; k++) {
      outRows = 0;
      for (int q = 0; q < c->size; q++) {
        if (c->groupOf[q] != c->group) {
          chunkRows(d, q, k, &s, &n);
          outRows += n;
        }
      }
      chunkRows(d, me, k, &s, &n);
      maxBuf = MAX(maxBuf, MAX(outRows * ds->norb[me], n * outOrbs));
    }
    d->sendbuf = allocateArray(maxBuf * sizeof(double));
    d->recvbuf = allocateArray(maxBuf * sizeof(double));
    d->scounts = malloc(c->size * sizeof(size_t));
    d->sdispls = malloc(c->size * sizeof(size_t));
    d->rcounts = malloc(c->size * sizeof(size_t));
    d->rdispls = malloc(c->size * sizeof(size_t));
  }

  size_t nn  = (size_t)norb * norb;
  d->overlap = commAllocShared((c->grank == 0) ? nn * sizeof(double) : 0, &d->winOverlap);
  d->xmat    = commAllocShared((c->grank == 0) ? nn * sizeof(double) : 0, &d->winXmat);
  d->eval    = commAllocShared((c->grank == 0) ? norb * sizeof(double) : 0, &d->winEval);
  commSyncGroup();

  // rotation work, the rank's rows in blocks of at most ROTATE_BLOCK
  d->wrows  = (commIsMpi()) ? MAX((size_t)1, MIN((size_t)ROTATE_BLOCK, ds->nrows[me]))
                            : ROTATE_BLOCK;
  d->wa     = allocateArray(d->wrows * norb * sizeof(double));
  d->segbuf = NULL;
  if (c->ngroups > 1 && c->grank < c->lmin) {
    d->segbuf = allocateArray(nn * sizeof(double) / c->lmin + norb * sizeof(double));
  }
}

static void freeState(DftType *d)
{
  if (d->mode != MODE_NONE) {
    commFreeShared(&d->winPsit, d->es.psit);
    free(d->sendbuf);
    free(d->recvbuf);
    free(d->scounts);
    free(d->sdispls);
    free(d->rcounts);
    free(d->rdispls);
  }
  commFreeShared(&d->winPsi, d->es.psi);
  commFreeShared(&d->winOverlap, d->overlap);
  commFreeShared(&d->winXmat, d->xmat);
  commFreeShared(&d->winEval, d->eval);
  free(d->choff);
  free(d->wa);
  free(d->segbuf);
}

/* First orbital of this rank in the group's orbital form */
static double *rankOrbitals(const DftType *d)
{
  const DistType *ds = &d->dist;
  int me             = d->comm->rank;
  return d->es.psi + (size_t)(ds->orb0[me] - ds->gorb0[d->comm->group]) * d->grid.ld;
}

/* Orbital form -> transposed form for the rows of chunk k. Rows of the own group are
 * copied in shared memory, rows of other groups are exchanged with one all-to-all.
 * Returns the number of rows of the group in the chunk. */
static size_t transposeChunk(DftType *d, int k, int phase)
{
  const CommType *c  = d->comm;
  const DistType *ds = &d->dist;
  int me = c->rank, g = c->group;
  int nmine        = ds->norb[me];
  size_t ld        = d->grid.ld;
  size_t ldt       = d->es.ldt;
  double *t        = d->es.psit;
  const double *pm = (nmine > 0) ? rankOrbitals(d) : NULL;
  size_t rows      = chunkOffsets(d, k);
  size_t s, n, pos = 0;

  for (int r = 0; r < c->gsize; r++) {
    int q = ds->members[r];
    chunkRows(d, q, k, &s, &n);
    for (int i = 0; i < nmine; i++) {
      memcpy(t + (size_t)(ds->orb0[me] + i) * ldt + d->choff[r],
          pm + i * ld + ds->row0[q] + s,
          n * sizeof(double));
    }
  }

  for (int q = 0; q < c->size; q++) {
    d->sdispls[q] = pos;
    d->scounts[q] = 0;
    if (c->groupOf[q] != g) {
      chunkRows(d, q, k, &s, &n);
      for (int i = 0; i < nmine; i++) {
        memcpy(d->sendbuf + pos, pm + i * ld + ds->row0[q] + s, n * sizeof(double));
        pos += n;
      }
      d->scounts[q] = n * nmine;
    }
  }
  chunkRows(d, me, k, &s, &n);
  pos = 0;
  for (int p = 0; p < c->size; p++) {
    d->rdispls[p] = pos;
    d->rcounts[p] = (c->groupOf[p] != g) ? n * ds->norb[p] : 0;
    pos += d->rcounts[p];
  }
  addBytes(d, phase, (double)pos * sizeof(double));
  commAlltoallv(d->sendbuf, d->scounts, d->sdispls, d->recvbuf, d->rcounts, d->rdispls);

  for (int p = 0; p < c->size; p++) {
    if (d->rcounts[p] == 0) {
      continue;
    }
    for (int i = 0; i < ds->norb[p]; i++) {
      memcpy(t + (size_t)(ds->orb0[p] + i) * ldt + d->choff[c->grank],
          d->recvbuf + d->rdispls[p] + i * n,
          n * sizeof(double));
    }
  }
  commSyncGroup();
  return rows;
}

/* Transposed form -> orbital form for the rows of chunk k, the reverse of
 * transposeChunk */
static void untransposeChunk(DftType *d, int k, int phase)
{
  const CommType *c  = d->comm;
  const DistType *ds = &d->dist;
  int me = c->rank, g = c->group;
  int nmine       = ds->norb[me];
  size_t ld       = d->grid.ld;
  size_t ldt      = d->es.ldt;
  const double *t = d->es.psit;
  double *pm      = (nmine > 0) ? rankOrbitals(d) : NULL;
  size_t s, n, pos = 0;

  chunkOffsets(d, k);
  for (int r = 0; r < c->gsize; r++) {
    int q = ds->members[r];
    chunkRows(d, q, k, &s, &n);
    for (int i = 0; i < nmine; i++) {
      memcpy(pm + i * ld + ds->row0[q] + s,
          t + (size_t)(ds->orb0[me] + i) * ldt + d->choff[r],
          n * sizeof(double));
    }
  }

  chunkRows(d, me, k, &s, &n);
  for (int p = 0; p < c->size; p++) {
    d->sdispls[p] = pos;
    d->scounts[p] = 0;
    if (c->groupOf[p] != g && n > 0) {
      for (int i = 0; i < ds->norb[p]; i++) {
        memcpy(d->sendbuf + pos,
            t + (size_t)(ds->orb0[p] + i) * ldt + d->choff[c->grank],
            n * sizeof(double));
        pos += n;
      }
      d->scounts[p] = n * ds->norb[p];
    }
  }
  pos = 0;
  for (int q = 0; q < c->size; q++) {
    d->rdispls[q] = pos;
    d->rcounts[q] = 0;
    if (c->groupOf[q] != g) {
      chunkRows(d, q, k, &s, &n);
      d->rcounts[q] = n * nmine;
    }
    pos += d->rcounts[q];
  }
  addBytes(d, phase, (double)pos * sizeof(double));
  commAlltoallv(d->sendbuf, d->scounts, d->sdispls, d->recvbuf, d->rcounts, d->rdispls);

  for (int q = 0; q < c->size; q++) {
    if (d->rcounts[q] > 0) {
      chunkRows(d, q, k, &s, &n);
      for (int i = 0; i < nmine; i++) {
        memcpy(pm + i * ld + ds->row0[q] + s,
            d->recvbuf + d->rdispls[q] + i * n,
            n * sizeof(double));
      }
    }
  }
  commSyncGroup();
}

static void placeOrbitals(DftType *d, const double *rat, double rgcut)
{
  const ElectronStateType *es = &d->es;
  const GridType *grid        = &d->grid;
  int o0                      = d->dist.orb0[d->comm->rank];
  int n                       = d->dist.norb[d->comm->rank];
  double *psi                 = (n > 0) ? rankOrbitals(d) : NULL;

  for (int i = 0; i < n; i++) {
    int iorb = o0 + i;
    putGtoSymOrthoOrb(&rat[3 * es->iat[iorb]],
        es->gw[iorb],
        rgcut,
        grid->origin,
        grid->ng,
        grid->h,
        psi + i * grid->ld);
  }
}

/* The FFT plans work on rank-private buffers and are executed on the orbitals */
static void kineticInit(KineticWorkType *kw, const GridType *grid)
{
  const int *ng = grid->ng;
  size_t nk     = (size_t)(ng[0] / 2 + 1) * ng[1] * ng[2];

  kw->rbuf      = allocateArray(grid->ngtot * sizeof(double));
  kw->cbuf      = allocateArray(nk * sizeof(fftw_complex));

  // x is stored contiguously, so it is the last (halved) dimension for FFTW
  kw->planf =
      fftw_plan_dft_r2c_3d(ng[2], ng[1], ng[0], kw->rbuf, kw->cbuf, FFTW_ESTIMATE);
  kw->planb =
      fftw_plan_dft_c2r_3d(ng[2], ng[1], ng[0], kw->cbuf, kw->rbuf, FFTW_ESTIMATE);

  for (int d = 0; d < 3; d++) {
    int n       = (d == 0) ? ng[0] / 2 + 1 : ng[d];
    double alat = ng[d] * grid->h[d];
    kw->k2[d]   = malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
      int m        = (i <= ng[d] / 2) ? i : i - ng[d];
      double k     = 2.0 * M_PI * m / alat;
      kw->k2[d][i] = k * k;
    }
  }
}

static void kineticFree(KineticWorkType *kw)
{
  fftw_destroy_plan(kw->planf);
  fftw_destroy_plan(kw->planb);
  free(kw->rbuf);
  free(kw->cbuf);
  for (int d = 0; d < 3; d++) {
    free(kw->k2[d]);
  }
}

/* <psi_i| -1/2 nabla^2 |psi_i> summed over all orbitals, using one FFT pair per orbital */
static double kineticEnergy(const KineticWorkType *kw, const DftType *d)
{
  const GridType *grid = &d->grid;
  const int *ng        = grid->ng;
  int nkx              = ng[0] / 2 + 1;
  double scale = 0.5 / (double)grid->ngtot; // 1/2 from T, 1/N from the c2r transform
  double ekin  = 0.0;
  double *tpsi = kw->rbuf;
  fftw_complex *cpsi = kw->cbuf;
  int n              = d->dist.norb[d->comm->rank];
  double *orbitals   = (n > 0) ? rankOrbitals(d) : NULL;

  for (int iorb = 0; iorb < n; iorb++) {
    double *psi = orbitals + iorb * grid->ld;
    double *in  = psi;
    // new-array execution needs the alignment of the planned buffer
    if (fftw_alignment_of(psi) != fftw_alignment_of(kw->rbuf)) {
      memcpy(kw->rbuf, psi, grid->ngtot * sizeof(double));
      in = kw->rbuf;
    }
    fftw_execute_dft_r2c(kw->planf, in, cpsi);
    for (int iz = 0; iz < ng[2]; iz++) {
      for (int iy = 0; iy < ng[1]; iy++) {
        double k2yz     = kw->k2[2][iz] + kw->k2[1][iy];
        fftw_complex *c = cpsi + ((size_t)iz * ng[1] + iy) * nkx;
        for (int ix = 0; ix < nkx; ix++) {
          double fac = scale * (k2yz + kw->k2[0][ix]);
          c[ix][0] *= fac;
          c[ix][1] *= fac;
        }
      }
    }
    fftw_execute_dft_c2r(kw->planb, cpsi, tpsi);
    double tt = 0.0;
    for (size_t i = 0; i < grid->ngtot; i++) {
      tt += psi[i] * tpsi[i];
    }
    ekin += tt * grid->voxel;
  }
  commAllreduce(&ekin, 1, COMM_SUM, COMM_WORLD);
  return ekin;
}

/* Upper triangle of S = voxel * Psi^T Psi (column-major, norb x norb) from the first m
 * rows of the transposed form. S is split into tiles, which the group's ranks compute
 * round-robin; beta = 1 accumulates. */
static void overlapTiles(DftType *d, size_t m, double beta)
{
  const ElectronStateType *es = &d->es;
  int norb                    = es->norb;
  int nt                      = (norb + OVERLAP_TILE - 1) / OVERLAP_TILE;
  int tile                    = 0;

  for (int jt = 0; jt < nt; jt++) {
    for (int it = 0; it <= jt; it++, tile++) {
      if (tile % d->comm->gsize != d->comm->grank) {
        continue;
      }
      int i0    = it * OVERLAP_TILE;
      int j0    = jt * OVERLAP_TILE;
      int ni    = MIN(OVERLAP_TILE, norb - i0);
      int nj    = MIN(OVERLAP_TILE, norb - j0);
      double *s = d->overlap + i0 + (size_t)j0 * norb;
      if (it == jt) {
        cblas_dsyrk(CblasColMajor,
            CblasUpper,
            CblasTrans,
            ni,
            (int)m,
            d->grid.voxel,
            es->psit + (size_t)i0 * es->ldt,
            (int)es->ldt,
            beta,
            s,
            norb);
      } else {
        cblas_dgemm(CblasColMajor,
            CblasTrans,
            CblasNoTrans,
            ni,
            nj,
            (int)m,
            d->grid.voxel,
            es->psit + (size_t)i0 * es->ldt,
            (int)es->ldt,
            es->psit + (size_t)j0 * es->ldt,
            (int)es->ldt,
            beta,
            s,
            norb);
      }
    }
  }
  commSyncGroup();
}

/* Sums the upper triangle of S over the groups. Lane l < lmin reduces the columns
 * holding the l-th share of the triangle, so the lanes work in parallel. */
static void reduceOverlap(DftType *d)
{
  const CommType *c = d->comm;
  int norb          = d->es.norb;

  if (c->grank < c->lmin) {
    double total = 0.5 * norb * (norb + 1.0);
    int j0 = 0, j1 = 0;
    while (j0 < norb && 0.5 * j0 * (j0 + 1.0) < total * c->grank / c->lmin) {
      j0++;
    }
    j1 = j0;
    while (j1 < norb && 0.5 * j1 * (j1 + 1.0) < total * (c->grank + 1) / c->lmin) {
      j1++;
    }
    size_t pos = 0;
    for (int j = j0; j < j1; j++) {
      memcpy(d->segbuf + pos, d->overlap + (size_t)j * norb, (j + 1) * sizeof(double));
      pos += j + 1;
    }
    commAllreduce(d->segbuf, pos, COMM_SUM, COMM_LANE);
    pos = 0;
    for (int j = j0; j < j1; j++) {
      memcpy(d->overlap + (size_t)j * norb, d->segbuf + pos, (j + 1) * sizeof(double));
      pos += j + 1;
    }
    // nominal volume of a ring allreduce
    addBytes(
        d, PH_REDUCE, 2.0 * (c->ngroups - 1) / c->ngroups * (double)pos * sizeof(double));
  }
  commSyncGroup();
}

/* Overlap phase incl. the transposes and the reduction: the upper triangle of S of
 * the current orbitals in every group's window */
static void overlapPass(DftType *d)
{
  double t0;
  if (d->mode == MODE_STREAMED) {
    commBarrier(COMM_WORLD);
    for (int k = 0; k < d->nchunks; k++) {
      t0          = phaseResume(d, PH_TRANSPOSE);
      size_t rows = transposeChunk(d, k, PH_TRANSPOSE);
      phaseStop(d, PH_TRANSPOSE, t0);
      t0 = phaseResume(d, PH_OVERLAP);
      overlapTiles(d, rows, k == 0 ? 0.0 : 1.0);
      phaseStop(d, PH_OVERLAP, t0);
    }
  } else {
    size_t rows = d->grid.ngtot;
    if (d->mode == MODE_STORED) {
      t0 = phaseStart(d, PH_TRANSPOSE);
      for (int k = 0; k < d->nchunks; k++) {
        transposeChunk(d, k, PH_TRANSPOSE);
      }
      phaseStop(d, PH_TRANSPOSE, t0);
      rows = d->dist.gnrows[d->comm->group];
    }
    t0 = phaseStart(d, PH_OVERLAP);
    overlapTiles(d, rows, 0.0);
    phaseStop(d, PH_OVERLAP, t0);
  }
  if (d->comm->ngroups > 1) {
    t0 = phaseStart(d, PH_REDUCE);
    reduceOverlap(d);
    phaseStop(d, PH_REDUCE, t0);
  }
}

/* Psi <- Psi * X on nrows rows of a, done in row blocks so that only a block-sized
 * buffer is needed */
static void rotateRows(DftType *d, double *a, size_t lda, size_t nrows)
{
  int norb = d->es.norb;
  for (size_t r0 = 0; r0 < nrows; r0 += d->wrows) {
    int m = (int)MIN(d->wrows, nrows - r0);
    cblas_dgemm(CblasColMajor,
        CblasNoTrans,
        CblasNoTrans,
        m,
        norb,
        norb,
        1.0,
        a + r0,
        (int)lda,
        d->xmat,
        norb,
        0.0,
        d->wa,
        m);
    for (int iorb = 0; iorb < norb; iorb++) {
      memcpy(a + iorb * lda + r0, d->wa + (size_t)iorb * m, m * sizeof(double));
    }
  }
}

/* Rotate phase incl. the transposes: every rank rotates its grid rows */
static void rotatePass(DftType *d)
{
  const DistType *ds = &d->dist;
  int me             = d->comm->rank;
  double t0;
  size_t s, n;

  if (d->mode == MODE_STREAMED) {
    // rotating chunk k reads only the old rows of chunk k: the update is in place
    commBarrier(COMM_WORLD);
    for (int k = 0; k < d->nchunks; k++) {
      t0 = phaseResume(d, PH_TRANSPOSE);
      transposeChunk(d, k, PH_TRANSPOSE);
      phaseStop(d, PH_TRANSPOSE, t0);
      t0 = phaseResume(d, PH_ROTATE);
      chunkRows(d, me, k, &s, &n);
      rotateRows(d, d->es.psit + d->choff[d->comm->grank], d->es.ldt, n);
      commSyncGroup();
      phaseStop(d, PH_ROTATE, t0);
      t0 = phaseResume(d, PH_UNTRANSPOSE);
      untransposeChunk(d, k, PH_UNTRANSPOSE);
      phaseStop(d, PH_UNTRANSPOSE, t0);
    }
    return;
  }

  t0 = phaseStart(d, PH_ROTATE);
  rotateRows(
      d, d->es.psit + ds->row0[me] - ds->grow0[d->comm->group], d->es.ldt, ds->nrows[me]);
  commSyncGroup();
  phaseStop(d, PH_ROTATE, t0);
  if (d->mode == MODE_STORED) {
    t0 = phaseStart(d, PH_UNTRANSPOSE);
    for (int k = 0; k < d->nchunks; k++) {
      untransposeChunk(d, k, PH_UNTRANSPOSE);
    }
    phaseStop(d, PH_UNTRANSPOSE, t0);
  }
}

static double maxDeviationFromIdentity(int norb, const double *overlap)
{
  double dev = 0.0;
  for (int j = 0; j < norb; j++) {
    for (int i = 0; i <= j; i++) {
      double tt = overlap[i + (size_t)j * norb] - (i == j ? 1.0 : 0.0);
      dev       = fmax(dev, fabs(tt));
    }
  }
  return dev;
}

static void printMatrix(const char *title, int norb, const double *overlap)
{
  commPrintf("%s\n", title);
  for (int i = 0; i < norb; i++) {
    for (int j = 0; j < norb; j++) {
      double tt =
          (i <= j) ? overlap[i + (size_t)j * norb] : overlap[j + (size_t)i * norb];
      commPrintf("%6.2f", fabs(tt) < 1.E-12 ? 0.0 : tt);
    }
    commPrintf("\n");
  }
}

static bool checkValue(const char *label, double value, double tol)
{
  bool ok = value <= tol;
  commPrintf("%-38s %12.4e  (tol %.0e)  %s\n", label, value, tol, ok ? "PASS" : "FAIL");
  return ok;
}

/* Checks of the freshly placed (non-orthogonal) orbitals */
static bool validateOrbitals(
    const DftType *d, const KineticWorkType *kw, const double *rat)
{
  const ElectronStateType *es = &d->es;
  const GridType *grid        = &d->grid;
  int o0                      = d->dist.orb0[d->comm->rank];
  int n                       = d->dist.norb[d->comm->rank];
  const double *orbitals      = (n > 0) ? rankOrbitals(d) : NULL;
  double res[3]               = { 0.0, 0.0, 0.0 }; // residual, max |psi|, norm error
  double ekinAnalytical       = 0.0;
  bool ok                     = true;

  for (int i = 0; i < n; i++) {
    int iorb          = o0 + i;
    const double *psi = orbitals + i * grid->ld;
    double tt         = testPutGtoSymOrthoOrb(
        &rat[3 * es->iat[iorb]], es->gw[iorb], grid->origin, grid->ng, grid->h, psi);
    double anorm = 0.0;
    double pmax  = 0.0;
    for (size_t j = 0; j < grid->ngtot; j++) {
      anorm += psi[j] * psi[j];
      pmax = fmax(pmax, fabs(psi[j]));
    }
    res[0] = fmax(res[0], tt);
    res[1] = fmax(res[1], pmax);
    res[2] = fmax(res[2], fabs(1.0 - anorm * grid->voxel));
  }
  commAllreduce(res, 3, COMM_MAX, COMM_WORLD);
  for (int iorb = 0; iorb < es->norb; iorb++) {
    // <T> = 3/(2 gw^2) for a normalized exp(-r^2/gw^2) orbital
    ekinAnalytical += 1.5 / (es->gw[iorb] * es->gw[iorb]);
  }
  double ekin = kineticEnergy(kw, d);

  ok &= checkValue(
      "GTO vs. reference (max rel. residual)", res[0] / res[1], TOL_GTO_RESIDUAL);
  ok &= checkValue("GTO normalization (max |1-<psi|psi>|)", res[2], TOL_NORM);
  commPrintf("Ekin numerical %20.10f  analytical %20.10f\n", ekin, ekinAnalytical);
  ok &= checkValue(
      "Ekin (rel. error)", fabs(ekin - ekinAnalytical) / ekinAnalytical, TOL_EKIN);
  return ok;
}

/* Compares a distributed eigensolver with the LAPACK path on world rank 0. sref is a
 * full copy of S taken before the solve and is overwritten. */
static bool checkEigenReference(
    int norb, double *sref, const double *eval, const double *xmat)
{
  bool ok = true;
  if (commRank() == 0) {
    size_t nn      = (size_t)norb * norb;
    double *evref  = allocateArray(norb * sizeof(double));
    double *xref   = allocateArray(nn * sizeof(double));
    double *wa     = allocateArray(nn * sizeof(double));
    double devEval = 0.0, devX = 0.0;
    loewdinLapack(norb, sref, evref, wa, xref);
    for (int i = 0; i < norb; i++) {
      devEval = fmax(devEval, fabs(eval[i] - evref[i]) / fabs(evref[i]));
    }
    for (size_t i = 0; i < nn; i++) {
      devX = fmax(devX, fabs(xmat[i] - xref[i]));
    }
    ok &= checkValue("Eigenvalues vs. LAPACK (max rel. dev.)", devEval, TOL_EIGEN_REF);
    ok &= checkValue("X vs. LAPACK (max abs. dev.)", devX, TOL_EIGEN_REF);
    free(evref);
    free(xref);
    free(wa);
  }
  return ok;
}

static void printMemory(const char *label, double bytes)
{
  if (bytes > GB) {
    commPrintf("%-14s: %.2f GB\n", label, bytes / GB);
  } else {
    commPrintf("%-14s: %.2f MB\n", label, bytes / (1024.0 * 1024.0));
  }
}

int runDft(const PoscarFileType *pf, const ConfigType *cfg, const CommType *comm)
{
  DftType d             = { .comm = comm, .timed = true };
  GridType *grid        = &d.grid;
  ElectronStateType *es = &d.es;
  KineticWorkType kw;
  double gw1, gw2, rgcut;
  bool pass = true;
  bool mpi  = commIsMpi();

  gw1       = 1.11 / BOHR2ANG * 0.5;
  gw2       = 1.11 / BOHR2ANG * 0.6;
  rgcut     = 6.0 * fmax(gw1, gw2);

  setupGrid(pf, cfg, grid);

  es->norb = 2 * pf->nat;
  int norb = es->norb;
  setupDistribution(&d);
  double budget = 0.0, groupMem = 0.0;
  d.mode    = MODE_NONE;
  d.nchunks = 1;
  if (mpi) {
    planMemory(&d, cfg, &budget, &groupMem);
  }

  size_t nk      = (size_t)(grid->ng[0] / 2 + 1) * grid->ng[1] * grid->ng[2];
  double memPsi  = (double)grid->ld * norb * sizeof(double);
  double memWork = (double)(grid->ngtot + 2 * nk) * sizeof(double) +
                   (3.0 * norb * norb + (double)ROTATE_BLOCK * norb) * sizeof(double);

  commPrintf("Input         : %s (%d atoms)\n", cfg->poscarFile, pf->nat);
  commPrintf("Cell (bohr)   : %10.5f %10.5f %10.5f\n",
      pf->cell[0][0],
      pf->cell[1][1],
      pf->cell[2][2]);
  commPrintf("Grid          : %d x %d x %d = %zu points\n",
      grid->ng[0],
      grid->ng[1],
      grid->ng[2],
      grid->ngtot);
  commPrintf(
      "hgrid (bohr)  : %10.5f %10.5f %10.5f\n", grid->h[0], grid->h[1], grid->h[2]);
  commPrintf("Orbitals      : %d\n", norb);
  printMemory("Orbitals mem", memPsi);
  if (mpi) {
    int gmin = comm->size, gmax = 0;
    for (int g = 0; g < comm->ngroups; g++) {
      gmin = MIN(gmin, comm->groupSize[g]);
      gmax = MAX(gmax, comm->groupSize[g]);
    }
    commPrintf("Ranks         : %d\n", comm->size);
    if (gmin == gmax) {
      commPrintf("Groups        : %d shared-memory group%s with %d ranks%s\n",
          comm->ngroups,
          comm->ngroups > 1 ? "s" : "",
          gmin,
          comm->ngroups > 1 ? " each" : "");
    } else {
      commPrintf("Groups        : %d shared-memory groups with %d to %d ranks\n",
          comm->ngroups,
          gmin,
          gmax);
    }
    if (d.mode == MODE_NONE) {
      commPrintf("Transpose     : %s\n", MODE_NAMES[d.mode]);
    } else {
      commPrintf("Transpose     : %s, %d chunk%s\n",
          MODE_NAMES[d.mode],
          d.nchunks,
          d.nchunks > 1 ? "s" : "");
    }
    printMemory("Group budget", budget);
    printMemory("Group mem", groupMem);
    memWork = (double)(grid->ngtot + 2 * nk) * sizeof(double) +
              (double)MAX((size_t)1, MIN((size_t)ROTATE_BLOCK, d.dist.nrows[0])) * norb *
                  sizeof(double) +
              (double)norb * norb * sizeof(double);
    printMemory("Work mem/rank", memWork);
  } else {
    printMemory("Work mem", memWork);
  }
  eigenInit(&d.eig, comm, cfg, norb);
  eigenPrintInfo(&d.eig);
  if (d.dist.idle > 0) {
    commPrintf("NOTICE: %d of %d ranks own no orbitals and idle in the gto and kinetic "
               "phases\n",
        d.dist.idle,
        comm->size);
  }
  commPrintf(HLINE);
  fflush(stdout);

  const int *phases = mpi ? MPI_PHASES : SEQUENTIAL_PHASES;
  int nphases       = mpi ? NPHASES : (int)(sizeof(SEQUENTIAL_PHASES) / sizeof(int));
  for (int i = 0; i < nphases; i++) {
    LIKWID_MARKER_REGISTER(PHASE_NAMES[phases[i]]);
  }

  es->iat = malloc(norb * sizeof(int));
  es->gw  = malloc(norb * sizeof(double));
  for (int iorb = 0; iorb < norb; iorb++) {
    es->iat[iorb] = iorb / 2;
    es->gw[iorb]  = (iorb % 2 == 0) ? gw1 : gw2;
  }
  allocateState(&d);
  kineticInit(&kw, grid);

  bool reference = cfg->validate && strcmp(eigenBackendName(), "LAPACK") != 0;
  double *sref   = NULL;
  if (reference && norb > EIGEN_REFERENCE_MAX) {
    reference = false;
    commPrintf("LAPACK reference check of the eigensolver skipped (%d > %d orbitals)\n",
        norb,
        EIGEN_REFERENCE_MAX);
  }

  double tmin[NPHASES], tsum[NPHASES];
  double ekin = 0.0;
  for (int ph = 0; ph < NPHASES; ph++) {
    tmin[ph] = 1.0E30;
    tsum[ph] = 0.0;
  }

  for (int rep = 0; rep < cfg->reps; rep++) {
    bool check = cfg->validate && rep == 0;
    double t0;

    for (int ph = 0; ph < NPHASES; ph++) {
      d.t[ph]     = 0.0;
      d.bytes[ph] = 0.0;
    }

    t0 = phaseStart(&d, PH_GTO);
    placeOrbitals(&d, pf->rat, rgcut);
    commSyncGroup();
    phaseStop(&d, PH_GTO, t0);

    if (check) {
      commPrintf("Validation\n");
      pass &= validateOrbitals(&d, &kw, pf->rat);
    }

    overlapPass(&d);

    if (cfg->verbose && rep == 0 && norb <= 64) {
      printMatrix("Overlap matrix:", norb, d.overlap);
    }
    if (check && reference && comm->rank == 0) {
      size_t nn = (size_t)norb * norb;
      sref      = allocateArray(nn * sizeof(double));
      memcpy(sref, d.overlap, nn * sizeof(double));
    }

    t0 = phaseStart(&d, PH_EIGEN);
    eigenSolve(&d.eig, d.overlap, d.eval, d.xmat);
    phaseStop(&d, PH_EIGEN, t0);

    t0 = phaseStart(&d, PH_GATHER);
    addBytes(&d, PH_GATHER, eigenGather(&d.eig, d.eval, d.xmat));
    phaseStop(&d, PH_GATHER, t0);

    if (check && reference) {
      pass &= checkEigenReference(norb, sref, d.eval, d.xmat);
      free(sref);
      sref = NULL;
    }

    rotatePass(&d);

    t0   = phaseStart(&d, PH_KINETIC);
    ekin = kineticEnergy(&kw, &d);
    phaseStop(&d, PH_KINETIC, t0);

    if (check) {
      d.timed = false;
      overlapPass(&d);
      d.timed = true;
      if (cfg->verbose && norb <= 64) {
        printMatrix("Overlap matrix after orthogonalization:", norb, d.overlap);
      }
      pass &= checkValue("Orthonormality (max |S-I|)",
          maxDeviationFromIdentity(norb, d.overlap),
          TOL_ORTHO);
      commPrintf(HLINE);
    }

    // a phase takes as long as its slowest rank
    commReduceRoot(d.t, NPHASES, COMM_MAX);
    for (int ph = 0; ph < NPHASES; ph++) {
      tmin[ph] = fmin(tmin[ph], d.t[ph]);
      tsum[ph] += d.t[ph];
    }
  }
  commReduceRoot(d.bytes, NPHASES, COMM_SUM);

  commPrintf(
      "Overlap eigenvalues: min %12.6e  max %12.6e\n", d.eval[0], d.eval[norb - 1]);
  if (cfg->verbose) {
    for (int iorb = 0; iorb < norb; iorb++) {
      commPrintf("%5d%14.6e\n", iorb, d.eval[iorb]);
    }
  }
  commPrintf("Ekin (orthonormalized orbitals) = %20.10f\n", ekin);
  commPrintf(HLINE);

  // nominal flop counts per repetition (FFT: 2.5 N log2 N per real transform)
  double n              = (double)grid->ngtot;
  double flops[NPHASES] = { 0.0 };
  flops[PH_OVERLAP]     = (double)norb * (norb + 1) * n;
  flops[PH_ROTATE]      = 2.0 * n * norb * norb;
  flops[PH_KINETIC]     = norb * (5.0 * n * log2(n) + 4.0 * n);

  double ttmin = 0.0, ttavg = 0.0;
  if (mpi) {
    commPrintf("%-11s %14s %14s %14s %14s\n",
        "Phase",
        "min time [s]",
        "avg time [s]",
        "GFLOP/s",
        "GB/s");
  } else {
    commPrintf(
        "%-10s %14s %14s %14s\n", "Phase", "min time [s]", "avg time [s]", "GFLOP/s");
  }
  for (int i = 0; i < nphases; i++) {
    int ph      = phases[i];
    double tavg = tsum[ph] / cfg->reps;
    ttmin += tmin[ph];
    ttavg += tavg;
    char rate[2][32] = { "-", "-" };
    if (flops[ph] > 0.0 && tmin[ph] > 0.0) {
      snprintf(rate[0], sizeof(rate[0]), "%.2f", flops[ph] / tmin[ph] * 1.0E-9);
    }
    if (d.bytes[ph] > 0.0 && tmin[ph] > 0.0) {
      snprintf(rate[1], sizeof(rate[1]), "%.2f", d.bytes[ph] / tmin[ph] * 1.0E-9);
    }
    if (mpi) {
      commPrintf("%-11s %14.6f %14.6f %14s %14s\n",
          PHASE_NAMES[ph],
          tmin[ph],
          tavg,
          rate[0],
          rate[1]);
    } else {
      commPrintf("%-10s %14.6f %14.6f %14s\n", PHASE_NAMES[ph], tmin[ph], tavg, rate[0]);
    }
  }
  commPrintf(
      mpi ? "%-11s %14.6f %14.6f\n" : "%-10s %14.6f %14.6f\n", "total", ttmin, ttavg);
  commPrintf(HLINE);
  pass = commAllAnd(pass);
  if (cfg->validate) {
    commPrintf("Validation %s\n", pass ? "PASSED" : "FAILED");
  }

  kineticFree(&kw);
  eigenFree(&d.eig);
  freeState(&d);
  freeDistribution(&d.dist);
  free(es->iat);
  free(es->gw);

  return pass ? EXIT_SUCCESS : EXIT_FAILURE;
}
