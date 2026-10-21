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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "allocate.h"
#include "blasLapack.h"
#include "constants.h"
#include "dft.h"
#include "gtoOnGrid.h"
#include "likwid-marker.h"
#include "poscarReader.h"
#include "timing.h"
#include "util.h"

enum { PH_GTO, PH_OVERLAP, PH_EIGEN, PH_ROTATE, PH_KINETIC, NPHASES };
static const char *const PHASE_NAMES[NPHASES] = {
  "gto", "overlap", "eigen", "rotate", "kinetic"
};

typedef struct {
  fftw_plan planf, planb;
  int nthreads;
  double **rbuf; /* per-thread T|psi> in real space */
  fftw_complex **cbuf; /* per-thread psi in reciprocal space */
  double *k2[3]; /* |k|^2 components along each direction */
} KineticWorkType;

static void *allocateArray(size_t bytes)
{
  return allocate(ARRAY_ALIGNMENT, bytes > 0 ? bytes : ARRAY_ALIGNMENT);
}

static int getNumThreads(void)
{
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}

static int getThreadNum(void)
{
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

static double phaseStart(int phase)
{
#ifdef LIKWID_PERFMON
#pragma omp parallel
  {
    LIKWID_MARKER_START(PHASE_NAMES[phase]);
  }
#else
  (void)phase;
#endif
  return getTimeStamp();
}

static double phaseStop(int phase, double t0)
{
  double t = getTimeStamp() - t0;
#ifdef LIKWID_PERFMON
#pragma omp parallel
  {
    LIKWID_MARKER_STOP(PHASE_NAMES[phase]);
  }
#else
  (void)phase;
#endif
  return t;
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

static void setupGrid(const PoscarFileType *pf, const ConfigType *cfg, GridType *grid)
{
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i != j && fabs(pf->cell[i][j]) > 1.0E-10 * fabs(pf->cell[i][i])) {
        fprintf(stderr, "ERROR: only orthogonal cells are supported\n");
        exit(EXIT_FAILURE);
      }
    }
    if (pf->cell[i][i] <= 0.0) {
      fprintf(stderr, "ERROR: cell vector %d must point in positive direction\n", i + 1);
      exit(EXIT_FAILURE);
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
    fprintf(stderr,
        "ERROR: %zu grid points exceed the 32 bit integer BLAS interface\n",
        grid->ngtot);
    exit(EXIT_FAILURE);
  }
}

static void placeOrbitals(
    ElectronStateType *es, const GridType *grid, const double *rat, double rgcut)
{
#pragma omp parallel for schedule(dynamic)
  for (int iorb = 0; iorb < es->norb; iorb++) {
    putGtoSymOrthoOrb(&rat[3 * es->iat[iorb]],
        es->gw[iorb],
        rgcut,
        grid->origin,
        grid->ng,
        grid->h,
        es->psi + iorb * grid->ld);
  }
}

static void kineticInit(KineticWorkType *kw, const GridType *grid, double *psi)
{
  const int *ng = grid->ng;
  size_t nk     = (size_t)(ng[0] / 2 + 1) * ng[1] * ng[2];

  kw->nthreads  = getNumThreads();
  kw->rbuf      = malloc(kw->nthreads * sizeof(double *));
  kw->cbuf      = malloc(kw->nthreads * sizeof(fftw_complex *));
  for (int t = 0; t < kw->nthreads; t++) {
    kw->rbuf[t] = allocateArray(grid->ngtot * sizeof(double));
    kw->cbuf[t] = allocateArray(nk * sizeof(fftw_complex));
  }

  // x is stored contiguously, so it is the last (halved) dimension for FFTW
  kw->planf = fftw_plan_dft_r2c_3d(ng[2], ng[1], ng[0], psi, kw->cbuf[0], FFTW_ESTIMATE);
  kw->planb =
      fftw_plan_dft_c2r_3d(ng[2], ng[1], ng[0], kw->cbuf[0], kw->rbuf[0], FFTW_ESTIMATE);

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
  for (int t = 0; t < kw->nthreads; t++) {
    free(kw->rbuf[t]);
    free(kw->cbuf[t]);
  }
  free(kw->rbuf);
  free(kw->cbuf);
  for (int d = 0; d < 3; d++) {
    free(kw->k2[d]);
  }
}

/* <psi_i| -1/2 nabla^2 |psi_i> summed over all orbitals, using one FFT pair per orbital */
static double kineticEnergy(
    const KineticWorkType *kw, const ElectronStateType *es, const GridType *grid)
{
  const int *ng = grid->ng;
  int nkx       = ng[0] / 2 + 1;
  double scale  = 0.5 / (double)grid->ngtot; // 1/2 from T, 1/N from the c2r transform
  double ekin   = 0.0;

#pragma omp parallel
  {
    double *tpsi       = kw->rbuf[getThreadNum()];
    fftw_complex *cpsi = kw->cbuf[getThreadNum()];

#pragma omp for schedule(dynamic) reduction(+ : ekin)
    for (int iorb = 0; iorb < es->norb; iorb++) {
      double *psi = es->psi + iorb * grid->ld;
      fftw_execute_dft_r2c(kw->planf, psi, cpsi);
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
  }
  return ekin;
}

/* Upper triangle of S = voxel * Psi^T Psi (column-major, norb x norb) */
static void computeOverlap(
    const ElectronStateType *es, const GridType *grid, double *overlap)
{
  cblas_dsyrk(CblasColMajor,
      CblasUpper,
      CblasTrans,
      es->norb,
      (int)grid->ngtot,
      grid->voxel,
      es->psi,
      (int)grid->ld,
      0.0,
      overlap,
      es->norb);
}

/* Computes X = S^{-1/2} = U diag(eval^{-1/2}) U^T. The overlap is overwritten by U. */
static void loewdinMatrix(int norb, double *overlap, double *eval, double *wa, double *x)
{
  lapack_int info = LAPACKE_dsyevd(LAPACK_COL_MAJOR, 'V', 'U', norb, overlap, norb, eval);
  if (info != 0) {
    fprintf(stderr,
        "ERROR: info=%d in dsyevd used for Loewdin orthogonalization\n",
        (int)info);
    exit(EXIT_FAILURE);
  }
  if (eval[0] <= 0.0) {
    fprintf(stderr,
        "ERROR: overlap matrix is not positive definite (smallest eigenvalue %e)\n",
        eval[0]);
    exit(EXIT_FAILURE);
  }
  for (int j = 0; j < norb; j++) {
    double fac = 1.0 / sqrt(eval[j]);
    for (int i = 0; i < norb; i++) {
      wa[i + (size_t)j * norb] = overlap[i + (size_t)j * norb] * fac;
    }
  }
  cblas_dgemm(CblasColMajor,
      CblasNoTrans,
      CblasTrans,
      norb,
      norb,
      norb,
      1.0,
      wa,
      norb,
      overlap,
      norb,
      0.0,
      x,
      norb);
}

/* Psi <- Psi * X, done in row blocks so that only a block-sized buffer is needed */
static void rotateOrbitals(
    ElectronStateType *es, const GridType *grid, const double *x, double *wa)
{
  int norb = es->norb;
  for (size_t r0 = 0; r0 < grid->ngtot; r0 += ROTATE_BLOCK) {
    int m = (int)MIN((size_t)ROTATE_BLOCK, grid->ngtot - r0);
    cblas_dgemm(CblasColMajor,
        CblasNoTrans,
        CblasNoTrans,
        m,
        norb,
        norb,
        1.0,
        es->psi + r0,
        (int)grid->ld,
        x,
        norb,
        0.0,
        wa,
        m);
#pragma omp parallel for schedule(static)
    for (int iorb = 0; iorb < norb; iorb++) {
      memcpy(es->psi + iorb * grid->ld + r0, wa + (size_t)iorb * m, m * sizeof(double));
    }
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
  printf("%s\n", title);
  for (int i = 0; i < norb; i++) {
    for (int j = 0; j < norb; j++) {
      double tt =
          (i <= j) ? overlap[i + (size_t)j * norb] : overlap[j + (size_t)i * norb];
      printf("%6.2f", fabs(tt) < 1.E-12 ? 0.0 : tt);
    }
    printf("\n");
  }
}

static bool checkValue(const char *label, double value, double tol)
{
  bool ok = value <= tol;
  printf("%-38s %12.4e  (tol %.0e)  %s\n", label, value, tol, ok ? "PASS" : "FAIL");
  return ok;
}

/* Checks of the freshly placed (non-orthogonal) orbitals */
static bool validateOrbitals(const ElectronStateType *es,
    const GridType *grid,
    const KineticWorkType *kw,
    const double *rat)
{
  double res = 0.0, amax = 0.0, normErr = 0.0;
  double ekinAnalytical = 0.0;
  bool ok               = true;

#pragma omp parallel for schedule(dynamic) reduction(max : res, amax, normErr)
  for (int iorb = 0; iorb < es->norb; iorb++) {
    const double *psi = es->psi + iorb * grid->ld;
    double tt         = testPutGtoSymOrthoOrb(
        &rat[3 * es->iat[iorb]], es->gw[iorb], grid->origin, grid->ng, grid->h, psi);
    double anorm = 0.0;
    double pmax  = 0.0;
    for (size_t i = 0; i < grid->ngtot; i++) {
      anorm += psi[i] * psi[i];
      pmax = fmax(pmax, fabs(psi[i]));
    }
    res     = fmax(res, tt);
    amax    = fmax(amax, pmax);
    normErr = fmax(normErr, fabs(1.0 - anorm * grid->voxel));
  }
  for (int iorb = 0; iorb < es->norb; iorb++) {
    // <T> = 3/(2 gw^2) for a normalized exp(-r^2/gw^2) orbital
    ekinAnalytical += 1.5 / (es->gw[iorb] * es->gw[iorb]);
  }
  double ekin = kineticEnergy(kw, es, grid);

  ok &= checkValue("GTO vs. reference (max rel. residual)", res / amax, TOL_GTO_RESIDUAL);
  ok &= checkValue("GTO normalization (max |1-<psi|psi>|)", normErr, TOL_NORM);
  printf("Ekin numerical %20.10f  analytical %20.10f\n", ekin, ekinAnalytical);
  ok &= checkValue(
      "Ekin (rel. error)", fabs(ekin - ekinAnalytical) / ekinAnalytical, TOL_EKIN);
  return ok;
}

static void printMemory(const char *label, double bytes)
{
  if (bytes > 1024.0 * 1024.0 * 1024.0) {
    printf("%-14s: %.2f GB\n", label, bytes / (1024.0 * 1024.0 * 1024.0));
  } else {
    printf("%-14s: %.2f MB\n", label, bytes / (1024.0 * 1024.0));
  }
}

int runDft(const PoscarFileType *pf, const ConfigType *cfg)
{
  GridType grid;
  ElectronStateType es;
  KineticWorkType kw;
  double gw1, gw2, rgcut;
  bool pass = true;

  gw1       = 1.11 / BOHR2ANG * 0.5;
  gw2       = 1.11 / BOHR2ANG * 0.6;
  rgcut     = 6.0 * fmax(gw1, gw2);

  setupGrid(pf, cfg, &grid);

  es.norb        = 2 * pf->nat;
  int norb       = es.norb;
  size_t nk      = (size_t)(grid.ng[0] / 2 + 1) * grid.ng[1] * grid.ng[2];
  double memPsi  = (double)grid.ld * norb * sizeof(double);
  double memWork = (double)getNumThreads() * (grid.ngtot + 2 * nk) * sizeof(double) +
                   (3.0 * norb * norb + (double)ROTATE_BLOCK * norb) * sizeof(double);

  printf("Input         : %s (%d atoms)\n", cfg->poscarFile, pf->nat);
  printf("Cell (bohr)   : %10.5f %10.5f %10.5f\n",
      pf->cell[0][0],
      pf->cell[1][1],
      pf->cell[2][2]);
  printf("Grid          : %d x %d x %d = %zu points\n",
      grid.ng[0],
      grid.ng[1],
      grid.ng[2],
      grid.ngtot);
  printf("hgrid (bohr)  : %10.5f %10.5f %10.5f\n", grid.h[0], grid.h[1], grid.h[2]);
  printf("Orbitals      : %d\n", norb);
  printMemory("Orbitals mem", memPsi);
  printMemory("Work mem", memWork);
  printf(HLINE);
  fflush(stdout);

  es.iat = malloc(norb * sizeof(int));
  es.gw  = malloc(norb * sizeof(double));
  es.psi = allocateArray((size_t)grid.ld * norb * sizeof(double));
  for (int iorb = 0; iorb < norb; iorb++) {
    es.iat[iorb] = iorb / 2;
    es.gw[iorb]  = (iorb % 2 == 0) ? gw1 : gw2;
  }

  // plan before the orbitals are filled, planning may overwrite the arrays
  kineticInit(&kw, &grid, es.psi);

  double *overlap = allocateArray((size_t)norb * norb * sizeof(double));
  double *xmat    = allocateArray((size_t)norb * norb * sizeof(double));
  double *wa      = allocateArray(
      (size_t)MAX((size_t)norb, (size_t)ROTATE_BLOCK) * norb * sizeof(double));
  double *eval = allocateArray(norb * sizeof(double));

  double tmin[NPHASES], tsum[NPHASES];
  double ekin = 0.0;
  for (int ph = 0; ph < NPHASES; ph++) {
    tmin[ph] = 1.0E30;
    tsum[ph] = 0.0;
  }

  for (int rep = 0; rep < cfg->reps; rep++) {
    bool check = cfg->validate && rep == 0;
    double t[NPHASES], t0;

    t0 = phaseStart(PH_GTO);
    placeOrbitals(&es, &grid, pf->rat, rgcut);
    t[PH_GTO] = phaseStop(PH_GTO, t0);

    if (check) {
      printf("Validation\n");
      pass &= validateOrbitals(&es, &grid, &kw, pf->rat);
    }

    t0 = phaseStart(PH_OVERLAP);
    computeOverlap(&es, &grid, overlap);
    t[PH_OVERLAP] = phaseStop(PH_OVERLAP, t0);

    if (cfg->verbose && rep == 0 && norb <= 64) {
      printMatrix("Overlap matrix:", norb, overlap);
    }

    t0 = phaseStart(PH_EIGEN);
    loewdinMatrix(norb, overlap, eval, wa, xmat);
    t[PH_EIGEN] = phaseStop(PH_EIGEN, t0);

    t0          = phaseStart(PH_ROTATE);
    rotateOrbitals(&es, &grid, xmat, wa);
    t[PH_ROTATE]  = phaseStop(PH_ROTATE, t0);

    t0            = phaseStart(PH_KINETIC);
    ekin          = kineticEnergy(&kw, &es, &grid);
    t[PH_KINETIC] = phaseStop(PH_KINETIC, t0);

    if (check) {
      computeOverlap(&es, &grid, overlap);
      if (cfg->verbose && norb <= 64) {
        printMatrix("Overlap matrix after orthogonalization:", norb, overlap);
      }
      pass &= checkValue("Orthonormality (max |S-I|)",
          maxDeviationFromIdentity(norb, overlap),
          TOL_ORTHO);
      printf(HLINE);
    }

    for (int ph = 0; ph < NPHASES; ph++) {
      tmin[ph] = fmin(tmin[ph], t[ph]);
      tsum[ph] += t[ph];
    }
  }

  printf("Overlap eigenvalues: min %12.6e  max %12.6e\n", eval[0], eval[norb - 1]);
  if (cfg->verbose) {
    for (int iorb = 0; iorb < norb; iorb++) {
      printf("%5d%14.6e\n", iorb, eval[iorb]);
    }
  }
  printf("Ekin (orthonormalized orbitals) = %20.10f\n", ekin);
  printf(HLINE);

  // nominal flop counts per repetition (FFT: 2.5 N log2 N per real transform)
  double n              = (double)grid.ngtot;
  double flops[NPHASES] = { 0.0 };
  flops[PH_OVERLAP]     = (double)norb * (norb + 1) * n;
  flops[PH_ROTATE]      = 2.0 * n * norb * norb;
  flops[PH_KINETIC]     = norb * (5.0 * n * log2(n) + 4.0 * n);

  double ttmin = 0.0, ttavg = 0.0;
  printf("%-10s %14s %14s %14s\n", "Phase", "min time [s]", "avg time [s]", "GFLOP/s");
  for (int ph = 0; ph < NPHASES; ph++) {
    double tavg = tsum[ph] / cfg->reps;
    ttmin += tmin[ph];
    ttavg += tavg;
    if (flops[ph] > 0.0) {
      printf("%-10s %14.6f %14.6f %14.2f\n",
          PHASE_NAMES[ph],
          tmin[ph],
          tavg,
          flops[ph] / tmin[ph] * 1.0E-9);
    } else {
      printf("%-10s %14.6f %14.6f %14s\n", PHASE_NAMES[ph], tmin[ph], tavg, "-");
    }
  }
  printf("%-10s %14.6f %14.6f\n", "total", ttmin, ttavg);
  printf(HLINE);
  if (cfg->validate) {
    printf("Validation %s\n", pass ? "PASSED" : "FAILED");
  }

  kineticFree(&kw);
  free(overlap);
  free(xmat);
  free(wa);
  free(eval);
  free(es.iat);
  free(es.gw);
  free(es.psi);

  return pass ? EXIT_SUCCESS : EXIT_FAILURE;
}
