/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of DFT-Bench.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "allocate.h"
#include "comm.h"
#include "constants.h"
#include "eigen.h"
#include "util.h"

/* ScaLAPACK backend, also used by the ELPA backend (EIGEN_ELPA) for everything but the
 * eigensolver call. The grid ranks copy their block-cyclic blocks of S from the
 * group-shared S, solve for U and the eigenvalues, form W = U diag(eval^-1/4) and
 * X = W W^T with pdsyrk. The gather phase redistributes X to world rank 0 with
 * pdgemr2d and broadcasts it over the group leaders. */

/* BLACS and ScaLAPACK interfaces of netlib and MKL */
// NOLINTBEGIN(readability-identifier-naming)
void Cblacs_get(int icontxt, int what, int *val);
void Cblacs_gridmap(int *icontxt, int *usermap, int ldumap, int nprow, int npcol);
void Cblacs_gridinfo(int icontxt, int *nprow, int *npcol, int *myrow, int *mycol);
void Cblacs_gridexit(int icontxt);
int numroc_(const int *n,
    const int *nb,
    const int *iproc,
    const int *isrcproc,
    const int *nprocs);
void descinit_(int *desc,
    const int *m,
    const int *n,
    const int *mb,
    const int *nb,
    const int *irsrc,
    const int *icsrc,
    const int *ictxt,
    const int *lld,
    int *info);
void pdsyevd_(const char *jobz,
    const char *uplo,
    const int *n,
    double *a,
    const int *ia,
    const int *ja,
    const int *desca,
    double *w,
    double *z,
    const int *iz,
    const int *jz,
    const int *descz,
    double *work,
    const int *lwork,
    int *iwork,
    const int *liwork,
    int *info);
void pdsyrk_(const char *uplo,
    const char *trans,
    const int *n,
    const int *k,
    const double *alpha,
    const double *a,
    const int *ia,
    const int *ja,
    const int *desca,
    const double *beta,
    double *c,
    const int *ic,
    const int *jc,
    const int *descc);
void pdgemr2d_(const int *m,
    const int *n,
    const double *a,
    const int *ia,
    const int *ja,
    const int *desca,
    double *b,
    const int *ib,
    const int *jb,
    const int *descb,
    const int *ictxt);
// NOLINTEND(readability-identifier-naming)

#ifdef EIGEN_ELPA
/* eigenElpa.c: eigenvalues w and eigenvectors z of the distributed matrix a */
void elpaSolve(int n,
    int nb,
    int mloc,
    int nloc,
    int myrow,
    int mycol,
    int fcomm,
    double *a,
    double *w,
    double *z);
#endif

typedef struct {
  int ctxt; /* BLACS context of the eigen grid, -1 on other ranks */
  int ctxt0; /* 1 x 1 context of world rank 0 */
  int myrow, mycol;
  int mloc, nloc, lld; /* local block-cyclic part */
  int desc[9]; /* descriptor of A, Z and X */
  int desc0[9]; /* descriptor of the full X on world rank 0 */
  double *a, *z, *x; /* local parts of S, of the eigenvectors and of X */
  double *w; /* eigenvalues, on every grid rank */
  double *work;
  int *iwork;
  int lwork, liwork;
  int fcomm; /* ELPA: communicator of the grid ranks */
  int gridGroups; /* groups that hold grid ranks */
  int gridPerGroup; /* most grid ranks in one group */
} ScalapackType;

const char *eigenBackendName(void)
{
#ifdef EIGEN_ELPA
  return "ELPA";
#else
  return "SCALAPACK";
#endif
}

/* Global index of local index l in a block-cyclic distribution */
static int globalIndex(int l, int nb, int iproc, int nprocs)
{
  return ((l / nb) * nprocs + iproc) * nb + l % nb;
}

void eigenInit(EigenType *eig, const CommType *comm, const ConfigType *cfg, int norb)
{
  ScalapackType *sp = calloc(1, sizeof(ScalapackType));
  int nproc         = comm->size;
  int nb            = EIGEN_NB;
  int zero = 0, one = 1, info;

  eig->comm   = comm;
  eig->norb   = norb;
  eig->spread = cfg->eigenSpread;
  eig->wa     = NULL;
  eig->priv   = sp;

  if (cfg->eigenGrid[0] > 0) {
    eig->grid[0] = cfg->eigenGrid[0];
    eig->grid[1] = cfg->eigenGrid[1];
    if (eig->grid[0] * eig->grid[1] > nproc) {
      dftAbort("ERROR: the eigen grid %d x %d needs %d ranks, the job has %d\n",
          eig->grid[0],
          eig->grid[1],
          eig->grid[0] * eig->grid[1],
          nproc);
    }
  } else {
    // square grid with about EIGEN_BLOCK_TARGET^2 elements per rank
    int q        = MAX(1, MIN(norb / EIGEN_BLOCK_TARGET, (int)sqrt((double)nproc)));
    eig->grid[0] = eig->grid[1] = q;
  }
  int npr = eig->grid[0], npc = eig->grid[1];
  int ngrid = npr * npc;
#ifdef EIGEN_ELPA
  // ELPA fails if a grid row or column owns no block of the matrix
  if ((norb + nb - 1) / nb < MAX(npr, npc)) {
    dftAbort("ERROR: ELPA needs at least one %d x %d block per grid row and column, "
             "the eigen grid %d x %d is too large for %d orbitals\n",
        nb,
        nb,
        npr,
        npc,
        norb);
  }
#endif

  // grid ranks: packed takes whole groups in group order, spread takes group rank 0
  // of every group, then group rank 1, ... World rank 0 is always the first one.
  int *order = malloc(nproc * sizeof(int)); // world ranks in group order
  int *first = calloc(comm->ngroups + 1, sizeof(int));
  int *ranks = malloc(ngrid * sizeof(int));
  int count  = 0;
  for (int g = 0; g < comm->ngroups; g++) {
    first[g + 1] = first[g] + comm->groupSize[g];
  }
  for (int p = 0; p < nproc; p++) {
    order[first[comm->groupOf[p]] + comm->grankOf[p]] = p;
  }
  if (eig->spread) {
    for (int r = 0; count < ngrid; r++) {
      for (int g = 0; g < comm->ngroups && count < ngrid; g++) {
        if (r < comm->groupSize[g]) {
          ranks[count++] = order[first[g] + r];
        }
      }
    }
  } else {
    memcpy(ranks, order, ngrid * sizeof(int));
  }
  free(order);
  free(first);
  int *perGroup = calloc(comm->ngroups, sizeof(int));
  for (int k = 0; k < ngrid; k++) {
    perGroup[comm->groupOf[ranks[k]]]++;
  }
  for (int g = 0; g < comm->ngroups; g++) {
    sp->gridGroups += (perGroup[g] > 0);
    sp->gridPerGroup = MAX(sp->gridPerGroup, perGroup[g]);
  }
  free(perGroup);
  int *usermap = malloc(ngrid * sizeof(int));
  for (int k = 0; k < ngrid; k++) {
    usermap[(k / npc) + (k % npc) * npr] = ranks[k]; // row-major grid positions
  }
  int nr, nc;
  Cblacs_get(-1, 0, &sp->ctxt);
  Cblacs_gridmap(&sp->ctxt, usermap, npr, npr, npc);
  int zeroRank = 0;
  Cblacs_get(-1, 0, &sp->ctxt0);
  Cblacs_gridmap(&sp->ctxt0, &zeroRank, 1, 1, 1);
  free(usermap);
  free(ranks);

  sp->myrow = sp->mycol = -1;
  bool member           = false;
  for (int k = 0; k < 9; k++) {
    sp->desc[k] = sp->desc0[k] = 0;
  }
  sp->desc[1] = sp->desc0[1] = -1;
  if (sp->ctxt >= 0) {
    Cblacs_gridinfo(sp->ctxt, &nr, &nc, &sp->myrow, &sp->mycol);
    member = sp->myrow >= 0 && sp->mycol >= 0;
  }
  sp->fcomm = commSubsetFortran(member);
  if (!member) {
    sp->ctxt = -1;
    return;
  }

  sp->mloc = numroc_(&norb, &nb, &sp->myrow, &zero, &npr);
  sp->nloc = numroc_(&norb, &nb, &sp->mycol, &zero, &npc);
  sp->lld  = MAX(1, sp->mloc);
  descinit_(sp->desc, &norb, &norb, &nb, &nb, &zero, &zero, &sp->ctxt, &sp->lld, &info);
  if (comm->rank == 0) {
    descinit_(sp->desc0, &norb, &norb, &nb, &nb, &zero, &zero, &sp->ctxt0, &norb, &info);
  }
  size_t loc = (size_t)sp->lld * MAX(1, sp->nloc);
  sp->a      = allocate(ARRAY_ALIGNMENT, loc * sizeof(double));
  sp->z      = allocate(ARRAY_ALIGNMENT, loc * sizeof(double));
  sp->x      = allocate(ARRAY_ALIGNMENT, loc * sizeof(double));
  sp->w      = allocate(ARRAY_ALIGNMENT, norb * sizeof(double));

#ifndef EIGEN_ELPA
  // workspace query
  double lw;
  int liw, query = -1;
  pdsyevd_("V",
      "U",
      &norb,
      sp->a,
      &one,
      &one,
      sp->desc,
      sp->w,
      sp->z,
      &one,
      &one,
      sp->desc,
      &lw,
      &query,
      &liw,
      &query,
      &info);
  // the query of netlib ScaLAPACK omits the workspace of the back transformation
  // (pdormtr), with the queried size the eigenvectors can be silently wrong
  sp->lwork  = (int)lw + 1 + (sp->mloc + sp->nloc + nb) * nb;
  sp->liwork = MAX(liw, 1);
  sp->work   = allocate(ARRAY_ALIGNMENT, sp->lwork * sizeof(double));
  sp->iwork  = malloc(sp->liwork * sizeof(int));
#else
  (void)one;
#endif
}

void eigenPrintInfo(const EigenType *eig)
{
  const ScalapackType *sp = eig->priv;
  commPrintf("Eigensolver   : %s on a %d x %d grid\n",
      eigenBackendName(),
      eig->grid[0],
      eig->grid[1]);
  commPrintf("Eigen ranks   : %s placement, %d group%s with up to %d grid ranks\n",
      eig->spread ? "spread" : "packed",
      sp->gridGroups,
      sp->gridGroups > 1 ? "s" : "",
      sp->gridPerGroup);
}

void eigenFree(EigenType *eig)
{
  ScalapackType *sp = eig->priv;
  if (sp->ctxt >= 0) {
    free(sp->a);
    free(sp->z);
    free(sp->x);
    free(sp->w);
#ifndef EIGEN_ELPA
    free(sp->work);
    free(sp->iwork);
#endif
    Cblacs_gridexit(sp->ctxt);
  }
  if (eig->comm->rank == 0) {
    Cblacs_gridexit(sp->ctxt0);
  }
  commFreeSubsetFortran(sp->fcomm);
  free(sp);
}

void eigenSolve(EigenType *eig, double *overlap, double *eval, double *xmat)
{
  ScalapackType *sp = eig->priv;
  int n = eig->norb, nb = EIGEN_NB, one = 1, info = 0;
  double done = 1.0, dzero = 0.0;
  (void)eval;
  (void)xmat;

  if (sp->ctxt < 0) {
    return;
  }
  // every group holds the full S: copy the local blocks, the lower triangle by symmetry
  for (int jl = 0; jl < sp->nloc; jl++) {
    int j = globalIndex(jl, nb, sp->mycol, eig->grid[1]);
    for (int il = 0; il < sp->mloc; il++) {
      int i = globalIndex(il, nb, sp->myrow, eig->grid[0]);
      sp->a[il + (size_t)jl * sp->lld] =
          (i <= j) ? overlap[i + (size_t)j * n] : overlap[j + (size_t)i * n];
    }
  }

#ifdef EIGEN_ELPA
  elpaSolve(
      n, nb, sp->mloc, sp->nloc, sp->myrow, sp->mycol, sp->fcomm, sp->a, sp->w, sp->z);
#else
  pdsyevd_("V",
      "U",
      &n,
      sp->a,
      &one,
      &one,
      sp->desc,
      sp->w,
      sp->z,
      &one,
      &one,
      sp->desc,
      sp->work,
      &sp->lwork,
      sp->iwork,
      &sp->liwork,
      &info);
#endif
  if (info != 0) {
    dftAbortRank("ERROR: info=%d in pdsyevd used for Loewdin orthogonalization\n", info);
  }
  if (sp->w[0] <= 0.0) {
    dftAbortRank(
        "ERROR: overlap matrix is not positive definite (smallest eigenvalue %e)\n",
        sp->w[0]);
  }

  // W = U diag(eval^-1/4), X = W W^T
  for (int jl = 0; jl < sp->nloc; jl++) {
    int j      = globalIndex(jl, nb, sp->mycol, eig->grid[1]);
    double fac = 1.0 / sqrt(sqrt(sp->w[j]));
    for (int il = 0; il < sp->mloc; il++) {
      sp->z[il + (size_t)jl * sp->lld] *= fac;
    }
  }
  pdsyrk_("U",
      "N",
      &n,
      &n,
      &done,
      sp->z,
      &one,
      &one,
      sp->desc,
      &dzero,
      sp->x,
      &one,
      &one,
      sp->desc);
}

double eigenGather(EigenType *eig, double *eval, double *xmat)
{
  ScalapackType *sp = eig->priv;
  const CommType *c = eig->comm;
  int n = eig->norb, one = 1;
  double bytes = 0.0;

  if (sp->ctxt >= 0) {
    pdgemr2d_(
        &n, &n, sp->x, &one, &one, sp->desc, xmat, &one, &one, sp->desc0, &sp->ctxt);
    if (c->group != 0) {
      bytes = (double)sp->mloc * sp->nloc * sizeof(double);
    }
  }
  if (c->rank == 0) {
    // pdsyrk computed the upper triangle
    for (int j = 0; j < n; j++) {
      for (int i = j + 1; i < n; i++) {
        xmat[i + (size_t)j * n] = xmat[j + (size_t)i * n];
      }
    }
    memcpy(eval, sp->w, n * sizeof(double));
  }
  return bytes + eigenBroadcast(c, n, eval, xmat);
}
