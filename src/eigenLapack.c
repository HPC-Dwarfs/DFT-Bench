/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of DFT-Bench.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "allocate.h"
#include "blasLapack.h"
#include "comm.h"
#include "constants.h"
#include "eigen.h"

void loewdinLapack(int norb, double *overlap, double *eval, double *wa, double *x)
{
  lapack_int info = LAPACKE_dsyevd(LAPACK_COL_MAJOR, 'V', 'U', norb, overlap, norb, eval);
  if (info != 0) {
    dftAbortRank(
        "ERROR: info=%d in dsyevd used for Loewdin orthogonalization\n", (int)info);
  }
  if (eval[0] <= 0.0) {
    dftAbortRank(
        "ERROR: overlap matrix is not positive definite (smallest eigenvalue %e)\n",
        eval[0]);
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

double eigenBroadcast(const CommType *comm, int norb, double *eval, double *xmat)
{
  double bytes = 0.0;
  if (comm->ngroups > 1 && comm->grank == 0) {
    commBcast(xmat, (size_t)norb * norb * sizeof(double), 0, COMM_LEADERS);
    commBcast(eval, (size_t)norb * sizeof(double), 0, COMM_LEADERS);
    if (comm->group != 0) {
      bytes = ((double)norb * norb + norb) * sizeof(double);
    }
  }
  commSyncGroup();
  return bytes;
}

#ifdef EIGEN_LAPACK

const char *eigenBackendName(void)
{
  return "LAPACK";
}

void eigenInit(EigenType *eig, const CommType *comm, const ConfigType *cfg, int norb)
{
  eig->comm    = comm;
  eig->norb    = norb;
  eig->grid[0] = eig->grid[1] = 1;
  eig->spread                 = false;
  eig->priv                   = NULL;
  eig->wa                     = NULL;
  if (cfg->eigenGrid[0] > 0) {
    commPrintf("NOTICE: -e has no effect with the LAPACK eigensolver\n");
  }
  if (cfg->eigenSpread) {
    commPrintf("NOTICE: -p has no effect with the LAPACK eigensolver\n");
  }
  if (comm->rank == 0) {
    eig->wa = allocate(ARRAY_ALIGNMENT, (size_t)norb * norb * sizeof(double));
  }
}

void eigenPrintInfo(const EigenType *eig)
{
  (void)eig;
  if (commIsMpi()) {
    commPrintf("Eigensolver   : LAPACK on rank 0\n");
  }
}

void eigenFree(EigenType *eig)
{
  free(eig->wa);
}

void eigenSolve(EigenType *eig, double *overlap, double *eval, double *xmat)
{
  if (eig->comm->rank == 0) {
    loewdinLapack(eig->norb, overlap, eval, eig->wa, xmat);
  }
}

double eigenGather(EigenType *eig, double *eval, double *xmat)
{
  return eigenBroadcast(eig->comm, eig->norb, eval, xmat);
}

#endif /* EIGEN_LAPACK */
