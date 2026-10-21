/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of DFT-Bench.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#ifndef __EIGEN_H_
#define __EIGEN_H_

#include <stdbool.h>

#include "comm.h"
#include "dft.h"

/* Backends of the eigen phase, selected with EIGEN_SOLVER at build time:
 * eigenLapack.c (LAPACK), eigenScalapack.c (SCALAPACK, and ELPA with eigenElpa.c).
 * All of them compute the eigenvalues of S and the Loewdin matrix X = S^{-1/2}. */
typedef struct {
  const CommType *comm;
  int norb;
  int grid[2]; /* process grid of the eigensolver */
  bool spread; /* grid ranks spread over the groups */
  double *wa; /* LAPACK: work array on world rank 0 */
  void *priv; /* backend data */
} EigenType;

const char *eigenBackendName(void);
/* Sets up the backend. Called by all ranks. */
void eigenInit(EigenType *eig, const CommType *comm, const ConfigType *cfg, int norb);
void eigenPrintInfo(const EigenType *eig);
void eigenFree(EigenType *eig);
/* Eigen phase. Called by all ranks with the group-shared upper triangle of S, which
 * the call may overwrite. On return world rank 0 holds the eigenvalues in ascending
 * order or the backend holds them distributed. */
void eigenSolve(EigenType *eig, double *overlap, double *eval, double *xmat);
/* Gather phase: stores the eigenvalues and the full X in the group-shared eval and
 * xmat of every group. Returns the bytes this rank received from other groups. */
double eigenGather(EigenType *eig, double *eval, double *xmat);

/* LAPACK path on the calling rank: X = S^{-1/2} = U diag(eval^{-1/2}) U^T. The overlap
 * is overwritten by U, wa holds norb x norb elements. */
void loewdinLapack(int norb, double *overlap, double *eval, double *wa, double *x);
/* Broadcasts eval and X from world rank 0 over the group leaders into every group's
 * windows and synchronizes the groups. Returns the bytes received by this rank. */
double eigenBroadcast(const CommType *comm, int norb, double *eval, double *xmat);

#endif // __EIGEN_H_
