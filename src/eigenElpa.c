/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of DFT-Bench.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#include <elpa/elpa.h>

#include "comm.h"

/* ELPA eigensolver of the ELPA backend. Grid, data distribution and the rest of the
 * eigen phase are shared with the ScaLAPACK backend in eigenScalapack.c. */

#define ELPA_API 20171201

static void elpaCheck(int error, const char *what)
{
  if (error != ELPA_OK) {
    dftAbortRank("ERROR: %s failed in ELPA (error %d)\n", what, error);
  }
}

/* Eigenvalues w and eigenvectors z of the block-cyclic matrix a (destroyed), with the
 * 2-stage solver */
void elpaSolve(int n,
    int nb,
    int mloc,
    int nloc,
    int myrow,
    int mycol,
    int fcomm,
    double *a,
    double *w,
    double *z)
{
  static elpa_t handle = NULL;
  int error;

  if (handle == NULL) {
    elpaCheck(elpa_init(ELPA_API), "elpa_init");
    handle = elpa_allocate(&error);
    elpaCheck(error, "elpa_allocate");
    elpa_set(handle, "na", n, &error);
    elpaCheck(error, "na");
    elpa_set(handle, "nev", n, &error);
    elpaCheck(error, "nev");
    elpa_set(handle, "local_nrows", mloc, &error);
    elpaCheck(error, "local_nrows");
    elpa_set(handle, "local_ncols", nloc, &error);
    elpaCheck(error, "local_ncols");
    elpa_set(handle, "nblk", nb, &error);
    elpaCheck(error, "nblk");
    elpa_set(handle, "mpi_comm_parent", fcomm, &error);
    elpaCheck(error, "mpi_comm_parent");
    elpa_set(handle, "process_row", myrow, &error);
    elpaCheck(error, "process_row");
    elpa_set(handle, "process_col", mycol, &error);
    elpaCheck(error, "process_col");
    elpaCheck(elpa_setup(handle), "elpa_setup");
    elpa_set(handle, "solver", ELPA_SOLVER_2STAGE, &error);
    elpaCheck(error, "solver");
  }
  elpa_eigenvectors(handle, a, w, z, &error);
  elpaCheck(error, "elpa_eigenvectors");
}
