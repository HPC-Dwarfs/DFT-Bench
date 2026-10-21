/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of DFT-Bench.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#ifndef __RUN_DFT_H_
#define __RUN_DFT_H_

#include <stdbool.h>
#include <stddef.h>

#include "poscarReader.h"

typedef struct {
  const char *poscarFile;
  double spacing; /* target grid spacing in bohr, used if ng[0] == 0 */
  int ng[3]; /* explicit grid size (0: derive from spacing) */
  int reps; /* number of timed repetitions */
  bool validate;
  bool verbose;
} ConfigType;

typedef struct {
  int ng[3]; /* grid points per direction, x runs fastest */
  double h[3]; /* grid spacing in bohr */
  double origin[3]; /* position of grid point (0,0,0) */
  size_t ngtot; /* ng[0] * ng[1] * ng[2] */
  size_t ld; /* leading dimension of the orbital block (ngtot padded) */
  double voxel;
} GridType;

typedef struct {
  int norb;
  int *iat; /* atom each orbital is centered on */
  double *gw; /* Gaussian width of each orbital in bohr */
  double *psi; /* column-major ld x norb block, orbital i at psi + i * ld */
} ElectronStateType;

/* Runs the benchmark, returns EXIT_SUCCESS or EXIT_FAILURE (failed validation). */
int runDft(const PoscarFileType *pf, const ConfigType *cfg);

#endif // __RUN_DFT_H_
