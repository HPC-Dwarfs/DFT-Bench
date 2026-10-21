/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of DFT-Bench.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#ifndef __RUN_DFT_H_
#define __RUN_DFT_H_

#include <stdbool.h>
#include <stddef.h>

#include "comm.h"
#include "poscarReader.h"

typedef struct {
  const char *poscarFile;
  double spacing; /* target grid spacing in bohr, used if ng[0] == 0 */
  int ng[3]; /* explicit grid size (0: derive from spacing) */
  int reps; /* number of timed repetitions */
  bool validate;
  bool verbose;
  int groupSize; /* MPI: max ranks per shared-memory group (0: node) */
  double memBudget; /* MPI: memory budget per group in GB (0: default) */
  int eigenGrid[2]; /* MPI: eigensolver process grid (0: default) */
  bool eigenSpread; /* MPI: spread the eigensolver grid over the groups */
} ConfigType;

typedef struct {
  int ng[3]; /* grid points per direction, x runs fastest */
  double h[3]; /* grid spacing in bohr */
  double origin[3]; /* position of grid point (0,0,0) */
  size_t ngtot; /* ng[0] * ng[1] * ng[2] */
  size_t ld; /* leading dimension of the orbital block (ngtot padded) */
  double voxel;
} GridType;

/* Decomposition of orbitals and grid rows. Both are split over all ranks in group
 * order (group 0 ranks 0..L-1, group 1, ...), so every group and every rank owns a
 * contiguous range. Remainders go to the first ranks, row counts are multiples of 8
 * except on the last rank. Arrays are indexed by world rank or by group. */
typedef struct {
  int *orb0, *norb; /* orbital range of every rank */
  size_t *row0, *nrows; /* grid-row range of every rank */
  int *gorb0, *gnorb; /* orbital range of every group */
  size_t *grow0, *gnrows; /* grid-row range of every group */
  int *members; /* world ranks of this rank's group, by group rank */
  int idle; /* ranks without orbitals */
} DistType;

typedef struct {
  int norb;
  int *iat; /* atom each orbital is centered on */
  double *gw; /* Gaussian width of each orbital in bohr */
  /* Orbital form, shared by the group: column-major ld x gnorb block of the group's
   * orbitals, orbital i at psi + (i - gorb0) * ld */
  double *psi;
  /* Transposed form, shared by the group: column-major ldt x norb block of the group's
   * grid rows of all orbitals. With a single group it is an alias of psi. */
  double *psit;
  size_t ldt;
} ElectronStateType;

/* Runs the benchmark, returns EXIT_SUCCESS or EXIT_FAILURE (failed validation). */
int runDft(const PoscarFileType *pf, const ConfigType *cfg, const CommType *comm);

#endif // __RUN_DFT_H_
