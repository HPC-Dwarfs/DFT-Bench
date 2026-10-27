/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of DFT-Bench.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#ifndef __CONSTANTS_H_
#define __CONSTANTS_H_

#define BOHR2ANG 0.529177210

#ifndef ARRAY_ALIGNMENT
#define ARRAY_ALIGNMENT 64
#endif

/* Defaults for the command line options */
#define DEFAULT_POSCAR "poscars/POSCAR-0008-mp-149"
#define DEFAULT_GRID_SPACING 0.3 /* bohr */
#define DEFAULT_REPETITIONS 3

/* Row block size of the blocked orbital rotation Psi <- Psi * S^{-1/2} */
#define ROTATE_BLOCK 4096

/* Edge length of the tiles of S that one rank computes in the overlap phase */
#ifndef OVERLAP_TILE
#define OVERLAP_TILE 128
#endif

/* Distributed eigensolvers: block size of the block-cyclic distribution, and the
 * default process grid holds about EIGEN_BLOCK_TARGET^2 elements per grid rank */
#ifndef EIGEN_NB
#define EIGEN_NB 32
#endif
#ifndef EIGEN_BLOCK_TARGET
#define EIGEN_BLOCK_TARGET 160
#endif
/* Largest matrix that is checked against the LAPACK reference with -c */
#define EIGEN_REFERENCE_MAX 4096

/* Default memory budget of the MPI build: fraction of the node's physical memory */
#define DEFAULT_MEMORY_FRACTION 0.8

/* Tolerances of the validation checks (-c) */
#define TOL_GTO_RESIDUAL 1.0E-12 /* max |psi - ref| / max |ref| */
#define TOL_NORM 1.0E-8 /* |1 - <psi|psi>| */
#define TOL_EKIN 1.0E-6 /* relative error of the kinetic energy */
#define TOL_ORTHO 1.0E-10 /* max |S - I| after orthogonalization */
#define TOL_EIGEN_REF 1.0E-10 /* distributed eigensolver vs. LAPACK reference */

#endif // __CONSTANTS_H_
