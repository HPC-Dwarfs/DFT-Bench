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

/* Tolerances of the validation checks (-c) */
#define TOL_GTO_RESIDUAL 1.0E-12 /* max |psi - ref| / max |ref| */
#define TOL_NORM 1.0E-8 /* |1 - <psi|psi>| */
#define TOL_EKIN 1.0E-6 /* relative error of the kinetic energy */
#define TOL_ORTHO 1.0E-10 /* max |S - I| after orthogonalization */

#endif // __CONSTANTS_H_
