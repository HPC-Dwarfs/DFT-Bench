/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of DFT-Bench.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#ifndef __GTO_ON_GRID_H_
#define __GTO_ON_GRID_H_

/* All routines assume an orthogonal cell. Grid point (ix,iy,iz) is located at
 * origin + (ix*h[0], iy*h[1], iz*h[2]) and stored at rho[ix + ng[0]*(iy + ng[1]*iz)].
 * The Gaussian is summed over all periodic images within the cutoff rgcut. */

/* Charge-normalized Gaussian: integrates to one. */
void putGtoSymOrtho(const double *rxyz,
    double gw,
    double rgcut,
    const double *origin,
    const int *ng,
    const double *h,
    double *rho);

/* Orbital-normalized Gaussian: its square integrates to one. */
void putGtoSymOrthoOrb(const double *rxyz,
    double gw,
    double rgcut,
    const double *origin,
    const int *ng,
    const double *h,
    double *orb);

/* Returns max |orb - reference| where the reference is evaluated independently
 * from 1D periodic image sums (valid for orthogonal cells). */
double testPutGtoSymOrthoOrb(const double *rxyz,
    double gw,
    const double *origin,
    const int *ng,
    const double *h,
    const double *orb);

#endif // __GTO_ON_GRID_H_
