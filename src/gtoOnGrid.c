/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of DFT-Bench.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gtoOnGrid.h"
#include "util.h"

static int modulo(int n, int m)
{
  return ((n % m) + m) % m;
}

/* Adds the box wa (nbx*nby*nbz points, lower corner at grid index i0) to the
 * periodic cell rho. The box may be larger than the cell. */
static void chargeBackToCell(const int *ng,
    int nbx,
    int nby,
    int nbz,
    const int *i0,
    const double *wa,
    double *rho)
{
  int iiz = modulo(i0[2], ng[2]);
  for (int kz = 0; kz < nbz; kz++) {
    size_t kkZ = (size_t)ng[0] * ng[1] * iiz;
    int iiy    = modulo(i0[1], ng[1]);
    for (int ky = 0; ky < nby; ky++) {
      size_t kkYz       = kkZ + (size_t)ng[0] * iiy;
      const double *row = wa + ((size_t)kz * nby + ky) * nbx;
      int iix           = modulo(i0[0], ng[0]);
      // split the row into contiguous segments so that the inner loop vectorizes
      for (int kx = 0; kx < nbx;) {
        int len = MIN(nbx - kx, ng[0] - iix);
        for (int i = 0; i < len; i++) {
          rho[kkYz + iix + i] += row[kx + i];
        }
        kx += len;
        iix = 0;
      }
      if (++iiy == ng[1])
        iiy = 0;
    }
    if (++iiz == ng[2])
      iiz = 0;
  }
}

static void putGto(const double *rxyz,
    double gw,
    double rgcut,
    const double *origin,
    const int *ng,
    const double *h,
    double fac,
    double *rho)
{
  int nb[3], i0[3];
  double *ww[3];
  double widthInv = 1.0 / gw;

  for (int d = 0; d < 3; d++) {
    double alat = ng[d] * h[d];
    double xat  = rxyz[d] - origin[d];
    xat -= alat * floor(xat / alat); // wrap atom into the cell
    int iato = (int)lround(xat / h[d]);
    nb[d]    = (int)(rgcut / h[d]) + 2;
    i0[d]    = iato - nb[d];
    ww[d]    = malloc((2 * nb[d] + 1) * sizeof(double));
    if (ww[d] == NULL) {
      fprintf(stderr, "ERROR: out of memory in putGtoSymOrtho\n");
      exit(EXIT_FAILURE);
    }
    for (int i = 0; i <= 2 * nb[d]; i++) {
      double tt1 = widthInv * ((i0[d] + i) * h[d] - xat);
      ww[d][i]   = exp(-tt1 * tt1);
    }
  }

  int nbx = 2 * nb[0] + 1, nby = 2 * nb[1] + 1, nbz = 2 * nb[2] + 1;
  double *wa = malloc((size_t)nbx * nby * nbz * sizeof(double));
  if (wa == NULL) {
    fprintf(stderr, "ERROR: out of memory in putGtoSymOrtho\n");
    exit(EXIT_FAILURE);
  }
  for (int iz = 0; iz < nbz; iz++) {
    double rhoz = fac * ww[2][iz];
    for (int iy = 0; iy < nby; iy++) {
      double rhoyz = rhoz * ww[1][iy];
      double *row  = wa + ((size_t)iz * nby + iy) * nbx;
      for (int ix = 0; ix < nbx; ix++) {
        row[ix] = rhoyz * ww[0][ix];
      }
    }
  }
  memset(rho, 0, (size_t)ng[0] * ng[1] * ng[2] * sizeof(double));
  chargeBackToCell(ng, nbx, nby, nbz, i0, wa, rho);
  free(ww[0]);
  free(ww[1]);
  free(ww[2]);
  free(wa);
}

void putGtoSymOrtho(const double *rxyz,
    double gw,
    double rgcut,
    const double *origin,
    const int *ng,
    const double *h,
    double *rho)
{
  double fac = 1.0 / pow(gw * sqrt(M_PI), 3);
  putGto(rxyz, gw, rgcut, origin, ng, h, fac, rho);
}

void putGtoSymOrthoOrb(const double *rxyz,
    double gw,
    double rgcut,
    const double *origin,
    const int *ng,
    const double *h,
    double *orb)
{
  double fac = 1.0 / pow(gw * sqrt(M_PI), 3) * sqrt(pow(gw * sqrt(2.0 * M_PI), 3));
  putGto(rxyz, gw, rgcut, origin, ng, h, fac, orb);
}

double testPutGtoSymOrthoOrb(const double *rxyz,
    double gw,
    const double *origin,
    const int *ng,
    const double *h,
    const double *orb)
{
  double *f[3];
  double gwsqinv = 1.0 / (gw * gw);
  double fac     = 1.0 / pow(gw * sqrt(M_PI), 3) * sqrt(pow(gw * sqrt(2.0 * M_PI), 3));

  for (int d = 0; d < 3; d++) {
    double alat = ng[d] * h[d];
    int nimg    = (int)ceil(10.0 * gw / alat) + 1;
    f[d]        = malloc(ng[d] * sizeof(double));
    if (f[d] == NULL) {
      fprintf(stderr, "ERROR: out of memory in testPutGtoSymOrthoOrb\n");
      exit(EXIT_FAILURE);
    }
    for (int i = 0; i < ng[d]; i++) {
      double s = 0.0;
      for (int n = -nimg; n <= nimg; n++) {
        double dx = rxyz[d] - origin[d] - i * h[d] - n * alat;
        s += exp(-dx * dx * gwsqinv);
      }
      f[d][i] = s;
    }
  }

  double res = 0.0;
  for (int iz = 0; iz < ng[2]; iz++) {
    for (int iy = 0; iy < ng[1]; iy++) {
      double fyz        = fac * f[2][iz] * f[1][iy];
      const double *row = orb + ((size_t)iz * ng[1] + iy) * ng[0];
      for (int ix = 0; ix < ng[0]; ix++) {
        res = fmax(res, fabs(row[ix] - fyz * f[0][ix]));
      }
    }
  }
  free(f[0]);
  free(f[1]);
  free(f[2]);
  return res;
}
