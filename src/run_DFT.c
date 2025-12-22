#include "constants.h"
#include "gto_on_grid.h"
#include "io_vasp.h"
#include "run_DFT.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void run_DFT(void)
{
  char filename[] = "poscars/POSCAR-0008-mp-149";
  int nat, ngx, ngy, ngz;
  double *rat = NULL;
  char (*sat)[5];
  double cellvec[3][3];
  double xyz111[3];
  double hgrid[3][3];
  double *orb;
  double gw, rgcut;
  char bc[] = "bulk";
  read_poscar(filename, &nat, &rat, &sat, &cellvec[0][0]);
  gw    = 1.11 / BOHR2ANG;
  rgcut = 6.0 * gw;
  ngx   = 128;
  ngy   = 128;
  ngz   = 128;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      hgrid[i][j] = 0.0;
  hgrid[0][0] = cellvec[0][0] / ngx;
  hgrid[1][1] = cellvec[1][1] / ngx;
  hgrid[2][2] = cellvec[2][2] / ngx;
  xyz111[0]   = 0.0;
  xyz111[1]   = 0.0;
  xyz111[2]   = 0.0;
  orb         = malloc(ngx * ngy * ngy * sizeof(double));
  put_gto_sym_ortho(bc, &rat[3 * 3], gw, rgcut, xyz111, ngx, ngy, ngz, &hgrid[0][0], orb);
  test_put_gto_sym_ortho(
      &rat[3 * 3], gw, rgcut, xyz111, ngx, ngy, ngz, &hgrid[0][0], orb);

  free(rat);
  free(sat);
  free(orb);
}
void test_put_gto_sym_ortho(double *rxyz,
    double gw,
    double rgcut,
    double *xyz111,
    int ngx,
    int ngy,
    int ngz,
    double *hgrid,
    double *orb_i)
{
  int kk, ll;
  double *wa_t, *wa, *orb;
  double gwsqinv, fac, pi, res;
  double hx, hy, hz, dx, dy, dz;
  int negx, negy, negz;
  negx = 5 * ngx;
  negy = 5 * ngy;
  negz = 5 * ngz;
  orb  = malloc(ngx * ngy * ngz * sizeof(double));
  wa_t = malloc(negx * negy * negz * sizeof(double));
  for (int i = 0; i < negx * negy * negz; i++)
    wa_t[i] = 0.0;
  wa      = (wa_t + negx * negy * ngz * 2 + negx * ngy * 2 + ngx * 2);
  hx      = hgrid[0];
  hy      = hgrid[4];
  hz      = hgrid[8];
  gwsqinv = 1.0 / (gw * gw);
  pi      = 4.0 * atan(1.0);
  fac     = 1.0 / pow(gw * sqrt(pi), 3);
  printf("fac= %lf\n", fac);
  printf("rxyz  %20.10lf  %20.10lf  %20.10lf\n", rxyz[0], rxyz[1], rxyz[2]);

  for (int iz = -2 * ngz; iz < 3 * ngz; iz++) {
    for (int iy = -2 * ngy; iy < 3 * ngy; iy++) {
      for (int ix = -2 * ngx; ix < 3 * ngx; ix++) {
        dx     = rxyz[0] - (ix - 0) * hx - xyz111[0];
        dy     = rxyz[1] - (iy - 0) * hy - xyz111[1];
        dz     = rxyz[2] - (iz - 0) * hz - xyz111[2];
        kk     = negx * negy * iz + negx * iy + ix;
        wa[kk] = fac * exp(-(dx * dx + dy * dy + dz * dz) * gwsqinv);
      }
    }
  }
  for (int i = 0; i < ngx * ngy * ngz; i++)
    orb[i] = 0.0;
  for (int icz = 0; icz < 5; icz++) {
    for (int icy = 0; icy < 5; icy++) {
      for (int icx = 0; icx < 5; icx++) {
        for (int iz = 0; iz < ngz; iz++) {
          for (int iy = 0; iy < ngy; iy++) {
            ll = negx * negy * (icz * ngz + iz) + negx * (icy * ngy + iy) + icx * ngx;
            for (int ix = 0; ix < ngx; ix++) {
              kk = ngx * ngy * iz + ngx * iy + ix;
              orb[kk] += wa_t[ll + ix];
            }
          }
        }
      }
    }
  }
  res = 0.0;
  for (int iz = 0; iz < ngz; iz++) {
    for (int iy = 0; iy < ngy; iy++) {
      for (int ix = 0; ix < ngx; ix++) {
        kk = ngx * ngy * iz + ngx * iy + ix;
        res += pow(orb[kk] - orb_i[kk], 2);
      }
    }
  }
  printf("res= %14.5E\n", sqrt(res));
  free(wa_t);
  free(orb);
}
