/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of DFT-Bench.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "constants.h"
#include "dft.h"
#include "gtoOnGrid.h"
#include "poscarReader.h"

void test_put_gto_sym_ortho(double *rxyz,
    double gw,
    double *xyz111,
    int ngx,
    int ngy,
    int ngz,
    double *hgrid,
    double *orb_i)
{
  int ll;
  double *waT, *wa, *orb;
  double gwsqinv, fac, res;
  double hx, hy, hz;
  double *expx, *expy, *expz;
  int negx, negy, negz;
  negx    = 5 * ngx;
  negy    = 5 * ngy;
  negz    = 5 * ngz;
  orb     = malloc(ngx * ngy * ngz * sizeof(double));
  waT     = calloc(negx * negy * negz, sizeof(double));

  wa      = (waT + negx * negy * ngz * 2 + negx * ngy * 2 + ngx * 2);
  hx      = hgrid[0];
  hy      = hgrid[4];
  hz      = hgrid[8];
  gwsqinv = 1.0 / (gw * gw);
  fac     = 1.0 / pow(gw * sqrt(M_PI), 3);
  fac     *= sqrt(pow(gw*sqrt(2.0*M_PI),3));
  printf("fac= %lf\n", fac);
  printf("rxyz  %20.10lf  %20.10lf  %20.10lf\n", rxyz[0], rxyz[1], rxyz[2]);

  // Precompute 1D exponentials (separable Gaussian)
  expx = malloc(negx * sizeof(double));
  expy = malloc(negy * sizeof(double));
  expz = malloc(negz * sizeof(double));
  for (int ix = -2 * ngx; ix < 3 * ngx; ix++) {
    double dx          = rxyz[0] - ix * hx - xyz111[0];
    expx[ix + 2 * ngx] = exp(-dx * dx * gwsqinv);
  }
  for (int iy = -2 * ngy; iy < 3 * ngy; iy++) {
    double dy          = rxyz[1] - iy * hy - xyz111[1];
    expy[iy + 2 * ngy] = exp(-dy * dy * gwsqinv);
  }
  for (int iz = -2 * ngz; iz < 3 * ngz; iz++) {
    double dz          = rxyz[2] - iz * hz - xyz111[2];
    expz[iz + 2 * ngz] = exp(-dz * dz * gwsqinv);
  }

#pragma omp parallel for schedule(static)
  for (int iz = -2 * ngz; iz < 3 * ngz; iz++) {
    double ez = fac * expz[iz + 2 * ngz];
    int kkZ   = negx * negy * iz;
    for (int iy = -2 * ngy; iy < 3 * ngy; iy++) {
      double eyz = ez * expy[iy + 2 * ngy];
      int kkYz   = kkZ + negx * iy;
      for (int ix = -2 * ngx; ix < 3 * ngx; ix++) {
        wa[kkYz + ix] = eyz * expx[ix + 2 * ngx];
      }
    }
  }

  memset(orb, 0, ngx * ngy * ngz * sizeof(double));

  for (int icz = 0; icz < 5; icz++) {
    for (int icy = 0; icy < 5; icy++) {
      for (int icx = 0; icx < 5; icx++) {
        for (int iz = 0; iz < ngz; iz++) {
          int kkZ = ngx * ngy * iz;
          int llZ = negx * negy * (icz * ngz + iz);
          for (int iy = 0; iy < ngy; iy++) {
            int kkYz = kkZ + ngx * iy;
            ll       = llZ + negx * (icy * ngy + iy) + icx * ngx;

            for (int ix = 0; ix < ngx; ix++) {
              orb[kkYz + ix] += waT[ll + ix];
            }
          }
        }
      }
    }
  }

  res = 0.0;
#pragma omp parallel for reduction(+ : res) schedule(static)
  for (int i = 0; i < ngx * ngy * ngz; i++) {
    double diff = orb[i] - orb_i[i];
    res += diff * diff;
  }

  printf("res= %14.5E\n", sqrt(res));
  free(expx);
  free(expy);
  free(expz);
  free(waT);
  free(orb);
}
void runDft(PoscarFileType *pf)
{
  int ngx, ngy, ngz;
  int igxt, igyt, igzt;
  double xyz111[3];
  double hgrid[3][3];
  double *orb, *orbk, *kinorb;
  double gw, rgcut, anorm;
  double akx, aky, akz, aknorm2, ekin, ekinAnalytical;
  double *rat     = *pf->rat_o;
  char (*sat)[5]  = *pf->sat_o;
  double (*cv)[3] = (double (*)[3])pf->cellvec;
  electronState es;
  fftw_plan planf, planb;

  gw              = 1.11 / BOHR2ANG * 0.5;
  rgcut           = 6.0 * gw;
  ngx             = 32;
  ngy             = 32;
  ngz             = 32;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      hgrid[i][j] = 0.0;
    }
  }

  hgrid[0][0] = cv[0][0] / ngx;
  hgrid[1][1] = cv[1][1] / ngy;
  hgrid[2][2] = cv[2][2] / ngz;
  printf("hgrid= %10.5f %10.5f %10.5f\n",hgrid[0][0],hgrid[1][1],hgrid[2][2]);

  xyz111[0]   = 0.0;
  xyz111[1]   = 0.0;
  xyz111[2]   = 0.0;

  es.norb=2*(*pf->nat_o);
  es.orbs=(double **) malloc(es.norb * sizeof(double *));
  for(int iorb=0;iorb<es.norb;iorb++) {
    es.orbs[iorb]=malloc(ngx * ngy * ngz * sizeof(double));
  }

  //orb         = malloc(ngx * ngy * ngz * sizeof(double));
  orb=es.orbs[0];
  orbk        = malloc((ngx+2) * ngy * ngz * sizeof(double));
  kinorb      = malloc(ngx * ngy * ngz * sizeof(double));

  putGtoSymOrtho(&rat[3 * 3], gw, rgcut, xyz111, ngx, ngy, ngz, &hgrid[0][0], orb);
  for(int i=0;i<ngx*ngy*ngz;i++) orb[i]=orb[i]*sqrt(pow(gw*sqrt(2.0*M_PI),3));
  anorm=0.0;
  for(int i=0;i<ngx*ngy*ngz;i++) anorm+=orb[i]*orb[i];
  anorm=sqrt(anorm*hgrid[0][0]*hgrid[1][1]*hgrid[2][2]);
  printf("anorm= %20.10f\n",anorm);

  planf=fftw_plan_dft_r2c_3d(ngx, ngy, ngz, orb, (fftw_complex *)orbk, FFTW_ESTIMATE);
  planb=fftw_plan_dft_c2r_3d (ngx, ngy, ngz, (fftw_complex *)orbk, kinorb, FFTW_ESTIMATE);

  fftw_execute(planf);
  for(int igz=0;igz<ngz;igz++) {
    for(int igy=0;igy<ngy;igy++) {
      for(int igx=0;igx<ngx;igx++) {
        igxt=igx;
        igyt=igy;
        igzt=igz;
        if(igy>ngy/2-1) igyt=ngy-igy;
        if(igz>ngz/2-1) igzt=ngz-igz;
        akx=2.0*M_PI*(igxt/2) / ((double) ngx * hgrid[0][0]);
        aky=2.0*M_PI*(igyt  ) / ((double) ngy * hgrid[1][1]);
        akz=2.0*M_PI*(igzt  ) / ((double) ngz * hgrid[2][2]);
        aknorm2=akx*akx+aky*aky+akz*akz;
        orbk[igz*(ngx+2)*ngy+igy*(ngx+2)+igx]*=0.5*aknorm2;
      }
    }
  }
  ekinAnalytical=3.0/(4.0*pow(gw,5)*M_PI*sqrt(2.0*M_PI))*(pow(gw*sqrt(2.0*M_PI),3));
  fftw_execute(planb);
  ekin=0.0;
  for(int i=0;i<ngx*ngy*ngz;i++) kinorb[i]*=1.0/(double)(ngx*ngy*ngz);
  for(int i=0;i<ngx*ngy*ngz;i++) ekin+=orb[i]*kinorb[i];
  ekin*=(hgrid[0][0]*hgrid[1][1]*hgrid[2][2]);
  printf("nat= %d\n",*pf->nat_o);
  printf("ekin= %20.10f  %20.10f  %14.5E\n",ekin,ekinAnalytical,ekin-ekinAnalytical);

  test_put_gto_sym_ortho(&rat[3 * 3], gw, xyz111, ngx, ngy, ngz, &hgrid[0][0], orb);

  fftw_destroy_plan(planf);
  fftw_destroy_plan(planb);

  for(int iorb=0;iorb<es.norb;iorb++) {
    free(es.orbs[iorb]);
  }
  free(es.orbs);
  //free(orb);
  free(orbk);
  free(kinorb);
}
