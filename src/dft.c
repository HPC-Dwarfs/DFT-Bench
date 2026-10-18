/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of DFT-Bench.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#include <fftw3.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "constants.h"
#include "dft.h"
#include "gtoOnGrid.h"
#include "mkl_cblas.h"
#include "poscarReader.h"

extern void dsyev(char *jobz,
    char *uplo,
    int *n,
    double *a,
    int *lda,
    double *w,
    double *work,
    int *lwork,
    int *info);
extern void dgemm(char *transa,
    char *transb,
    int *m,
    int *n,
    int *k,
    double *alpha,
    const double *a,
    int *lda,
    const double *b,
    int *ldb,
    double *beta,
    double *c,
    int *ldc);

void testPutGtoSymOrthoOrb(double *rxyz,
    double gw,
    double *xyz111,
    int ngx,
    int ngy,
    int ngz,
    double *hgrid,
    const double *orb_i)
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
  fac *= sqrt(pow(gw * sqrt(2.0 * M_PI), 3));
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
double calKineticEnergyOrbSingle(int ngx,
    int ngy,
    int ngz,
    double (*hgrid)[3],
    const double *orb,
    double *orbk,
    double *kinorb,
    fftw_plan *planf,
    fftw_plan *planb)
{
  int igxt, igyt, igzt;
  double akx, aky, akz, aknorm2, ekin;
  fftw_execute(*planf);
  for (int igz = 0; igz < ngz; igz++) {
    for (int igy = 0; igy < ngy; igy++) {
      for (int igx = 0; igx < ngx; igx++) {
        igxt = igx;
        igyt = igy;
        igzt = igz;
        if (igy > ngy / 2 - 1)
          igyt = ngy - igy;
        if (igz > ngz / 2 - 1)
          igzt = ngz - igz;
        akx     = 2.0 * M_PI * (igxt / 2) / ((double)ngx * hgrid[0][0]);
        aky     = 2.0 * M_PI * (igyt) / ((double)ngy * hgrid[1][1]);
        akz     = 2.0 * M_PI * (igzt) / ((double)ngz * hgrid[2][2]);
        aknorm2 = akx * akx + aky * aky + akz * akz;
        orbk[igz * (ngx + 2) * ngy + igy * (ngx + 2) + igx] *= 0.5 * aknorm2;
      }
    }
  }
  fftw_execute(*planb);
  for (int i = 0; i < ngx * ngy * ngz; i++)
    kinorb[i] *= 1.0 / (double)(ngx * ngy * ngz);
  ekin = 0.0;
  for (int i = 0; i < ngx * ngy * ngz; i++)
    ekin += orb[i] * kinorb[i];
  ekin *= (hgrid[0][0] * hgrid[1][1] * hgrid[2][2]);
  return ekin;
}
void getOverlap(electronState *es, double voxel, int ngtot, double *overlap)
{
  for (int iorb = 0; iorb < es->norb; iorb++) {
    for (int jorb = 0; jorb <= iorb; jorb++) {
      double ttdot;
      ttdot = cblas_ddot(ngtot, es->orbs[iorb], 1, es->orbs[jorb], 1);
      overlap[iorb * es->norb + jorb] = ttdot * voxel;
    }
  }
  for (int iorb = 0; iorb < es->norb; iorb++) {
    for (int jorb = 0; jorb <= iorb; jorb++) {
      overlap[jorb * es->norb + iorb] = overlap[iorb * es->norb + jorb];
    }
  }
}
void orthogonalizeLoewdin(electronState *es, double voxel, int ngx, int ngy, int ngz)
{
  double *overlap;
  overlap   = (double *)malloc(es->norb * es->norb * sizeof(double));
  int ngtot = ngx * ngy * ngz;
  getOverlap(es, voxel, ngtot, overlap);
  for (int iorb = 0; iorb < es->norb; iorb++) {
    for (int jorb = 0; jorb < es->norb; jorb++) {
      printf("%6.2f", overlap[iorb * es->norb + jorb]);
    }
    printf("\n");
  }
  int info, nwa;
  double *wa1, *wa2, *eval;
  nwa       = es->norb * es->norb + 100;
  wa1       = (double *)malloc(nwa * sizeof(double));
  wa2       = (double *)malloc(nwa * sizeof(double));
  eval      = (double *)malloc(es->norb * sizeof(double));
  char jobz = 'V';
  char uplo = 'U';
  dsyev(&jobz, &uplo, &es->norb, overlap, &es->norb, eval, wa1, &nwa, &info);
  //for (int iorb = 0; iorb < es->norb; iorb++) {
  //  printf("evec %4d",iorb+1);
  //  for (int jorb = 0; jorb < es->norb; jorb++) {
  //    printf("%6.2f",overlap[iorb*es->norb+jorb]);
  //  }
  //  printf("\n");
  //}
  if (info != 0) {
    printf("ERROR: info=%d in dsyev used for Loewdin orthogonalization", info);
    exit(1);
  }
  for (int iorb = 0; iorb < es->norb; iorb++) {
    printf("%5d%8.3f\n", iorb, eval[iorb]);
  }
  for (int i = 0; i < es->norb * es->norb; i++) {
    wa1[i] = overlap[i];
  }
  for (int iorb = 0; iorb < es->norb; iorb++) {
    for (int jorb = 0; jorb < es->norb; jorb++) {
      overlap[iorb * es->norb + jorb] = wa1[jorb * es->norb + iorb];
    }
  }
  for (int iorb = 0; iorb < es->norb; iorb++) {
    for (int jorb = 0; jorb < es->norb; jorb++) {
      wa1[iorb * es->norb + jorb] =
          overlap[iorb * es->norb + jorb] * sqrt(1.0 / eval[jorb]);
    }
  }
  char noTrans = 'N';
  char doTrans = 'T';
  double one   = 1.0;
  double zero  = 0.0;
  dgemm(&doTrans,
      &noTrans,
      &es->norb,
      &es->norb,
      &es->norb,
      &one,
      wa1,
      &es->norb,
      overlap,
      &es->norb,
      &zero,
      wa2,
      &es->norb);
  double *orbst1, *orbst2;
  orbst1 = (double *)malloc((ngx * ngy * ngz * es->norb) * sizeof(double));
  orbst2 = (double *)malloc((ngx * ngy * ngz * es->norb) * sizeof(double));
  for (int iorb = 0; iorb < es->norb; iorb++) {
    for (int i = 0; i < ngx * ngy * ngz; i++) {
      orbst1[i * es->norb + iorb] = es->orbs[iorb][i];
    }
  }
  dgemm(&doTrans,
      &noTrans,
      &es->norb,
      &ngtot,
      &es->norb,
      &one,
      wa2,
      &es->norb,
      orbst1,
      &es->norb,
      &zero,
      orbst2,
      &es->norb);
  for (int iorb = 0; iorb < es->norb; iorb++) {
    for (int i = 0; i < ngx * ngy * ngz; i++) {
      es->orbs[iorb][i] = orbst2[i * es->norb + iorb];
    }
  }
  //for (int iorb = 0; iorb < es->norb; iorb++) {
  //  printf("orbst2 %4d\n",iorb+1);
  //  for (int i = 0; i < ngx * ngy * ngz; i++) {
  //    printf("%15.5f\n",orbst2[i * es->norb + iorb]);
  //  }
  //}
  printf("-----------------------------------------------------------------\n");
  getOverlap(es, voxel, ngtot, overlap);
  for (int iorb = 0; iorb < es->norb; iorb++) {
    for (int jorb = 0; jorb < es->norb; jorb++) {
      double tt;
      tt = overlap[iorb * es->norb + jorb];
      if (fabs(tt) < 1.E-12)
        tt = 0.0;
      printf("%6.2f", tt);
    }
    printf("\n");
  }
  free(orbst1);
  free(orbst2);
  free(overlap);
  free(wa1);
  free(wa2);
  free(eval);
}
void runDft(PoscarFileType *pf)
{
  int ngx, ngy, ngz;
  double xyz111[3];
  double hgrid[3][3];
  double **orbk, **kinorb;
  double gw1, gw2, rgcut, anorm;
  double ekin, ekinAnalytical;
  double *rat     = *pf->rat_o;
  char (*sat)[5]  = *pf->sat_o;
  double (*cv)[3] = (double (*)[3])pf->cellvec;
  electronState es;
  fftw_plan *planf, *planb;

  gw1   = 1.11 / BOHR2ANG * 0.5;
  gw2   = 1.11 / BOHR2ANG * 0.6;
  rgcut = 6.0 * fmax(gw1, gw2);
  ngx   = 32;
  ngy   = 32;
  ngz   = 32;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      hgrid[i][j] = 0.0;
    }
  }

  hgrid[0][0] = cv[0][0] / ngx;
  hgrid[1][1] = cv[1][1] / ngy;
  hgrid[2][2] = cv[2][2] / ngz;
  printf("hgrid= %10.5f %10.5f %10.5f\n", hgrid[0][0], hgrid[1][1], hgrid[2][2]);

  xyz111[0] = 0.0;
  xyz111[1] = 0.0;
  xyz111[2] = 0.0;

  es.norb   = 2 * (*pf->nat_o);
  planf     = (fftw_plan *)malloc(es.norb * sizeof(fftw_plan));
  planb     = (fftw_plan *)malloc(es.norb * sizeof(fftw_plan));
  es.gwall  = (double *)malloc(es.norb * sizeof(double));
  es.orbs   = (double **)malloc(es.norb * sizeof(double *));
  for (int iorb = 0; iorb < es.norb; iorb++) {
    es.orbs[iorb] = malloc(ngx * ngy * ngz * sizeof(double));
  }
  for (int iorb = 0; iorb < es.norb; iorb++) {
    if (iorb % 2 == 0) {
      es.gwall[iorb] = gw1;
    } else {
      es.gwall[iorb] = gw2;
    }
  }

  orbk = (double **)malloc(es.norb * sizeof(double *));
  for (int iorb = 0; iorb < es.norb; iorb++) {
    orbk[iorb] = malloc((ngx + 2) * ngy * ngz * sizeof(double));
  }
  kinorb = (double **)malloc(es.norb * sizeof(double *));
  for (int iorb = 0; iorb < es.norb; iorb++) {
    kinorb[iorb] = malloc(ngx * ngy * ngz * sizeof(double));
  }

  putGtoSymOrthoOrb(
      &rat[3 * 3], es.gwall[0], rgcut, xyz111, ngx, ngy, ngz, &hgrid[0][0], es.orbs[0]);
  anorm = 0.0;
  for (int i = 0; i < ngx * ngy * ngz; i++)
    anorm += es.orbs[0][i] * es.orbs[0][i];
  anorm = sqrt(anorm * hgrid[0][0] * hgrid[1][1] * hgrid[2][2]);
  printf("anorm= %20.10f\n", anorm);
  testPutGtoSymOrthoOrb(
      &rat[3 * 3], es.gwall[0], xyz111, ngx, ngy, ngz, &hgrid[0][0], es.orbs[0]);

  for (int iorb = 0; iorb < es.norb; iorb++) {
    putGtoSymOrthoOrb(&rat[iorb / 2],
        es.gwall[iorb],
        rgcut,
        xyz111,
        ngx,
        ngy,
        ngz,
        &hgrid[0][0],
        es.orbs[iorb]);
  }
  for (int iorb = 0; iorb < es.norb; iorb++) {
    planf[iorb] = fftw_plan_dft_r2c_3d(
        ngx, ngy, ngz, es.orbs[iorb], (fftw_complex *)orbk[iorb], FFTW_ESTIMATE);
    planb[iorb] = fftw_plan_dft_c2r_3d(
        ngx, ngy, ngz, (fftw_complex *)orbk[iorb], kinorb[iorb], FFTW_ESTIMATE);
  }

  ekinAnalytical = 0.0;
  ekin           = 0.0;
  for (int iorb = 0; iorb < es.norb; iorb++) {
    double ekinSingle;
    ekinSingle = calKineticEnergyOrbSingle(ngx,
        ngy,
        ngz,
        hgrid,
        es.orbs[iorb],
        orbk[iorb],
        kinorb[iorb],
        &planf[iorb],
        &planb[iorb]);
    ekin += ekinSingle;
    ekinAnalytical += 3.0 / (4.0 * pow(es.gwall[iorb], 5) * M_PI * sqrt(2.0 * M_PI)) *
                      (pow(es.gwall[iorb] * sqrt(2.0 * M_PI), 3));
  }
  printf("nat= %d\n", *pf->nat_o);
  printf("ekin= %20.10f  %20.10f  %14.5E\n", ekin, ekinAnalytical, ekin - ekinAnalytical);

  double voxel = hgrid[0][0] * hgrid[1][1] * hgrid[2][2];
  orthogonalizeLoewdin(&es, voxel, ngx, ngy, ngz);

  for (int iorb = 0; iorb < es.norb; iorb++) {
    fftw_destroy_plan(planf[iorb]);
    fftw_destroy_plan(planb[iorb]);
  }

  free(es.gwall);
  for (int iorb = 0; iorb < es.norb; iorb++) {
    free(es.orbs[iorb]);
  }
  free(es.orbs);
  for (int iorb = 0; iorb < es.norb; iorb++) {
    free(orbk[iorb]);
  }
  free(orbk);
  for (int iorb = 0; iorb < es.norb; iorb++) {
    free(kinorb[iorb]);
  }
  free(kinorb);
  free(planf);
  free(planb);
}
