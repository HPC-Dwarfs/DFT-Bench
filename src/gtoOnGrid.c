#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int modulo(int n, int m)
{
  return ((n % m) + m) % m;
}

static void chargeBackToCell(
    int ngx, int ngy, int ngz, int nagx, int nagy, int nagz, double *wa_i, double *rho)
{
  int nnx     = ngx + 2 * nagx;
  int nny     = ngy + 2 * nagy;
  double *wa  = wa_i + nnx * nny * (nagz - 1) + nnx * (nagy - 1) + nagx - 1;
  int iizInit = modulo(-nagz, ngz) - 1;
  int iiyInit = modulo(-nagy, ngy) - 1;
  int iixInit = modulo(-nagx, ngx) - 1;
  int iiz     = iizInit;
  for (int igz = 1 - nagz; igz <= ngz + nagz; igz++) {
    iiz++;
    if (iiz == ngz)
      iiz = 0;
    int kkZ = ngx * ngy * iiz;
    int llZ = nnx * nny * igz;
    int iiy = iiyInit;
    for (int igy = 1 - nagy; igy <= ngy + nagy; igy++) {
      iiy++;
      if (iiy == ngy)
        iiy = 0;
      int kkYz = kkZ + ngx * iiy;
      int llYz = llZ + nnx * igy;
      int iix  = iixInit;
      for (int igx = 1 - nagx; igx <= ngx + nagx; igx++) {
        iix++;
        if (iix == ngx)
          iix = 0;
        rho[kkYz + iix] += wa[llYz + igx];
      }
    }
  }
}

void putGtoSymOrtho(double *rxyz,
    double gw,
    double rgcut,
    double *xyz111,
    int ngx,
    int ngy,
    int ngz,
    double *hgrid,
    double *rho)
{
  double hgxinv, hgyinv, hgzinv;
  double widthInv, widthInvXyz[3];
  double widthInvHhh[3];
  double xat, yat, zat, fac;
  double hx, hy, hz, tt1;
  int iatox, iatoy, iatoz;
  int nwa;
  int nbgx, nbgy, nbgz, nagx, nagy, nagz;
  int nnx, nny;
  double *wa;
  bool ortho = !(hgrid[1] || hgrid[2] || hgrid[3] || hgrid[5] || hgrid[6] || hgrid[7]);
  if (!ortho) {
    printf("ERROR: this routine is only for orthogonal cell\n");
    printf("%24.15lf%24.15lf%24.15lf\n", hgrid[0], hgrid[1], hgrid[2]);
    printf("%24.15lf%24.15lf%24.15lf\n", hgrid[3], hgrid[4], hgrid[5]);
    printf("%24.15lf%24.15lf%24.15lf\n", hgrid[6], hgrid[7], hgrid[8]);
    exit(0);
  }
  hx          = hgrid[0];
  hy          = hgrid[4];
  hz          = hgrid[8];
  nbgx        = (int)(rgcut / hx) + 2;
  nbgy        = (int)(rgcut / hy) + 2;
  nbgz        = (int)(rgcut / hz) + 2;
  nagx        = nbgx + 1;
  nagy        = nbgy + 1;
  nagz        = nbgz + 1;
  nnx         = ngx + 2 * nagx;
  nny         = ngy + 2 * nagy;
  nwa         = nnx * nny * (ngz + 2 * nagz);
  wa          = calloc(nwa, sizeof(double));
  double *wwx = malloc((2 * nbgx + 1) * sizeof(double));
  double *wwy = malloc((2 * nbgy + 1) * sizeof(double));
  double *wwz = malloc((2 * nbgz + 1) * sizeof(double));
  hgxinv      = 1.0 / hx;
  hgyinv      = 1.0 / hy;
  hgzinv      = 1.0 / hz;
  iatox       = (int)lround((rxyz[0] - xyz111[0]) * hgxinv);
  iatoy       = (int)lround((rxyz[1] - xyz111[1]) * hgyinv);
  iatoz       = (int)lround((rxyz[2] - xyz111[2]) * hgzinv);
  xat         = rxyz[0] - (iatox - 1) * hx - xyz111[0];
  yat         = rxyz[1] - (iatoy - 1) * hy - xyz111[1];
  zat         = rxyz[2] - (iatoz - 1) * hz - xyz111[2];
  if (iatox < 1 - nagx || iatox > ngx + nagx) {
    printf("ERROR: charge of atom outside box in x-direction!\n");
    exit(0);
  }
  if (iatoy < 1 - nagy || iatoy > ngy + nagy) {
    printf("ERROR: charge of atom outside box in y-direction!\n");
    exit(0);
  }
  if (iatoz < 1 - nagz || iatoz > ngz + nagz) {
    printf("ERROR: charge of atom outside box in z-direction!\n");
    exit(0);
  }
  widthInv       = 1.0 / gw;
  fac            = 1.0 / pow(gw * sqrt(M_PI), 3);
  widthInvHhh[0] = widthInv * hx;
  widthInvHhh[1] = widthInv * hy;
  widthInvHhh[2] = widthInv * hz;
  widthInvXyz[0] = widthInv * xat;
  widthInvXyz[1] = widthInv * yat;
  widthInvXyz[2] = widthInv * zat;
  for (int ix = -nbgx; ix <= nbgx; ix++) {
    tt1            = widthInvHhh[0] * ix - widthInvXyz[0];
    wwx[ix + nbgx] = exp(-tt1 * tt1);
  }
  for (int iy = -nbgy; iy <= nbgy; iy++) {
    tt1            = widthInvHhh[1] * iy - widthInvXyz[1];
    wwy[iy + nbgy] = exp(-tt1 * tt1);
  }
  for (int iz = -nbgz; iz <= nbgz; iz++) {
    tt1            = widthInvHhh[2] * iz - widthInvXyz[2];
    wwz[iz + nbgz] = exp(-tt1 * tt1);
  }
#pragma omp parallel for schedule(static)
  for (int iz = -nbgz; iz <= nbgz; iz++) {
    double rhoz = fac * wwz[iz + nbgz];
    for (int iy = -nbgy; iy <= nbgy; iy++) {
      double rhoyz = rhoz * wwy[iy + nbgy];
      int ii       = nnx * nny * (iz + nbgz) + nnx * (iy + nbgy);
      for (int ix = -nbgx; ix <= nbgx; ix++) {
        wa[ii + ix + nbgx] += rhoyz * wwx[ix + nbgx];
      }
    }
  }
  memset(rho, 0, ngx * ngy * ngz * sizeof(double));
  chargeBackToCell(ngx, ngy, ngz, nagx, nagy, nagz, wa, rho);
  free(wwx);
  free(wwy);
  free(wwz);
  free(wa);
}
