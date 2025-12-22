#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int max_int(int a, int b);

int modulo(int n, int m);

void charge_back_to_cell(int ngx,
    int ngy,
    int ngz,
    int nagx,
    int nagy,
    int nagz,
    int ibcx,
    double *wa,
    double *rho);

void put_gto_sym_ortho(char *bc,
    double *rxyz,
    double gw,
    double rgcut,
    double *xyz111,
    int ngx,
    int ngy,
    int ngz,
    double *hgrid,
    double *rho)
{
  //work array that is bigger than rho array, big enough to include of
  //grid points that are outside of box.
  //work arrays to save the values of one dimensional gaussian function.
  double *ww;
  double rhoz, rhoyz, pi, qat;
  double hgxinv, hgyinv, hgzinv;
  double width_inv, width_inv_xyz[3];
  double width_inv_hhh[3];
  double xat, yat, zat, facqiat, fac;
  double hx, hy, hz, tt1;
  int iatox, iatoy, iatoz, jx, jy, jz;
  int ii, nwa;
  int nbgx, nbgy, nbgz, nagx, nagy, nagz;
  int nbgmax;
  int ibcx, ibcy, ibcz, nnx, nny;
  double *wa;
  //int *mboundg;
  bool ortho;
  qat = 1.0;
  //ortho=true;
  //if(hgrid[1] > 0.0 || hgrid[1] < 0.0) ortho=false;
  //if(hgrid[2] > 0.0 || hgrid[2] < 0.0) ortho=false;
  //if(hgrid[3] > 0.0 || hgrid[3] < 0.0) ortho=false;
  //if(hgrid[5] > 0.0 || hgrid[5] < 0.0) ortho=false;
  //if(hgrid[6] > 0.0 || hgrid[6] < 0.0) ortho=false;
  //if(hgrid[7] > 0.0 || hgrid[7] < 0.0) ortho=false;
  ortho = (hgrid[1] > 0.0 || hgrid[1] < 0.0) ? false : true;
  ortho = (hgrid[2] > 0.0 || hgrid[2] < 0.0) ? false : ortho;
  ortho = (hgrid[3] > 0.0 || hgrid[3] < 0.0) ? false : ortho;
  ortho = (hgrid[5] > 0.0 || hgrid[5] < 0.0) ? false : ortho;
  ortho = (hgrid[6] > 0.0 || hgrid[6] < 0.0) ? false : ortho;
  ortho = (hgrid[7] > 0.0 || hgrid[7] < 0.0) ? false : ortho;
  if (!ortho) {
    printf("ERROR: this routine is only for orthogonal cell\n");
    printf("%24.15lf%24.15lf%24.15lf\n", hgrid[0], hgrid[1], hgrid[2]);
    printf("%24.15lf%24.15lf%24.15lf\n", hgrid[3], hgrid[4], hgrid[5]);
    printf("%24.15lf%24.15lf%24.15lf\n", hgrid[6], hgrid[7], hgrid[8]);
    exit(0);
  }
  hx   = hgrid[0];
  hy   = hgrid[4];
  hz   = hgrid[8];
  nbgx = (int)(rgcut / hx) + 2;
  nbgy = (int)(rgcut / hy) + 2;
  nbgz = (int)(rgcut / hz) + 2;
  nagx = nbgx + 1;
  nagy = nbgy + 1;
  nagz = nbgz + 1;
  ibcx = 0; //zero means periodic
  ibcy = 0; //zero means periodic
  ibcz = 0; //zero means periodic
  //printf("ngx,ngy,ngz = %d  %d  %d\n",ngx,ngy,ngz);
  //printf("nagx,nagy,nagz = %d  %d  %d\n",nagx,nagy,nagz);
  //mboundg=malloc((2*nbgy+1)*(2*nbgz+1)*2 * sizeof(int));
  //call get_glimitsphere(hx,hy,hz,nbgx,nbgy,nbgz,mboundg)
  nnx    = ngx + 2 * nagx;
  nny    = ngy + 2 * nagy;
  nwa    = nnx * nny * (ngz + 2 * nagz);
  wa     = malloc(nwa * sizeof(double));
  nbgmax = max_int(max_int(nbgx, nbgy), nbgz);
  ww     = malloc((2 * nbgmax + 1) * 3 * sizeof(double));
  pi     = 4.0 * atan(1.0);
  hgxinv = 1.0 / hx;
  hgyinv = 1.0 / hy;
  hgzinv = 1.0 / hz;
  for (int i = 0; i < nwa; i++)
    wa[i] = 0.0;
  //shift the gaussian centers
  iatox = (int)lround((rxyz[0] - xyz111[0]) * hgxinv) + 0;
  iatoy = (int)lround((rxyz[1] - xyz111[1]) * hgyinv) + 0;
  iatoz = (int)lround((rxyz[2] - xyz111[2]) * hgzinv) + 0;
  xat   = rxyz[0] - (iatox - 1) * hx - xyz111[0];
  yat   = rxyz[1] - (iatoy - 1) * hy - xyz111[1];
  zat   = rxyz[2] - (iatoz - 1) * hz - xyz111[2];
  if (iatox - nbgx * ibcx < 1 - nagx || iatox + nbgx * ibcx > ngx + nagx) {
    printf("ERROR: charge of atom outside box in x-direction!\n");
    exit(0);
  }
  if (iatoy - nbgy * ibcy < 1 - nagy || iatoy + nbgy * ibcy > ngy + nagy) {
    printf("ERROR: charge of atom outside box in y-direction!\n");
    exit(0);
  }
  if (iatoz - nbgz * ibcz < 1 - nagz || iatoz + nbgz * ibcz > ngz + nagz) {
    printf("ERROR: charge of atom outside box in z-direction!\n");
    exit(0);
  }
  //construct the one-dimensional gaussians
  width_inv        = 1.0 / gw;
  fac              = 1.0 / pow(gw * sqrt(pi), 3);
  width_inv_hhh[0] = width_inv * hx;
  width_inv_hhh[1] = width_inv * hy;
  width_inv_hhh[2] = width_inv * hz;
  width_inv_xyz[0] = width_inv * xat;
  width_inv_xyz[1] = width_inv * yat;
  width_inv_xyz[2] = width_inv * zat;
  for (int ii = 0; ii < 3; ii++) {
    for (int iww = -nbgmax; iww <= nbgmax; iww++) {
      tt1 = width_inv_hhh[ii] * iww - width_inv_xyz[ii];
      ww[(2 * nbgmax + 1) * ii + iww + nbgmax] = exp(-tt1 * tt1);
    }
  }
  facqiat = fac * qat;
  for (int iz = -nbgz; iz <= nbgz; iz++) {
    rhoz = facqiat * ww[(2 * nbgmax + 1) * 2 + iz + nbgmax];
    jz   = iatoz + iz;
    for (int iy = -nbgy; iy <= nbgy; iy++) {
      rhoyz = rhoz * ww[(2 * nbgmax + 1) + iy + nbgmax];
      jy    = iatoy + iy;
      //for(int ix=mboundg(1,iy,iz),mboundg(2,iy,iz)
      for (int ix = -nbgx; ix <= nbgx; ix++) {
        jx = iatox + ix;
        //write(*,'(5i5)') iat,iatox,ix,iatoy,iy
        ii                 = nnx * nny * (iz + nbgz) + nnx * (iy + nbgy);
        wa[ii + ix + nbgx] = wa[ii + ix + nbgx] + rhoyz * ww[ix + nbgmax];
      }
    }
  }
  //int kk;
  //for(int iz=nagz;iz<ngz+nagz;iz++)
  //{
  //  for(int iy=nagy;iy<ngy+nagy;iy++)
  //  {
  //    for(int ix=nagx;ix<ngx+nagx;ix++)
  //    {
  //      //kk= negx * negy * iz + negx * iy + ix;
  //      //kk      = nnx * nny * (nagz + iz + 0) + nnx * (nagy + iy + 0) + nagx + ix + 0;
  //      kk      = nnx * nny * iz + nnx * iy + ix;
  //      printf("AAA %5d %5d %5d %20.10lf\n",ix-nagx,iy-nagy,iz-nagz,wa[kk]);
  //    }
  //  }
  //}
  for (int i = 0; i < ngx * ngy * ngz; i++)
    rho[i] = 0.0;
  charge_back_to_cell(ngx, ngy, ngz, nagx, nagy, nagz, ibcx, wa, rho);
  //  !do iz=1,ngz
  //  !    do iy=1,ngy
  //  !        do ix=1,ngx
  //  !            write(61,'(3i4,es20.10)') ix,iy,iz,rho(ix,iy,iz)
  //  !        enddo
  //  !    enddo
  //  !enddo
  free(ww);
  free(wa);
  //  call f_free(mboundg)
}
void charge_back_to_cell(int ngx,
    int ngy,
    int ngz,
    int nagx,
    int nagy,
    int nagz,
    int ibcx,
    double *wa_i,
    double *rho)
{
  //real(8), intent(in):: wa(1-nagx:ngx+nagx,1-nagy:ngy+nagy,1-nagz:ngz+nagz)
  //real(8), intent(inout):: rho(ngx,ngy,ngz)
  //work arrays to save the values of one dimensional gaussian function.
  int jgx;
  int iix, iiy, iiz;
  int kk, ll;
  int ifinalx, igxs, igxf;
  int istartx, istarty, istartz;
  int nnx, nny;
  double *wa;
  nnx     = ngx + 2 * nagx;
  nny     = ngy + 2 * nagy;
  wa      = (wa_i + nnx * nny * (nagz - 1) + nnx * (nagy - 1) + nagx - 1);
  istartx = modulo(1 - nagx - 1, ngx) + 1;
  ifinalx = modulo(ngx + nagx - 1, ngx) + 1;
  igxs    = 1 - nagx + (ngx - istartx + 1);
  igxf    = ngx + nagx - ifinalx + 1 + ibcx * ngx;
  istarty = modulo(-nagy, ngy) + 1;
  istartz = modulo(-nagz, ngz) + 1;
  iiz     = istartz - 1;
  for (int igz = 1 - nagz; igz <= ngz + nagz; igz++) {
    iiy = istarty - 1;
    iiz = iiz + 1;
    if (iiz == ngz + 1)
      iiz = 1;
    for (int igy = 1 - nagy; igy <= ngy + nagy; igy++) {
      iiy = iiy + 1;
      if (iiy == ngy + 1)
        iiy = 1;
      iix = istartx - 1;
      for (int igx = 1 - nagx; igx <= igxs - 1; igx++) {
        iix     = iix + 1;
        kk      = ngx * ngy * (iiz - 1) + ngx * (iiy - 1) + iix - 1;
        ll      = nnx * nny * igz + nnx * igy + igx;
        rho[kk] = rho[kk] + wa[ll];
      }
      for (int igx = igxs; igx <= igxf - 1; igx += ngx) {
        for (int iix = 1; iix <= ngx; iix++) {
          jgx     = igx + iix - 1;
          kk      = ngx * ngy * (iiz - 1) + ngx * (iiy - 1) + iix - 1;
          ll      = nnx * nny * igz + nnx * igy + jgx;
          rho[kk] = rho[kk] + wa[ll];
        }
      }
      iix = 0;
      for (int igx = igxf; igx <= ngx + nagx; igx++) {
        iix     = iix + 1;
        kk      = ngx * ngy * (iiz - 1) + ngx * (iiy - 1) + iix - 1;
        ll      = nnx * nny * igz + nnx * igy + igx;
        rho[kk] = rho[kk] + wa[ll];
      }
    }
  }
}
int max_int(int a, int b)
{
  return (a > b) ? a : b;
}
int modulo(int n, int m)
{
  return ((n % m) + m) % m;
}
