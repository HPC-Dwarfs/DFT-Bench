#include <stdio.h>
#include <stdlib.h>

void read_poscar(char *filename,int *nat_o,double **rat,char (**sat)[5],double *latvec_o);

void run_DFT(void) {
  char filename[]="poscars/POSCAR-0008-mp-149";
  int nat;
  double *rat=NULL;
  char (*sat)[5];
  double cellvec[9];
  read_poscar(filename,&nat,&rat,&sat,cellvec);
  free(rat);
  free(sat);
}
