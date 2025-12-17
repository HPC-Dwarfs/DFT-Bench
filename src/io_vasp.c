#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BUFFER_SIZE 256

void adjustl(char *str);

void read_poscar(
    char *filename, int *nat_o, double **rat_o, char (**sat_o)[5], double *cellvec)
{
  double bohr2ang = 0.529177210;
  int ntypat_tmp, natin, istat, nat_tmp, ityp;
  int *nitype;
  double *pos;
  double *rat;
  double scaling, cv[3][3];
  char all_line[256], all_line_tmp[256];
  char hlines[8][BUFFER_SIZE] = { 0 };
  bool reduced, selective_dynamics;
  int count, offset, chars_scanned;
  char *chptr_tmp = NULL;
  char (*char_type)[5], (*sat)[5];
  FILE *fptr;
  fptr = fopen(filename, "r");
  if (fptr != NULL) { }
  //first line is in principle a comment but usually contains the atomic elements
  for (int il = 0; il < 8; il++) {
    chptr_tmp = fgets(hlines[il], BUFFER_SIZE, fptr);
    if (chptr_tmp == NULL) {
      printf("ERROR: cannot read line %d of file %s\n", il + 1, filename);
    }
  }
  istat = sscanf(hlines[0], "%s", all_line_tmp);
  istat = sscanf(hlines[1], "%lf", &scaling);
  if (istat != 1)
    printf("ERROR: 2nd line of POSCAR must but one real number!");
  istat = sscanf(hlines[2], "%lf%lf%lf", &cellvec[0], &cellvec[1], &cellvec[2]);
  if (istat != 3)
    printf("ERROR: 3rd line of POSCAR must but three real numbers!");
  istat = sscanf(hlines[3], "%lf%lf%lf", &cellvec[3], &cellvec[4], &cellvec[5]);
  if (istat != 3)
    printf("ERROR: 4th line of POSCAR must but three real numbers!");
  istat = sscanf(hlines[4], "%lf%lf%lf", &cellvec[6], &cellvec[7], &cellvec[8]);
  if (istat != 3)
    printf("ERROR: 5th line of POSCAR must but three real numbers!");
  printf("%s\n", all_line_tmp);
  printf("%lf\n", scaling);
  printf("%21.16lf %21.16lf %21.16lf\n", cellvec[0], cellvec[1], cellvec[2]);
  printf("%21.16lf %21.16lf %21.16lf\n", cellvec[3], cellvec[4], cellvec[5]);
  printf("%21.16lf %21.16lf %21.16lf\n", cellvec[6], cellvec[7], cellvec[8]);
  strcpy(all_line, " ");
  strcat(all_line, hlines[5]);
  ntypat_tmp = 0;
  for (int i = 0; i < BUFFER_SIZE - 2; i++) {
    if (all_line[i] == ' ' && all_line[i + 1] != ' ') {
      ntypat_tmp = ntypat_tmp + 1;
    }
    if (all_line[i] == '\n')
      break;
  }
  nitype        = malloc(ntypat_tmp * sizeof(int));
  char_type     = malloc(ntypat_tmp * sizeof(*char_type));
  count         = 0;
  offset        = 0;
  chars_scanned = 0;
  while (count < 100 &&
         sscanf(hlines[5] + offset, "%19s%n", char_type[count], &chars_scanned) == 1) {
    offset += chars_scanned;
    count++;
  }
  printf("%s %s\n", char_type[0], char_type[1]);
  count         = 0;
  offset        = 0;
  chars_scanned = 0;
  while (count < 100 &&
         sscanf(hlines[6] + offset, "%19d%n", &nitype[count], &chars_scanned) == 1) {
    offset += chars_scanned;
    count++;
  }
  printf("%d %d\n", nitype[0], nitype[1]);
  //Compute total number of atoms and types
  natin = 0;
  for (int i = 0; i < ntypat_tmp; i++)
    natin = natin + nitype[i];
  pos    = malloc(3 * natin * sizeof(double));
  *rat_o = malloc(3 * natin * sizeof(double));
  *sat_o = malloc(natin * sizeof(char[5]));
  rat    = *rat_o;
  sat    = *sat_o;
  //Read positions
  strcpy(all_line, " ");
  strcpy(all_line, hlines[7]);
  selective_dynamics = false;
  adjustl(all_line);
  if (all_line[0] == 'S' || all_line[0] == 's')
    selective_dynamics = true;
  if (selective_dynamics) {
    chptr_tmp = fgets(all_line, BUFFER_SIZE, fptr);
  } else {
    strcpy(all_line, hlines[7]);
  }
  adjustl(all_line);
  if (all_line[0] == 'D' || all_line[0] == 'd') {
    reduced = true;
  } else if (all_line[0] == 'C' || all_line[0] == 'c') {
    reduced = false;
  } else {
    printf("ERROR: coordinates must be either Direct or Cartesian but %s\n", all_line);
  }
  printf("%s", all_line);
  ityp    = 0;
  nat_tmp = nitype[0];
  for (int iat = 0; iat < natin; iat++) {
    chptr_tmp = fgets(all_line, BUFFER_SIZE, fptr);
    if (chptr_tmp == NULL) {
      printf(
          "ERROR: cannot read coordinates of %d-th atom in file %s\n", iat + 1, filename);
    }
    istat = sscanf(
        all_line, "%lf%lf%lf", &pos[3 * iat], &pos[3 * iat + 1], &pos[3 * iat + 2]);
    if (iat == nat_tmp) {
      ityp++;
      nat_tmp += nitype[ityp];
    }
    for (int i = 0; i < 5; i++)
      sat[iat][i] = char_type[ityp][i];
  }
  for (int i = 0; i < 9; i++)
    cellvec[i] = cellvec[i] * (scaling / bohr2ang);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      cv[i][j] = cellvec[3 * i + j];
  if (reduced) {
    for (int iat = 0; iat < natin; iat++) {
      double x, y, z;
      x = pos[3 * iat + 0];
      y = pos[3 * iat + 1];
      z = pos[3 * iat + 2];
      printf("%21.16lf %21.16lf %21.16lf %s\n", x, y, z, sat[iat]);
      x = pos[3 * iat] * cv[0][0] + pos[3 * iat + 1] * cv[1][0] +
          pos[3 * iat + 2] * cv[2][0];
      y = pos[3 * iat] * cv[0][1] + pos[3 * iat + 1] * cv[1][1] +
          pos[3 * iat + 2] * cv[2][1];
      z = pos[3 * iat] * cv[0][2] + pos[3 * iat + 1] * cv[1][2] +
          pos[3 * iat + 2] * cv[2][2];
      rat[3 * iat + 0] = x;
      rat[3 * iat + 1] = y;
      rat[3 * iat + 2] = z;
    }
  } else {
    for (int iat = 0; iat < natin; iat++) {
      rat[3 * iat + 0] = pos[3 * iat + 0] / bohr2ang;
      rat[3 * iat + 1] = pos[3 * iat + 1] / bohr2ang;
      rat[3 * iat + 2] = pos[3 * iat + 2] / bohr2ang;
    }
  }
  //for(int iat=0;iat<natin;iat++) {
  //  double x, y, z;
  //  x=rat[0][3*iat+0]/cv[0][0];
  //  y=rat[0][3*iat+1]/cv[1][1];
  //  z=rat[0][3*iat+2]/cv[2][2];
  //  x=pos[3*iat+0];
  //  y=pos[3*iat+1];
  //  z=pos[3*iat+2];
  //  printf("%21.16lf %21.16lf %21.16lf\n",x,y,z);
  //}
  fclose(fptr);
  free(char_type);
  free(nitype);
  free(pos);
}
void adjustl(char *str)
{
  int nblanks;
  nblanks = 0;
  for (int i = 0; i < BUFFER_SIZE; i++) {
    if (str[i] == ' ')
      nblanks += 1;
    else
      break;
  }
  if (nblanks > 0) {
    for (int i = nblanks; i < BUFFER_SIZE; i++)
      str[i - nblanks] = str[i];
  }
}
