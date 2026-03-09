/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of DFT-Bench.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "constants.h"
#include "poscarReader.h"

#define BUFFER_SIZE 256

static void adjustl(char *str)
{
  size_t skip = strspn(str, " ");
  if (skip > 0)
    memmove(str, str + skip, strlen(str + skip) + 1);
}

void readPoscarFile(char *filename, PoscarFileType *pf)
{
  int ntypatTmp, natin, istat, natTmp, ityp;
  int *nitype;
  double *pos;
  double *rat;
  double scaling;
  char allLine[BUFFER_SIZE];
  char hlines[8][BUFFER_SIZE] = { 0 };
  bool reduced, selectiveDynamics;
  int count, offset, charsScanned;
  char *chptrTmp = NULL;
  char (*charType)[5], (*sat)[5];
  FILE *fptr;
  fptr = fopen(filename, "re");
  if (fptr == NULL) {
    printf("ERROR: cannot open file %s\n", filename);
    return;
  }
  //first line is in principle a comment but usually contains the atomic elements
  for (int il = 0; il < 8; il++) {
    chptrTmp = fgets(hlines[il], BUFFER_SIZE, fptr);
    if (chptrTmp == NULL) {
      printf("ERROR: cannot read line %d of file %s\n", il + 1, filename);
    }
  }
  istat = sscanf(hlines[1], "%lf", &scaling);
  if (istat != 1)
    printf("ERROR: 2nd line of POSCAR must but one real number!");
  for (int row = 0; row < 3; row++) {
    istat = sscanf(hlines[2 + row],
        "%lf%lf%lf",
        &pf->cellvec[3 * row],
        &pf->cellvec[3 * row + 1],
        &pf->cellvec[3 * row + 2]);
    if (istat != 3)
      printf("ERROR: line %d of POSCAR must have three real numbers!", row + 3);
  }
  // Count type names by scanning line 6
  ntypatTmp    = 0;
  count        = 0;
  offset       = 0;
  charsScanned = 0;
  char tmpWord[5];
  while (sscanf(hlines[5] + offset, "%4s%n", tmpWord, &charsScanned) == 1) {
    ntypatTmp++;
    offset += charsScanned;
  }
  nitype       = malloc(ntypatTmp * sizeof(int));
  charType     = malloc(ntypatTmp * sizeof(*charType));
  count        = 0;
  offset       = 0;
  charsScanned = 0;
  while (count < ntypatTmp &&
         sscanf(hlines[5] + offset, "%4s%n", charType[count], &charsScanned) == 1) {
    offset += charsScanned;
    count++;
  }
  count        = 0;
  offset       = 0;
  charsScanned = 0;
  while (count < ntypatTmp &&
         sscanf(hlines[6] + offset, "%d%n", &nitype[count], &charsScanned) == 1) {
    offset += charsScanned;
    count++;
  }
  //Compute total number of atoms and types
  natin = 0;
  for (int i = 0; i < ntypatTmp; i++)
    natin = natin + nitype[i];
  pos        = malloc(3 * natin * sizeof(double));
  *pf->rat_o = malloc(3 * natin * sizeof(double));
  *pf->sat_o = malloc(natin * sizeof(char[5]));
  *pf->nat_o = natin;
  rat        = *pf->rat_o;
  sat        = *pf->sat_o;
  //Read positions
  strcpy(allLine, hlines[7]);
  selectiveDynamics = false;
  adjustl(allLine);
  if (allLine[0] == 'S' || allLine[0] == 's')
    selectiveDynamics = true;
  if (selectiveDynamics) {
    chptrTmp = fgets(allLine, BUFFER_SIZE, fptr);
  } else {
    strcpy(allLine, hlines[7]);
  }
  adjustl(allLine);
  if (allLine[0] == 'D' || allLine[0] == 'd') {
    reduced = true;
  } else if (allLine[0] == 'C' || allLine[0] == 'c') {
    reduced = false;
  } else {
    printf("ERROR: coordinates must be either Direct or Cartesian but %s\n", allLine);
  }
  ityp   = 0;
  natTmp = nitype[0];
  for (int iat = 0; iat < natin; iat++) {
    chptrTmp = fgets(allLine, BUFFER_SIZE, fptr);
    if (chptrTmp == NULL) {
      printf(
          "ERROR: cannot read coordinates of %d-th atom in file %s\n", iat + 1, filename);
    }
    istat =
        sscanf(allLine, "%lf%lf%lf", &pos[3 * iat], &pos[3 * iat + 1], &pos[3 * iat + 2]);
    if (iat == natTmp) {
      ityp++;
      natTmp += nitype[ityp];
    }
    for (int i = 0; i < 5; i++) {
      sat[iat][i] = charType[ityp][i];
    }
  }
  for (int i = 0; i < 9; i++) {
    pf->cellvec[i] = pf->cellvec[i] * (scaling / BOHR2ANG);
  }

  double (*cv)[3] = (double (*)[3])pf->cellvec;
  if (reduced) {
    for (int iat = 0; iat < natin; iat++) {
      double x, y, z;
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
      rat[3 * iat + 0] = pos[3 * iat + 0] / BOHR2ANG;
      rat[3 * iat + 1] = pos[3 * iat + 1] / BOHR2ANG;
      rat[3 * iat + 2] = pos[3 * iat + 2] / BOHR2ANG;
    }
  }
  fclose(fptr);
  free(charType);
  free(nitype);
  free(pos);
}
