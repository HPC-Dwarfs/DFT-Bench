/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of DFT-Bench.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "constants.h"
#include "poscarReader.h"

#define BUFFER_SIZE 256

static void adjustl(char *str)
{
  size_t skip = strspn(str, " \t");
  if (skip > 0)
    memmove(str, str + skip, strlen(str + skip) + 1);
}

static void poscarError(const char *filename, const char *fmt, ...)
{
  va_list args;
  fprintf(stderr, "ERROR reading POSCAR file %s: ", filename);
  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  va_end(args);
  fprintf(stderr, "\n");
  exit(EXIT_FAILURE);
}

static void *checkedMalloc(const char *filename, size_t bytes)
{
  void *ptr = malloc(bytes);
  if (ptr == NULL)
    poscarError(filename, "out of memory");
  return ptr;
}

void readPoscarFile(const char *filename, PoscarFileType *pf)
{
  int ntypatTmp, natin, istat, natTmp, ityp;
  int *nitype;
  double *pos;
  double scaling;
  char allLine[BUFFER_SIZE];
  char hlines[8][BUFFER_SIZE] = { 0 };
  bool reduced, selectiveDynamics;
  int count, offset, charsScanned;
  char (*charType)[5];
  FILE *fptr;
  fptr = fopen(filename, "re");
  if (fptr == NULL) {
    poscarError(filename, "cannot open file");
  }
  //first line is in principle a comment but usually contains the atomic elements
  for (int il = 0; il < 8; il++) {
    if (fgets(hlines[il], BUFFER_SIZE, fptr) == NULL) {
      poscarError(filename, "cannot read line %d", il + 1);
    }
  }
  istat = sscanf(hlines[1], "%lf", &scaling);
  if (istat != 1)
    poscarError(filename, "line 2 must contain one real number");
  for (int row = 0; row < 3; row++) {
    istat = sscanf(hlines[2 + row],
        "%lf%lf%lf",
        &pf->cell[row][0],
        &pf->cell[row][1],
        &pf->cell[row][2]);
    if (istat != 3)
      poscarError(filename, "line %d must contain three real numbers", row + 3);
  }
  // A negative scaling factor is the cell volume in VASP
  if (scaling < 0.0) {
    double (*c)[3] = pf->cell;
    double volume  = fabs(c[0][0] * (c[1][1] * c[2][2] - c[1][2] * c[2][1]) -
                          c[0][1] * (c[1][0] * c[2][2] - c[1][2] * c[2][0]) +
                          c[0][2] * (c[1][0] * c[2][1] - c[1][1] * c[2][0]));
    scaling        = cbrt(-scaling / volume);
  }
  // Count type names by scanning line 6
  ntypatTmp    = 0;
  offset       = 0;
  charsScanned = 0;
  char tmpWord[5];
  while (sscanf(hlines[5] + offset, "%4s%n", tmpWord, &charsScanned) == 1) {
    ntypatTmp++;
    offset += charsScanned;
  }
  if (ntypatTmp == 0)
    poscarError(filename, "line 6 must contain the element symbols (VASP5 format)");
  nitype       = checkedMalloc(filename, ntypatTmp * sizeof(int));
  charType     = checkedMalloc(filename, ntypatTmp * sizeof(*charType));
  count        = 0;
  offset       = 0;
  charsScanned = 0;
  while (count < ntypatTmp &&
         sscanf(hlines[5] + offset, "%4s%n", charType[count], &charsScanned) == 1) {
    offset += charsScanned;
    count++;
  }
  if (charType[0][0] >= '0' && charType[0][0] <= '9')
    poscarError(filename, "line 6 must contain the element symbols (VASP5 format)");
  count        = 0;
  offset       = 0;
  charsScanned = 0;
  while (count < ntypatTmp &&
         sscanf(hlines[6] + offset, "%d%n", &nitype[count], &charsScanned) == 1) {
    offset += charsScanned;
    count++;
  }
  if (count != ntypatTmp)
    poscarError(filename, "line 7 must contain %d atom counts", ntypatTmp);
  //Compute total number of atoms and types
  natin = 0;
  for (int i = 0; i < ntypatTmp; i++) {
    if (nitype[i] <= 0)
      poscarError(filename, "atom counts on line 7 must be positive");
    natin = natin + nitype[i];
  }
  pos     = checkedMalloc(filename, 3 * natin * sizeof(double));
  pf->rat = checkedMalloc(filename, 3 * natin * sizeof(double));
  pf->sat = checkedMalloc(filename, natin * sizeof(char[5]));
  pf->nat = natin;
  //Read positions
  strcpy(allLine, hlines[7]);
  adjustl(allLine);
  selectiveDynamics = (allLine[0] == 'S' || allLine[0] == 's');
  if (selectiveDynamics) {
    if (fgets(allLine, BUFFER_SIZE, fptr) == NULL)
      poscarError(filename, "cannot read coordinate mode line");
    adjustl(allLine);
  }
  if (allLine[0] == 'D' || allLine[0] == 'd') {
    reduced = true;
  } else if (allLine[0] == 'C' || allLine[0] == 'c' || allLine[0] == 'K' ||
             allLine[0] == 'k') {
    reduced = false;
  } else {
    poscarError(
        filename, "coordinates must be either Direct or Cartesian but are %s", allLine);
  }
  ityp   = 0;
  natTmp = nitype[0];
  for (int iat = 0; iat < natin; iat++) {
    if (fgets(allLine, BUFFER_SIZE, fptr) == NULL) {
      poscarError(filename, "cannot read coordinates of atom %d", iat + 1);
    }
    istat =
        sscanf(allLine, "%lf%lf%lf", &pos[3 * iat], &pos[3 * iat + 1], &pos[3 * iat + 2]);
    if (istat != 3)
      poscarError(filename, "coordinates of atom %d must be three real numbers", iat + 1);
    if (iat == natTmp) {
      ityp++;
      natTmp += nitype[ityp];
    }
    memcpy(pf->sat[iat], charType[ityp], sizeof(char[5]));
  }
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      pf->cell[i][j] *= scaling / BOHR2ANG;
    }
  }

  double (*cv)[3] = pf->cell;
  double *rat     = pf->rat;
  if (reduced) {
    for (int iat = 0; iat < natin; iat++) {
      double x, y, z;
      x                = pos[3 * iat] * cv[0][0] + pos[3 * iat + 1] * cv[1][0] +
                         pos[3 * iat + 2] * cv[2][0];
      y                = pos[3 * iat] * cv[0][1] + pos[3 * iat + 1] * cv[1][1] +
                         pos[3 * iat + 2] * cv[2][1];
      z                = pos[3 * iat] * cv[0][2] + pos[3 * iat + 1] * cv[1][2] +
                         pos[3 * iat + 2] * cv[2][2];
      rat[3 * iat + 0] = x;
      rat[3 * iat + 1] = y;
      rat[3 * iat + 2] = z;
    }
  } else {
    // Cartesian coordinates are scaled with the universal scaling factor as well
    for (int i = 0; i < 3 * natin; i++) {
      rat[i] = pos[i] * (scaling / BOHR2ANG);
    }
  }
  fclose(fptr);
  free(charType);
  free(nitype);
  free(pos);
}

void freePoscarFile(PoscarFileType *pf)
{
  free(pf->rat);
  free(pf->sat);
  pf->rat = NULL;
  pf->sat = NULL;
  pf->nat = 0;
}
