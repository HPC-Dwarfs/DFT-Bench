/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of DFT-Bench.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#ifndef __POSCARREADER_H_
#define __POSCARREADER_H_

typedef struct {
  int nat; /* number of atoms */
  double *rat; /* Cartesian atom positions in bohr, 3*nat */
  char (*sat)[5]; /* element symbol of each atom */
  double cell[3][3]; /* cell vectors (rows) in bohr */
} PoscarFileType;

/* Reads a VASP5 POSCAR file. Terminates the program on any error. */
void readPoscarFile(const char *filename, PoscarFileType *pf);
void freePoscarFile(PoscarFileType *pf);

#endif
