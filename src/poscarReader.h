/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of DFT-Bench.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#ifndef __POSCARREADER_H_
#define __POSCARREADER_H_

typedef struct {
  int *nat_o;
  double **rat_o;
  char (**sat_o)[5];
  double *cellvec;
} PoscarFileType;

void readPoscarFile(char *filename, PoscarFileType *pf);

#endif
