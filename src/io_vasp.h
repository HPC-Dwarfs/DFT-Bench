/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of DFT-Bench.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#ifndef __IO_VASP_H_
#define __IO_VASP_H_

void read_poscar(
    char *filename, int *nat_o, double **rat, char (**sat)[5], double *latvec_o);

#endif // __IO_VASP_H_
