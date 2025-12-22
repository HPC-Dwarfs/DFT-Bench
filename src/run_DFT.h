/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of DFT-Bench.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#ifndef __RUN_DFT_H_
#define __RUN_DFT_H_

void test_put_gto_sym_ortho(double *rxyz,
    double gw,
    double rgcut,
    double *xyz111,
    int ngx,
    int ngy,
    int ngz,
    double *hgrid,
    double *orb);

void run_DFT(void);

#endif // __RUN_DFT_H_
