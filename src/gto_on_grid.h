/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of DFT-Bench.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#ifndef __GTO_ON_GRID_H_
#define __GTO_ON_GRID_H_

void put_gto_sym_ortho(char *bc,
    double *rxyz,
    double gw,
    double rgcut,
    double *xyz111,
    int ngx,
    int ngy,
    int ngz,
    double *hgrid,
    double *rho);

#endif // __GTO_ON_GRID_H_
