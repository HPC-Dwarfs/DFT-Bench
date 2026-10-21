/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of DFT-Bench.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#ifndef __BLAS_LAPACK_H_
#define __BLAS_LAPACK_H_

/* Selects the CBLAS/LAPACKE headers of the configured backend (see BLAS in config.mk) */
#ifdef USE_MKL
#include <mkl_cblas.h>
#include <mkl_lapacke.h>
#else
#include <cblas.h>
#include <lapacke.h>
#endif

#endif // __BLAS_LAPACK_H_
