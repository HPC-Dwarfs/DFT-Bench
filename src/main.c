/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of DFT-Bench.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#include <limits.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#ifdef _OPENMP
#include "affinity.h"
#include <omp.h>
#endif

#include "constants.h"

int main(const int argc, char **argv) {
  const size_t bytesPerWord = sizeof(double);

#ifdef _OPENMP
  const size_t numThreads = omp_get_max_threads();
#else
  const size_t numThreads = 1;
#endif

#ifdef _OPENMP
  printf("OpenMP enabled, running with %zu threads\n", numThreads);

#ifdef VERBOSE_AFFINITY
#pragma omp parallel
  {
    int i = omp_get_thread_num();
#pragma omp critical
    {
      printf("Thread %d running on processor %d\n", i,
             affinity_getProcessorId());
      affinity_getmask();
    }
  }
#endif
#endif

  printf("Hi there, I'm running with %zu threads\n", numThreads);

  return EXIT_SUCCESS;
}
