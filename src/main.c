/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of DFT-Bench.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "constants.h"
#include "dft.h"
#include "likwid-marker.h"
#include "poscarReader.h"
#include "timing.h"
#include "util.h"

static void usage(const char *prog)
{
  printf("Usage: %s [options]\n"
         "  -f <file>     POSCAR input file (default: %s)\n"
         "  -g <h>        target grid spacing in bohr (default: %.2f)\n"
         "  -n <x,y,z>    explicit grid size, overrides -g\n"
         "  -r <reps>     number of timed repetitions (default: %d)\n"
         "  -c            run validation checks (exit code 1 on failure)\n"
         "  -v            verbose output\n"
         "  -h            print this help\n",
      prog,
      DEFAULT_POSCAR,
      DEFAULT_GRID_SPACING,
      DEFAULT_REPETITIONS);
}

int main(int argc, char **argv)
{
  ConfigType cfg = {
    .poscarFile = DEFAULT_POSCAR,
    .spacing    = DEFAULT_GRID_SPACING,
    .ng         = { 0, 0, 0 },
    .reps       = DEFAULT_REPETITIONS,
    .validate   = false,
    .verbose    = false
  };
  int opt;

  while ((opt = getopt(argc, argv, "f:g:n:r:cvh")) != -1) {
    switch (opt) {
    case 'f':
      cfg.poscarFile = optarg;
      break;
    case 'g':
      cfg.spacing = atof(optarg);
      if (cfg.spacing <= 0.0) {
        fprintf(stderr, "ERROR: grid spacing must be positive\n");
        return EXIT_FAILURE;
      }
      break;
    case 'n':
      if (sscanf(optarg, "%d,%d,%d", &cfg.ng[0], &cfg.ng[1], &cfg.ng[2]) != 3 ||
          cfg.ng[0] < 2 || cfg.ng[1] < 2 || cfg.ng[2] < 2) {
        fprintf(stderr, "ERROR: -n expects three grid sizes >= 2, e.g. -n 32,32,48\n");
        return EXIT_FAILURE;
      }
      break;
    case 'r':
      cfg.reps = atoi(optarg);
      if (cfg.reps < 1) {
        fprintf(stderr, "ERROR: number of repetitions must be at least 1\n");
        return EXIT_FAILURE;
      }
      break;
    case 'c':
      cfg.validate = true;
      break;
    case 'v':
      cfg.verbose = true;
      break;
    case 'h':
      usage(argv[0]);
      return EXIT_SUCCESS;
    default:
      usage(argv[0]);
      return EXIT_FAILURE;
    }
  }

  printf(HLINE);
  printf("DFT-Bench\n");
  printf(HLINE);
#ifdef VERBOSE_TIMER
  printf("Timer resolution %.2e s\n", getTimeResolution());
#endif

  LIKWID_MARKER_INIT;
  LIKWID_MARKER_REGISTER("gto");
  LIKWID_MARKER_REGISTER("overlap");
  LIKWID_MARKER_REGISTER("eigen");
  LIKWID_MARKER_REGISTER("rotate");
  LIKWID_MARKER_REGISTER("kinetic");

  fflush(stdout);

  PoscarFileType pf;
  readPoscarFile(cfg.poscarFile, &pf);

  int status = runDft(&pf, &cfg);

  freePoscarFile(&pf);
  LIKWID_MARKER_CLOSE;

  return status;
}
