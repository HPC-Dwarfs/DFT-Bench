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

#include "comm.h"
#include "constants.h"
#include "dft.h"
#include "likwid-marker.h"
#include "poscarReader.h"
#include "timing.h"
#include "util.h"

#ifdef _MPI
#define OPTSTRING "f:g:n:r:cvhs:m:e:p:"
#else
#define OPTSTRING "f:g:n:r:cvh"
#endif

static void usage(const char *prog)
{
  commPrintf("Usage: %s [options]\n"
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
#ifdef _MPI
  commPrintf("MPI options:\n"
             "  -s <n>        at most n ranks per shared-memory group (default: node)\n"
             "  -m <GB>       memory budget per group (default: %.0f%% of node memory)\n"
             "  -e <pr>x<pc>  eigensolver process grid (default: about %d^2 elements "
             "per rank)\n"
             "  -p <mode>     eigensolver rank placement: packed or spread (default: "
             "packed)\n",
      100.0 * DEFAULT_MEMORY_FRACTION,
      EIGEN_BLOCK_TARGET);
#endif
}

int main(int argc, char **argv)
{
  ConfigType cfg = {
    .poscarFile  = DEFAULT_POSCAR,
    .spacing     = DEFAULT_GRID_SPACING,
    .ng          = { 0, 0, 0 },
    .reps        = DEFAULT_REPETITIONS,
    .validate    = false,
    .verbose     = false,
    .groupSize   = 0,
    .memBudget   = 0.0,
    .eigenGrid   = { 0, 0 },
    .eigenSpread = false
  };
  int opt;

  commInit(&argc, &argv);
  // only rank 0 reports unknown options
  opterr = commRank() == 0;

  while ((opt = getopt(argc, argv, OPTSTRING)) != -1) {
    switch (opt) {
    case 'f':
      cfg.poscarFile = optarg;
      break;
    case 'g':
      cfg.spacing = atof(optarg);
      if (cfg.spacing <= 0.0) {
        dftAbort("ERROR: grid spacing must be positive\n");
      }
      break;
    case 'n':
      if (sscanf(optarg, "%d,%d,%d", &cfg.ng[0], &cfg.ng[1], &cfg.ng[2]) != 3 ||
          cfg.ng[0] < 2 || cfg.ng[1] < 2 || cfg.ng[2] < 2) {
        dftAbort("ERROR: -n expects three grid sizes >= 2, e.g. -n 32,32,48\n");
      }
      break;
    case 'r':
      cfg.reps = atoi(optarg);
      if (cfg.reps < 1) {
        dftAbort("ERROR: number of repetitions must be at least 1\n");
      }
      break;
    case 'c':
      cfg.validate = true;
      break;
    case 'v':
      cfg.verbose = true;
      break;
#ifdef _MPI
    case 's':
      cfg.groupSize = atoi(optarg);
      if (cfg.groupSize < 1) {
        dftAbort("ERROR: -s expects at least 1 rank per group\n");
      }
      break;
    case 'm':
      cfg.memBudget = atof(optarg);
      if (cfg.memBudget <= 0.0) {
        dftAbort("ERROR: -m expects a positive memory budget in GB\n");
      }
      break;
    case 'e': {
      char x;
      if (sscanf(optarg, "%d%c%d", &cfg.eigenGrid[0], &x, &cfg.eigenGrid[1]) != 3 ||
          (x != 'x' && x != 'X') || cfg.eigenGrid[0] < 1 || cfg.eigenGrid[1] < 1) {
        dftAbort("ERROR: -e expects a process grid <pr>x<pc>, e.g. -e 2x3\n");
      }
      break;
    }
    case 'p':
      if (strcmp(optarg, "packed") == 0) {
        cfg.eigenSpread = false;
      } else if (strcmp(optarg, "spread") == 0) {
        cfg.eigenSpread = true;
      } else {
        dftAbort("ERROR: -p expects packed or spread\n");
      }
      break;
#endif
    case 'h':
      usage(argv[0]);
      commFinalize();
      return EXIT_SUCCESS;
    default:
      usage(argv[0]);
      dftAbort("");
    }
  }

  CommType comm;
  commSetupGroups(cfg.groupSize, &comm);

  commPrintf(HLINE);
  commPrintf("DFT-Bench\n");
  commPrintf(HLINE);
#ifdef VERBOSE_TIMER
  commPrintf("Timer resolution %.2e s\n", getTimeResolution());
#endif

  LIKWID_MARKER_INIT;

  fflush(stdout);

  // fail once instead of on every rank if the input is missing
  FILE *fp = fopen(cfg.poscarFile, "re");
  if (fp != NULL) {
    fclose(fp);
  }
  if (!commAllAnd(fp != NULL)) {
    dftAbort("ERROR reading POSCAR file %s: cannot open file\n", cfg.poscarFile);
  }

  PoscarFileType pf;
  readPoscarFile(cfg.poscarFile, &pf);

  int status = runDft(&pf, &cfg, &comm);

  freePoscarFile(&pf);
  LIKWID_MARKER_CLOSE;
  commFreeGroups(&comm);
  commFinalize();

  return status;
}
