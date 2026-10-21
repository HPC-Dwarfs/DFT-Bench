/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of DFT-Bench.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#ifndef __COMM_H_
#define __COMM_H_

#include <stdbool.h>
#include <stddef.h>

/* Communication layer. All MPI calls of the benchmark live in comm.c. Without _MPI the
 * same interface is implemented for a single rank: collectives are no-ops or identities
 * and shared windows are plain allocations. */

typedef enum { COMM_WORLD, COMM_GROUP, COMM_LANE, COMM_LEADERS } CommScopeType;
typedef enum { COMM_SUM, COMM_MAX, COMM_MIN } CommOpType;

/* Ranks are organized in shared-memory groups (a node, or -s ranks of a node). Lane l
 * holds the ranks with group rank l of all groups, the leaders are lane 0. Group 0 is
 * the group of world rank 0. */
typedef struct {
  int rank, size; /* world */
  int group, ngroups; /* group of this rank, number of groups */
  int grank, gsize; /* rank and size within the group */
  int lmin; /* smallest group size */
  int nodeGroups; /* groups on the node of this rank */
  int *groupOf; /* group of every world rank */
  int *grankOf; /* group rank of every world rank */
  int *groupSize; /* ranks of every group */
} CommType;

typedef struct {
  void *win; /* MPI_Win of the MPI build, NULL otherwise */
} CommWindowType;

/* True in the MPI build */
bool commIsMpi(void);

void commInit(int *argc, char ***argv);
void commFinalize(void);
int commRank(void);
int commSize(void);

/* Builds groups of at most maxGroupSize ranks per node (0: whole node) */
void commSetupGroups(int maxGroupSize, CommType *comm);
void commFreeGroups(CommType *comm);

/* printf on world rank 0 only */
void commPrintf(const char *fmt, ...);

/* Terminates the job after an error that every rank detects alike: rank 0 prints the
 * message, all ranks exit with EXIT_FAILURE. */
void dftAbort(const char *fmt, ...);
/* Terminates the job after an error detected by the calling rank only */
void dftAbortRank(const char *fmt, ...);

void commBarrier(CommScopeType scope);
void commAllreduce(double *buf, size_t n, CommOpType op, CommScopeType scope);
/* In-place reduction over all ranks into world rank 0 */
void commReduceRoot(double *buf, size_t n, CommOpType op);
bool commAllAnd(bool value);
void commBcast(void *buf, size_t bytes, int root, CommScopeType scope);

/* All-to-all exchange of doubles over all ranks, counts and displacements in elements */
void commAlltoallv(const double *sendbuf,
    const size_t *scounts,
    const size_t *sdispls,
    double *recvbuf,
    const size_t *rcounts,
    const size_t *rdispls);

/* Allocates a window shared by the ranks of a group. Every rank contributes `bytes`,
 * the segments are contiguous in group rank order: rank r's memory starts after the
 * bytes of ranks 0..r-1. Returns the start of the group's memory, aligned to
 * ARRAY_ALIGNMENT. Rank boundaries are not page aligned, so a column-major block keeps
 * its stride across ranks. */
void *commAllocShared(size_t bytes, CommWindowType *win);
void commFreeShared(CommWindowType *win, void *ptr);
/* Makes writes to all shared windows visible within the group (sync + barrier) */
void commSyncGroup(void);

/* Communicator of the ranks with member = true, as a Fortran handle for libraries
 * such as ELPA (-1 on the other ranks). Called by all ranks. */
int commSubsetFortran(bool member);
void commFreeSubsetFortran(int handle);

#endif // __COMM_H_
