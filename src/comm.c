/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of DFT-Bench.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#include <limits.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "allocate.h"
#include "comm.h"
#include "constants.h"

#ifdef _MPI
#include <mpi.h>

#define MAX_WINDOWS 16
#define COMM_CHUNK ((size_t)1 << 30) /* elements per MPI call */

static MPI_Comm GroupComm = MPI_COMM_NULL;
static MPI_Comm LaneComm  = MPI_COMM_NULL;
static int GroupRank      = 0;
static MPI_Win Windows[MAX_WINDOWS];
static int NWindows = 0;

bool commIsMpi(void)
{
  return true;
}

void commInit(int *argc, char ***argv)
{
  MPI_Init(argc, argv);
}

void commFinalize(void)
{
  if (LaneComm != MPI_COMM_NULL) {
    MPI_Comm_free(&LaneComm);
  }
  if (GroupComm != MPI_COMM_NULL) {
    MPI_Comm_free(&GroupComm);
  }
  MPI_Finalize();
}

int commRank(void)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
}

int commSize(void)
{
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  return size;
}

void commSetupGroups(int maxGroupSize, CommType *comm)
{
  MPI_Comm node;
  int nrank, nsize;

  comm->rank = commRank();
  comm->size = commSize();
  MPI_Comm_split_type(
      MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, comm->rank, MPI_INFO_NULL, &node);
  MPI_Comm_rank(node, &nrank);
  MPI_Comm_size(node, &nsize);
  int color = (maxGroupSize > 0) ? nrank / maxGroupSize : 0;
  MPI_Comm_split(node, color, nrank, &GroupComm);
  MPI_Comm_free(&node);
  comm->nodeGroups = (maxGroupSize > 0) ? (nsize + maxGroupSize - 1) / maxGroupSize : 1;

  MPI_Comm_rank(GroupComm, &comm->grank);
  MPI_Comm_size(GroupComm, &comm->gsize);
  GroupRank = comm->grank;

  // lane 0 holds the group leaders, ordered by world rank: group 0 contains rank 0
  MPI_Comm_split(MPI_COMM_WORLD, comm->grank, comm->rank, &LaneComm);
  int ids[2] = { 0, 0 };
  if (comm->grank == 0) {
    MPI_Comm_rank(LaneComm, &ids[0]);
    MPI_Comm_size(LaneComm, &ids[1]);
  }
  MPI_Bcast(ids, 2, MPI_INT, 0, GroupComm);
  comm->group   = ids[0];
  comm->ngroups = ids[1];

  int mine[2]   = { comm->group, comm->grank };
  int *all      = malloc(2 * comm->size * sizeof(int));
  MPI_Allgather(mine, 2, MPI_INT, all, 2, MPI_INT, MPI_COMM_WORLD);
  comm->groupOf   = malloc(comm->size * sizeof(int));
  comm->grankOf   = malloc(comm->size * sizeof(int));
  comm->groupSize = calloc(comm->ngroups, sizeof(int));
  for (int p = 0; p < comm->size; p++) {
    comm->groupOf[p] = all[2 * p];
    comm->grankOf[p] = all[2 * p + 1];
    comm->groupSize[all[2 * p]]++;
  }
  free(all);
  comm->lmin = comm->size;
  for (int g = 0; g < comm->ngroups; g++) {
    comm->lmin = (comm->groupSize[g] < comm->lmin) ? comm->groupSize[g] : comm->lmin;
  }
}

void commPrintf(const char *fmt, ...)
{
  if (commRank() == 0) {
    va_list args;
    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);
  }
}

void dftAbort(const char *fmt, ...)
{
  if (commRank() == 0) {
    va_list args;
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);
    fflush(stderr);
  }
  fflush(stdout);
  MPI_Finalize();
  exit(EXIT_FAILURE);
}

void dftAbortRank(const char *fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  va_end(args);
  fflush(stderr);
  fflush(stdout);
  MPI_Abort(MPI_COMM_WORLD, 1);
}

static MPI_Comm scopeComm(CommScopeType scope)
{
  switch (scope) {
  case COMM_GROUP:
    return GroupComm;
  case COMM_LANE:
    return LaneComm;
  case COMM_LEADERS:
    if (GroupRank != 0) {
      dftAbortRank("ERROR: leader communication on a non-leader rank\n");
    }
    return LaneComm;
  default:
    return MPI_COMM_WORLD;
  }
}

static MPI_Op mpiOp(CommOpType op)
{
  switch (op) {
  case COMM_MAX:
    return MPI_MAX;
  case COMM_MIN:
    return MPI_MIN;
  default:
    return MPI_SUM;
  }
}

void commBarrier(CommScopeType scope)
{
  MPI_Barrier(scopeComm(scope));
}

void commAllreduce(double *buf, size_t n, CommOpType op, CommScopeType scope)
{
  MPI_Comm c = scopeComm(scope);
  for (size_t i = 0; i < n; i += COMM_CHUNK) {
    int m = (int)((n - i < COMM_CHUNK) ? n - i : COMM_CHUNK);
    MPI_Allreduce(MPI_IN_PLACE, buf + i, m, MPI_DOUBLE, mpiOp(op), c);
  }
}

void commReduceRoot(double *buf, size_t n, CommOpType op)
{
  bool root = commRank() == 0;
  for (size_t i = 0; i < n; i += COMM_CHUNK) {
    int m = (int)((n - i < COMM_CHUNK) ? n - i : COMM_CHUNK);
    MPI_Reduce(root ? MPI_IN_PLACE : buf + i,
        buf + i,
        m,
        MPI_DOUBLE,
        mpiOp(op),
        0,
        MPI_COMM_WORLD);
  }
}

bool commAllAnd(bool value)
{
  int v = value ? 1 : 0;
  MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
  return v != 0;
}

void commBcast(void *buf, size_t bytes, int root, CommScopeType scope)
{
  MPI_Comm c = scopeComm(scope);
  for (size_t i = 0; i < bytes; i += COMM_CHUNK) {
    int m = (int)((bytes - i < COMM_CHUNK) ? bytes - i : COMM_CHUNK);
    MPI_Bcast((char *)buf + i, m, MPI_BYTE, root, c);
  }
}

static int *toIntCounts(const size_t *v, int n)
{
  int *iv = malloc(n * sizeof(int));
  for (int i = 0; i < n; i++) {
    if (v[i] > INT_MAX) {
      dftAbortRank("ERROR: all-to-all count %zu exceeds the MPI int range\n", v[i]);
    }
    iv[i] = (int)v[i];
  }
  return iv;
}

void commAlltoallv(const double *sendbuf,
    const size_t *scounts,
    const size_t *sdispls,
    double *recvbuf,
    const size_t *rcounts,
    const size_t *rdispls)
{
  int size = commSize();
  int *sc  = toIntCounts(scounts, size);
  int *sd  = toIntCounts(sdispls, size);
  int *rc  = toIntCounts(rcounts, size);
  int *rd  = toIntCounts(rdispls, size);
  MPI_Alltoallv(sendbuf, sc, sd, MPI_DOUBLE, recvbuf, rc, rd, MPI_DOUBLE, MPI_COMM_WORLD);
  free(sc);
  free(sd);
  free(rc);
  free(rd);
}

void *commAllocShared(size_t bytes, CommWindowType *win)
{
  MPI_Win w;
  void *mine;
  int gsize;

  if (NWindows == MAX_WINDOWS) {
    dftAbortRank("ERROR: too many shared windows\n");
  }
  // MPI does not guarantee the alignment: group rank 0 allocates ARRAY_ALIGNMENT bytes
  // more and the group's memory starts at the first aligned address
  MPI_Info info;
  char align[16];
  snprintf(align, sizeof(align), "%d", ARRAY_ALIGNMENT);
  MPI_Info_create(&info);
  MPI_Info_set(info, "mpi_minimum_memory_alignment", align);
  size_t extra = (GroupRank == 0) ? ARRAY_ALIGNMENT : 0;
  MPI_Win_allocate_shared((MPI_Aint)(bytes + extra), 1, info, GroupComm, &mine, &w);
  MPI_Info_free(&info);
  MPI_Win_lock_all(MPI_MODE_NOCHECK, w);
  Windows[NWindows++]  = w;
  win->win             = malloc(sizeof(MPI_Win));
  *(MPI_Win *)win->win = w;

  // the segments must form one contiguous block in group rank order
  MPI_Comm_size(GroupComm, &gsize);
  char *start   = NULL;
  size_t offset = 0;
  for (int r = 0; r < gsize; r++) {
    MPI_Aint sz;
    int disp;
    char *ptr;
    MPI_Win_shared_query(w, r, &sz, &disp, &ptr);
    if (r == 0) {
      start = ptr;
    } else if (sz > 0 && ptr != start + offset) {
      dftAbortRank("ERROR: shared window segments are not contiguous\n");
    }
    offset += (size_t)sz;
  }
  return start + (ARRAY_ALIGNMENT - (uintptr_t)start % ARRAY_ALIGNMENT) % ARRAY_ALIGNMENT;
}

void commFreeShared(CommWindowType *win, void *ptr)
{
  MPI_Win w = *(MPI_Win *)win->win;
  (void)ptr;
  for (int i = 0; i < NWindows; i++) {
    if (Windows[i] == w) {
      Windows[i] = Windows[--NWindows];
      break;
    }
  }
  MPI_Win_unlock_all(w);
  MPI_Win_free(&w);
  free(win->win);
  win->win = NULL;
}

void commSyncGroup(void)
{
  for (int i = 0; i < NWindows; i++) {
    MPI_Win_sync(Windows[i]);
  }
  MPI_Barrier(GroupComm);
}

int commSubsetFortran(bool member)
{
  MPI_Comm sub;
  MPI_Comm_split(MPI_COMM_WORLD, member ? 0 : MPI_UNDEFINED, commRank(), &sub);
  return (sub == MPI_COMM_NULL) ? -1 : (int)MPI_Comm_c2f(sub);
}

void commFreeSubsetFortran(int handle)
{
  if (handle != -1) {
    MPI_Comm sub = MPI_Comm_f2c((MPI_Fint)handle);
    MPI_Comm_free(&sub);
  }
}

#else /* _MPI */

bool commIsMpi(void)
{
  return false;
}

void commInit(int *argc, char ***argv)
{
  (void)argc;
  (void)argv;
}

void commFinalize(void) { }

int commRank(void)
{
  return 0;
}

int commSize(void)
{
  return 1;
}

void commSetupGroups(int maxGroupSize, CommType *comm)
{
  (void)maxGroupSize;
  comm->rank = comm->group = comm->grank = 0;
  comm->size = comm->ngroups = comm->gsize = comm->lmin = comm->nodeGroups = 1;
  comm->groupOf      = malloc(sizeof(int));
  comm->grankOf      = malloc(sizeof(int));
  comm->groupSize    = malloc(sizeof(int));
  comm->groupOf[0]   = 0;
  comm->grankOf[0]   = 0;
  comm->groupSize[0] = 1;
}

void commPrintf(const char *fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  vprintf(fmt, args);
  va_end(args);
}

void dftAbort(const char *fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  va_end(args);
  exit(EXIT_FAILURE);
}

void dftAbortRank(const char *fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  va_end(args);
  exit(EXIT_FAILURE);
}

void commBarrier(CommScopeType scope)
{
  (void)scope;
}

void commAllreduce(double *buf, size_t n, CommOpType op, CommScopeType scope)
{
  (void)buf;
  (void)n;
  (void)op;
  (void)scope;
}

void commReduceRoot(double *buf, size_t n, CommOpType op)
{
  (void)buf;
  (void)n;
  (void)op;
}

bool commAllAnd(bool value)
{
  return value;
}

void commBcast(void *buf, size_t bytes, int root, CommScopeType scope)
{
  (void)buf;
  (void)bytes;
  (void)root;
  (void)scope;
}

void commAlltoallv(const double *sendbuf,
    const size_t *scounts,
    const size_t *sdispls,
    double *recvbuf,
    const size_t *rcounts,
    const size_t *rdispls)
{
  (void)rcounts;
  memcpy(recvbuf + rdispls[0], sendbuf + sdispls[0], scounts[0] * sizeof(double));
}

void *commAllocShared(size_t bytes, CommWindowType *win)
{
  win->win = NULL;
  return allocate(ARRAY_ALIGNMENT, bytes > 0 ? bytes : ARRAY_ALIGNMENT);
}

void commFreeShared(CommWindowType *win, void *ptr)
{
  (void)win;
  free(ptr);
}

void commSyncGroup(void) { }

int commSubsetFortran(bool member)
{
  return member ? 0 : -1;
}

void commFreeSubsetFortran(int handle)
{
  (void)handle;
}

#endif /* _MPI */

void commFreeGroups(CommType *comm)
{
  free(comm->groupOf);
  free(comm->grankOf);
  free(comm->groupSize);
}
