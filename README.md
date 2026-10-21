# DFT-Bench

A planewave quantum chemistry inspired benchmark application.

DFT-Bench mimics the core operations of a planewave / real-space DFT code on a
set of trial orbitals for an atomistic structure read from a VASP POSCAR file:

1. **gto** – place two normalized Gaussian orbitals per atom on a periodic real-space
   grid (separable 1D Gaussians, folded back into the periodic cell).
2. **overlap** – compute the overlap matrix S = ΔV · Ψᵀ Ψ (BLAS `dsyrk`).
3. **eigen** – diagonalize S and build S^{-1/2} (LAPACK `dsyevd` + `dgemm`).
4. **rotate** – Löwdin orthogonalization Ψ ← Ψ · S^{-1/2} (blocked BLAS `dgemm`).
5. **kinetic** – kinetic energy ⟨ψ|−½∇²|ψ⟩ of all orbitals via a 3D real-to-complex
   FFT, multiplication by ½|k|² and a complex-to-real FFT (FFTW).

The grid spacing is fixed, so the number of grid points N grows linearly with the
number of atoms. The cost of the overlap and rotation phases is therefore O(N_atoms³),
as in production DFT codes.

This is the sequential version of the benchmark: it runs on a single core and
contains no threading.

## Getting Started

To build and run **DFT-Bench**, you need a C compiler, GNU make, FFTW3, and a
BLAS/LAPACK library with CBLAS and LAPACKE interfaces. Intel oneMKL provides both
BLAS/LAPACK and the FFTW3 interface.

1. **Install a supported compiler**
   - GCC
   - Clang
   - Intel ICX

2. **Clone the repository**

   ```sh
   git clone https://github.com/HPC-Dwarfs/DFT-Bench.git

   cd DFT-Bench
   ```

3. **Adjust configuration**

   The first `make` creates `config.mk` from `mk/config-default.mk`. Edit it to select
   the tool chain, the BLAS backend, etc.

4. **Build**

   ```sh
   make
   ```

5. **Run**

   ```sh
   ./dftbench-<TOOLCHAIN> -f poscars/POSCAR-0064-mp-149 -c
   ```

## Build

Configure the tool chain and additional options in `config.mk`:

```make
# Supported: GCC, CLANG, ICX
TOOLCHAIN ?= GCC
ENABLE_LIKWID ?= false

# BLAS/LAPACK and FFT backend. Supported: MKL, OPENBLAS
BLAS ?= MKL
# Installation prefixes for BLAS=OPENBLAS (leave empty for system paths)
OPENBLAS_ROOT ?=
FFTW_ROOT ?=

#Feature options
OPTIONS +=  -DARRAY_ALIGNMENT=64
#OPTIONS +=  -DVERBOSE_TIMER
```

- `BLAS=MKL` requires `MKLROOT` to be set (e.g. by sourcing oneAPI `setvars.sh`) and
  uses MKL's FFTW3 interface and its sequential layer.
- `BLAS=OPENBLAS` links `-lopenblas -lfftw3`. Any BLAS library that provides `cblas.h`
  and `lapacke.h` can be used this way, e.g. with a symlink or by adjusting
  `mk/include_BLAS.mk`. If the OpenBLAS library is multi-threaded, run with
  `OPENBLAS_NUM_THREADS=1`. Example for Homebrew on macOS:

  ```sh
  make TOOLCHAIN=CLANG BLAS=OPENBLAS \
       OPENBLAS_ROOT=/opt/homebrew/opt/openblas FFTW_ROOT=/opt/homebrew/opt/fftw
  ```

- `ENABLE_LIKWID=true` enables LIKWID marker regions named after the phases above
  (`gto`, `overlap`, `eigen`, `rotate`, `kinetic`):

  ```sh
  likwid-perfctr -C 0 -g FLOPS_DP -m ./dftbench-GCC -f poscars/POSCAR-0216-mp-149
  ```

- `-DVERBOSE_TIMER` prints the timer resolution.

You can build multiple tool chains in the same directory, but notice that the
Makefile is only acting on the one currently set. Intermediate build results are
located in the `./build/<TOOLCHAIN>` directory.

Clean up intermediate build results for the active tool chain with `make clean`, and
for all tool chains with `make distclean`. `make asm` generates assembler files in
`./build/<TOOLCHAIN>`, `make format` reformats all sources with `clang-format`.

## Usage

```
Usage: ./dftbench-<TOOLCHAIN> [options]
  -f <file>     POSCAR input file (default: poscars/POSCAR-0008-mp-149)
  -g <h>        target grid spacing in bohr (default: 0.30)
  -n <x,y,z>    explicit grid size, overrides -g
  -r <reps>     number of timed repetitions (default: 3)
  -c            run validation checks (exit code 1 on failure)
  -v            verbose output
  -h            print this help
```

The grid size per direction is the smallest even, FFT friendly (2, 3, 5 smooth) number
of points that achieves the requested spacing. Only orthogonal cells are supported.

The program prints the minimum and average time of every phase over all repetitions
and a nominal GFLOP/s rate. The rate is based on n(n+1)·N flops for the overlap,
2·N·n² for the rotation and 5·N·log₂N + 4·N per orbital for the kinetic energy, where
n is the number of orbitals and N the number of grid points.

### Validation

With `-c` the first repetition additionally checks:

- the placed orbitals against an independently computed periodic Gaussian,
- the orbital normalization,
- the numerical kinetic energy against the analytical value 3/(2·gw²) per orbital
  (fails if the grid is too coarse for the Gaussians),
- the orthonormality max|S − I| after the Löwdin orthogonalization.

The exit code is 1 if any check fails.

### Inputs

`poscars/` contains silicon diamond supercells (Materials Project mp-149). Grid size
and orbital memory at the default spacing of 0.3 bohr:

| Atoms | Orbitals | Grid  | Orbital memory |
| ----: | -------: | ----- | -------------: |
| 8     | 16       | 36³   | 0.01 GB        |
| 64    | 128      | 72³   | 0.36 GB        |
| 216   | 432      | 108³  | 4.05 GB        |
| 512   | 1024     | 144³  | 22.8 GB        |
| 1000  | 2000     | 180³  | 86.9 GB        |
| 1728  | 3456     | 216³  | 259 GB         |

The larger inputs (2744 to 8000 atoms) need several TB at this spacing. They are
intended for a future distributed-memory version or can be run with a coarser grid
(`-g`), which fails the kinetic-energy validation.
