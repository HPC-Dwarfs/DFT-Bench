# DFT-Bench

A planewave quantum chemistry inspired benchmark application.

## Getting Started

To build and run **DFT-Bench**, you only need a compiler and GNU
make.

1. **Install a supported compiler**
   - GCC
   - Clang
   - Intel ICC/ICX

2. **Clone the repository**

   ```sh
   git clone https://github.com/HPC-Dwarfs/DFT-Bench.git

   cd DFT-Bench
   ```

3. **(Optional) Adjust configuration**

   Edit `config.mk` to change the default problem size, enable OpenMP, etc.

4. **Build**

   ```sh
   make
   ```

   See the full [Build](#build) section for more details.

5. **Usage**

   ```sh
   ./dftbench-<TOOLCHAIN>
   ```

   See the full [Usage](#usage) section for more details.

   Get _Help_ on command line arguments:

   ```sh
   ./bwBench-<TOOLCHAIN> -h
   ```

## Build

### CPU Build

1. Configure the tool chain and additional options in `config.mk`:

```make
# Supported: GCC, CLANG, ICC, ICX
TOOLCHAIN ?= GCC
ENABLE_OPENMP ?= true
ENABLE_LIKWID ?= false

#Feature options
OPTIONS +=  -DARRAY_ALIGNMENT=64
#OPTIONS +=  -DVERBOSE_AFFINITY
#OPTIONS +=  -DVERBOSE_DATASIZE
#OPTIONS +=  -DVERBOSE_TIMER
```

The verbosity options enable detailed output about affinity settings, allocation
sizes, and timer resolution. If you uncomment `-DVERBOSE_AFFINITY` the processor
id every thread is currently scheduled on and the complete affinity mask for
every thread is printed.

- Build with:

```sh
make
```

You can build multiple tool chains in the same directory, but notice that the
Makefile is only acting on the one currently set. Intermediate build results are
located in the `./build/<TOOLCHAIN>` directory.

- Clean up intermediate build results for active tool chain, data files and plots with:

```sh
make clean
```

Clean all build results for all tool chains:

```sh
make distclean
```

- Optional targets:

Generate assembler:

```sh
make asm
```

The assembler files will also be located in the `./build/<TOOLCHAIN>` directory.

Reformat all source files using `clang-format` (only works if `clang-format` is
in your path):

```sh
make format
```
