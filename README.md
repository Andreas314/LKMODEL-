# LKMODEL

**LKMODEL** is a parallel C++ program for diagonalizing the Hamiltonian of the 8-band Luttinger–Kohn (LK) model for zinc-blend structures, a key framework in semiconductor physics. It computes the electronic band structure and the quantum interference (QI) tensor for detailed analysis of multiband quantum systems.

---

## Features

- Diagonalization of the 8-band LK Hamiltonian
- Parallel execution with **OpenMP**
- Output options:
  - `BAND_LKMODEL.dat` — Electronic band structure
  - `QI_LKMODEL.dat` — Quantum Interference tensor
- Includes:
  - `.Plot.gn` — Gnuplot script for band structure visualization
  - `Print_col.sh` — Shell script to print the N-th column of the QI tensor

---

## Requirements

To compile and run LKMODEL, ensure the following dependencies are installed:

- `g++` (supporting **C++17** or newer)
- `OpenMP`
- `LAPACK`
- `BLAS`

## Compilation
- Download with `git clone git@github.com:Andreas314/LKMODEL-.git`
- Compile with `make`
## Usage
- run `./lkmodel -h` to get the available options
- lkmodel reads `INPUT.txt` file, which contains information about the calculations
- see `INPUT_ex.txt` for explanation of the INPUT.txt file
- run `./Plot.gn` to plot the band structure from `BAND_LKMODEL.dat`
- run `./Print_col.sh` to show N-th column of the resulting tensor

