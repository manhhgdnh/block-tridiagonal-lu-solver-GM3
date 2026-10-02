# Block-Tridiagonal LU Solver

Numerical linear algebra project developed during GM3 at INSA Rouen Normandie. The project studies and implements LU factorization without pivoting for **tridiagonal** and **block-tridiagonal** linear systems, with an emphasis on numerical accuracy, conditioning and stability.

## Overview

The repository contains two related implementations:

- **Scalar tridiagonal systems**: specialized LU factorization followed by forward and backward substitution.
- **Block-tridiagonal systems**: block LU factorization using custom matrix utilities.

The project also includes representative input datasets, generated numerical results, a small Python helper for condition-number calculations, and the academic report.

## Numerical methods

For a tridiagonal or block-tridiagonal matrix \(A\), the programs compute a factorization

\[
A = LU,
\]

then solve

\[
Ly=b, \qquad Ux=y.
\]

The experiments are designed to study the effect of matrix size and conditioning on the numerical solution.

## Repository structure

```text
.
├── src/
│   ├── scalar/
│   │   ├── tridiag_lu_scalaire_etude.c
│   │   └── tridiag_lu_scalaire_general.c
│   └── block/
│       ├── tridiag_lu_blocs_etude.c
│       ├── tridiag_lu_blocs_general.c
│       └── matrices_tools.c
├── include/
│   └── matrices_tools.h
├── data/
│   ├── scalar/
│   └── block/
├── results/
│   ├── scalar/
│   └── block/
├── scripts/
│   └── condition_number.py
├── docs/
│   └── report.pdf
├── Makefile
├── .gitignore
└── README.md
```

## Build

Requirements:

- GCC
- GNU Make
- Python 3 with NumPy for the condition-number helper

Compile all C programs from the repository root:

```bash
make
```

This creates:

```text
bin/scalar_study
bin/scalar_general
bin/block_study
bin/block_general
```

Remove compiled binaries with:

```bash
make clean
```

## Run examples

### Scalar — study case

```bash
./bin/scalar_study \
  data/scalar/diffusion_n=6.d \
  results/scalar/scalar_study_output.txt
```

### Scalar — general case

```bash
./bin/scalar_general \
  data/scalar/diffusion_general_n=20.d \
  results/scalar/scalar_general_output.txt
```

### Block — study case

```bash
./bin/block_study \
  'data/block/laplace_N=3,n=3.d' \
  results/block/block_study_output.txt
```

### Block — general case

```bash
./bin/block_general \
  data/block/laplace_general_N=4_n=4.d \
  results/block/block_general_output.txt
```

### Condition-number helper

```bash
python3 scripts/condition_number.py
```

The current script demonstrates the condition-number calculation on a predefined tridiagonal matrix using NumPy.

## Main topics

- LU factorization
- Tridiagonal matrices
- Block-tridiagonal matrices
- Forward / backward substitution
- Numerical conditioning
- Stability and numerical error
- Scientific programming in C

## Authors

- Manh Hung Nguyen
- Tan Minh Duy Ngo

Academic project — GM3, INSA Rouen Normandie.
