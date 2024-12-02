# ParSolGen (Parallel Solvers Generator) - main page

The repository contains ParSolGen alpha release and some examples: 
vectors sum (vec/), CSR sparse matrix-vector product (csr/) and 
LLT dense matrix decomposition algorithm (llt/).

## Install & build

To install ParSolGen (linux-based OS) simply extract the release 
archive. ParSolGen release consists of three directories:
1. `include/` contains headers necessary to compile parallel programs constructed by ParSolGen compiler
2. `bin/` contains ParSolGen compiler binary
3. `lib/` contains ParSolGen runtime library that should be linked to the constructed programs

The following dependencies should be pr-installed:
1. `mpich-4.2.3` 
2. `boost-1.79`
3. `libc++6`
4. `libc6`
5. `libm6`
6. `libomp`
7. `libpthread`
8. `libudev`
9. `GkLib`
10. `BLAS/LAPACK` implementation
11. `Intel MKL` (optional, necessary for building CSR examples)

## Build examples

1. `<path-to-parsolgen-bin>/sgc <example.sol>` - compile ParSolGen algorithm description into a .cpp file (`out.cpp` is the default name, to override use `-o` option)
2. `mpicxx out.cpp -L<path-to-parsolgen>/lib -I<path-to-parsolgen>/include -lPSRts -l<blas_lapack> -lGKlib -fopenmp`

The CSR-based examples should also be linked to IntelMKL Sparse BLAS library.

## Future plans

Distributed BiCG solver example is comming soon.