#pragma once

#include <mpi.h>
#include <cstdint>
#include <iostream>
#include <unistd.h>

extern "C" {
void dgemm_(const char* TRANSA, const char* TRANSB, const int* M, const int* N, const int* K, const double* alpha,
    const double* a, const int* ldA, const double* b, const int* ldB, const double* beta, double* c, const int* ldC);

void dsyrk_(const char* UPLO, const char* TRANS, const int* N, const int* K, const double* alpha, const double* a,
    const int* ldA, const double* beta, double* c, const int* ldC);

void dtrsm_(const char* SIDE, const char* UPLO, const char* TRANS, const char* DIAG, const int* M, const int* N,
    const double* alpha, const double* a, const int* ldA, double* b, const int* ldB);

void dpotrf_(const char* uplo, const int* n, double* a, const int* lda, int* info);
}

void hello()
{
	std::cout << "Works fine\n";
}

void init_block(
	int block_i,
	int block_j,
	int block_n,
	double* a)
{
    std::cout << "init(" << block_i << ", " << block_j << ")\n";
    for (std::size_t i = 0; i < block_n; ++i) {
        for (std::size_t j = 0; j < block_n; ++j) {
            if (block_i == block_j && i == j) {
                a[j * block_n + i] = 10e+10;
            }
            else {
                a[j * block_n + i] = i + j;
            }
        }
    }
}

void potrf(int n, const double* a_old, double* a)
{
    (void) a_old;
    int info = 0;
    const auto lda = n;

    std::cout << "potrf\n";
    dpotrf_("L", &n, a, &lda, &info);
}

void trsm(int n, const double* a, const double* in_a, double* b)
{
    const auto one = 1.;
    const auto lda = n;
    const auto ldb = n;

    dtrsm_("R", "L", "T", "N", &n, &n, &one, a, &lda, b, &ldb);
}

void syrk(int n, const double* a, const double* in_c, double* c)
{
    const auto minus_one = -1.;
    const auto lda = n;
    const auto ldc = n;
    const auto one = 1.;
    (void)in_c;

    dsyrk_("L", "N", &n, &n, &minus_one, a, &lda, &one, c, &ldc);
}

void gemm(int n, const double* a, const double* b, const double* in_c, double* c)
{
    const auto minus_one = -1.;
    const auto lda = n;
    const auto ldc = n;
    const auto ldb = n;
    const auto one = 1.;
    (void)in_c;

    dgemm_("N", "N", &n, &n, &n, &minus_one, a, &lda, b, &ldb, &one, c, &ldc);
}

void time(double& t)
{
    std::cout << "time()\n";
    t = MPI_Wtime();
}

void print_time(double t)
{
    std::cout << "Time taken: " << t << "\n";
}
