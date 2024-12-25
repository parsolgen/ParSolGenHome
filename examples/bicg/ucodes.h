#pragma once

#include "CgGate/rts.h"

#include <mkl_spblas.h>
#include <mpi.h>

#include <cmath>
#include <stdlib.h>
#include <vector>
#include <optional>
#include <fstream>

std::vector<std::optional<sparse_matrix_t>> mats;
std::vector<std::optional<sparse_matrix_t>> t_mats;

std::ifstream ifs;

int mat_n = -1;
int mat_nnz = -1;

int* mat_ia = nullptr;
int* mat_ja = nullptr;
double* mat_a = nullptr;

double t1 = 0.;

void init_mat_desc(int nb)
{
	mats.resize(nb * nb);
	t_mats.resize(nb * nb);
}

void read_matrix(const std::string& file_name)
{
   	std::ifstream ifs{file_name, std::ios::binary};
   	if (!ifs.is_open()) {
       std::cerr << "Error opening matrix file\n";
       rts::abort();
       return;
   	}
	ifs.read(reinterpret_cast<char*>(&mat_n), sizeof(int));
	ifs.read(reinterpret_cast<char*>(&mat_nnz), sizeof(int));

 	mat_ia = new int[mat_n+1];
	mat_ja = new int[mat_nnz];
	mat_a = new double[mat_nnz];

	ifs.read(reinterpret_cast<char*>(mat_ia), sizeof(int) * (mat_n + 1));
	ifs.read(reinterpret_cast<char*>(mat_ja), sizeof(int) * (mat_nnz));
	ifs.read(reinterpret_cast<char*>(mat_a), sizeof(double) * (mat_nnz));
}

void init_params(int nb, int& n, int& nnz, int& blockn)
{
	const auto args = rts::get_command_line_args();
	if (args.size() < 2) {
		std::cerr << "Bad arguments\n";
		rts::abort();
		return;
	}

	const auto file_path = args[1];
	read_matrix(file_path);

	if (mat_nnz == -1 || mat_n == -1) {
		std::cerr << "Incorrect matrix file\n";
		rts::abort();
		return;
	}

	n = mat_n;
	nnz = mat_nnz;
	blockn = n / nb;
}

void init_mat_arrays(int n, int nnz, int* ia, int* ja, double* a)
{
	std::copy(mat_ia, mat_ia + n + 1, ia);
	std::copy(mat_ja, mat_ja + nnz, ja);
	std::copy(mat_a, mat_a + nnz, a);
}

void print(int value)
{
	std::cout << value << "\n";
}

void print_d(double value)
{
	std::cout << "d = " << value << "\n";
}

void init_vec(int iblock, int blockn, double* block)
{
	for (int i = 0; i < blockn; ++i)
	{
		block[i] = 1;
	}
}

void init_zero_vec(int iblock, int blockn, double* block)
{
	for (int i = 0; i < blockn; ++i)
	{
		block[i] = 0;
	}
}

void csr_mv(
	int nb,
    int blockn,
	const rts::CSRMatrixBlock<int, double>& block,
	const double* x,
    const double* ysrc,
	double* ydst)
{
    if (!mats[nb * block.i + block.j]) {
      	sparse_matrix_t mat_block_desc;
      	mkl_sparse_d_create_csr(
            &mat_block_desc,
            SPARSE_INDEX_BASE_ZERO,
            block.n,
            block.m,
            block.ia,
            &block.ia[1],
            block.ja,
            block.a);
        mats[nb * block.i + block.j] = std::move(mat_block_desc);
    }
    matrix_descr desc;
    desc.type = SPARSE_MATRIX_TYPE_GENERAL;
    desc.diag = SPARSE_DIAG_NON_UNIT;
    mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1., mats[nb * block.i + block.j].value(), desc,
        x, 1., ydst);
}

void csr_mv_t(
	int nb,
    int blockn,
	const rts::CSRMatrixBlock<int, double>& block,
	const double* x,
    const double* ysrc,
	double* ydst)
{
    if (!t_mats[nb * block.i + block.j]) {
      	sparse_matrix_t mat_block_desc;
      	mkl_sparse_d_create_csr(
            &mat_block_desc,
            SPARSE_INDEX_BASE_ZERO,
            block.n,
            block.m,
            block.ia,
            &block.ia[1],
            block.ja,
            block.a);
        t_mats[nb * block.i + block.j] = std::move(mat_block_desc);
    }
    matrix_descr desc;
    desc.type = SPARSE_MATRIX_TYPE_GENERAL;
    desc.diag = SPARSE_DIAG_NON_UNIT;
    mkl_sparse_d_mv(SPARSE_OPERATION_CONJUGATE_TRANSPOSE, 1., t_mats[nb * block.i + block.j].value(), desc,
        x, 1., ydst);
}

void print_vec(int blocki, int blockn, const double* vec)
{
	std::cout << "i = " << blocki << ", val = [";
	for (int i = 0; i < blockn; ++i) {
		std::cout << vec[i] << ", ";
	}
	std::cout << "]\n";
}

void vec_sum(int blockn, const double* src1, const double* src2, double* dst)
{
	for (int i = 0; i < blockn; ++i) {
		dst[i] = src1[i] + src2[i];
	}
}

void sqrt(double a, double& b)
{
	b = std::sqrt(a);
}

void time_start()
{
	t1 = MPI_Wtime();
}

void print_time()
{
	std::cout << "Time taken: " << MPI_Wtime() - t1 << "\n";
}
