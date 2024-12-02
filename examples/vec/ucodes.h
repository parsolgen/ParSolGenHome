#pragma once

#include <iostream>
#include <stdlib.h>


void print(int value)
{
	std::cout << value << "\n";
}

void print_dbl(double value)
{
	std::cout << value << "\n";
}

void init_vec(int iblock, int blockn, double* block)
{
	for (int i = 0; i < blockn; ++i)
	{
		block[i] = iblock;
	}
}

void init_zero_vec(int iblock, int blockn, double* block)
{
	for (int i = 0; i < blockn; ++i)
	{
		block[i] = 0;
	}
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

void hello_world()
{
	std::cout << "Hello world\n";
}

void set_var(int& a)
{
	a = 1;
}
