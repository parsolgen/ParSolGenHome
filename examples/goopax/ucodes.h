#pragma once

#include <iostream>
#include <goopax>

void print_el(int ii, int jj, std::size_t n, std::size_t m, const double* a)
{
	std::cout << "A[" << ii << "][" << jj << "]:\n";
	for (std::size_t i = 0; i < n; ++i)
	{
		for (std::size_t j = 0; j < m; ++j)
		{
			std::cout << a[i * m +j] << ' ';
		}
		std::cout << "\n";
	}
}

void init_el(int ii, int jj, std::size_t n, std::size_t m, double* a)
{
	for (std::size_t i = 0; i < n; ++i)
	{
		for (std::size_t j = 0; j < m; ++j)
		{
			a[i * m + j] = i + j;
		}
	}
}

void calc_main(int n, int m, const double* prev, double* curr)
{
	for (std::size_t i = 0; i < n; ++i)
	{
		for (std::size_t j = 0; j < m; ++j)
		{
			curr[i * m + j] = prev[i * m + j] * 2;
		}
	}
}

void calc_main_goopax(
	goopax::goopax_device dev,
	int n,
	int m, 
	const goopax::buffer<double>& prev,
	goopax::buffer<double>& curr)
{
   	goopax::kernel foo(dev, [](const goopax::resource<double>& A, goopax::resource<double>& B) {
   		B[0] = 1;
   	}, 0, 16);
   	foo(prev, curr);
}
