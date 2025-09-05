#pragma once

#include <iostream>

void init_el(int ii, int jj, std::size_t in_n, std::size_t in_m, double* a)
{
	for (int i = 0; i < in_n; ++i)
	{
		for (int j = 0; j < in_m; ++j)
		{
			a[i * in_m + j] = ii * i + jj * j;
		}
	}
}

void print(double val)
{
	std::cout << "V = " << val << "\n";
}
