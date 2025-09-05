#pragma once

#include <iostream>
#include <fstream>
#include <mpi.h>
#include <cassert>

double t = 0;

class MeasurementPoint final
{
public:
	MeasurementPoint()
	{
		time_ = MPI_Wtime();
	}
	~MeasurementPoint()
	{
		t += (MPI_Wtime() - time_);
	}

private:
	double time_ = 0;
};

void print_el(int ii, int jj, int in_n, int in_m, const double* a)
{
	MeasurementPoint p1{};

	const auto n = in_n + 2;
	const auto m = in_m + 2;
	std::ofstream of{"of_" + std::to_string(ii) + "_" + std::to_string(jj), std::ios::binary};
	of.write(reinterpret_cast<char*>(&in_n), sizeof(int));
	of.write(reinterpret_cast<char*>(&in_m), sizeof(int));
	for (int i = 1; i <= in_n; ++i)
	{
		for (int j = 1; j <= in_m; ++j)
		{
			of.write(reinterpret_cast<const char*>(&a[i * m + j]), sizeof(double));
		}
	}
}

void init_el(int ii, int jj, int in_n, int in_m, int nb, int mb, double* a)
{
	MeasurementPoint p1{};

	const auto n = in_n + 2;
	const auto m = in_m + 2;

	const double t_disc = 1.0;
	const double t_area = 65.0;

	const auto glob_n = in_n * nb;
	const auto glob_m = in_m * mb;

    const double radius = glob_n / 6.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
	    	const auto ind = i * m + j;
            /* Distance of point i, j from the origin */
            const auto glob_i = ii * n + i;
            const auto glob_j = jj * m + j;
            const auto dx = glob_i - glob_n / 2 + 1;
            const auto dy = glob_j - glob_m / 2 + 1;
            if (dx * dx + dy * dy < radius * radius) {
                a[ind] = t_disc;
            } else {
                a[ind] = t_area;
            }
        }
    }

}

void calc_main(
	int in_n,
	int in_m,
	double a,
	double dt,
	double dx,
	double dy,
	const double* prev,
	double* curr)
{
	MeasurementPoint p1{};

	const auto n = in_n + 2;
	const auto m = in_m + 2;
	const auto dx2 = dx * dx;
	const auto dy2 = dy * dy;

	for (std::size_t i = 0; i < m; ++i) {
		curr[i] = prev[i];
		curr[(n - 1) * m + i] = prev[(n - 1) * m + i];
	}
	for (std::size_t i = 0; i < n; ++i)
	{
		curr[i * m] = prev[i * m];
		curr[i * m + m -1] = prev[i * m + m -1];
	}
	for (int i = 1; i <= in_n; i++) {
		for (int j = 1; j <= in_m; j++) {
			int ind = i * m + j;
			int ip = (i + 1) * m + j;
			int im = (i - 1) * m + j;
			int jp = i * m + j + 1;
			int jm = i * m + j - 1;
			curr[ind] = prev[ind] + a * dt *
		    	((prev[ip] - 2.0 * prev[ind] + prev[im]) / dx2 +
		    	(prev[jp] - 2.0 * prev[ind] + prev[jm]) / dy2);
		}
	}
}

void print(double val)
{
	std::cout << "V = " << val << "\n";
}

void read_ghost(int in_n, int in_m, int ghost_edge_len, const double* a, auto* ghost_edge)
{
	MeasurementPoint p1{};

	// LEFT
	std::size_t g_i = 0;
	const auto n = in_n + 2;
	const auto m = in_m + 2;

	for (std::size_t i = 1; i <= in_n; ++i, ++g_i)
	{
		assert(i < ghost_edge_len);
		ghost_edge[g_i] = a[i * m + 1];
	}
	// RIGHT
	for (std::size_t i = 1; i <= in_n; ++i, ++g_i)
	{
		assert(i < ghost_edge_len);
		ghost_edge[g_i] = a[i * m + m - 2];
	}
	// UP
	for (std::size_t j = 1; j <= in_m; ++j, ++g_i)
	{
		assert(j < ghost_edge_len);
		ghost_edge[g_i] = a[m + j];
	}
	// DOWN
	for (std::size_t j = 1; j <= in_m; ++j, ++g_i)
	{
		assert(j < ghost_edge_len);
		ghost_edge[g_i] = a[(n - 2) * m + j];
	}
}

void write_ghost(
	int in_n,
	int in_m,
	int ghost_edge_len,
	const double* left,
	const double* right,
	const double* up,
	const double* down,
	const double* a,
	double* b)
{
	MeasurementPoint p1{};

	const auto n = in_n + 2;
	const auto m = in_m + 2;
	// LEFT
	for (std::size_t i = 0; i < in_n; ++i)
	{
		// right -> left
		b[(i + 1) * m] = left[in_n + i];
	}
	// RIGHT
	for (std::size_t i = 0; i < in_n; ++i)
	{
		// left -> right
		b[(i + 1) * m + m - 1] = right[i];
	}
	// UP
	for (std::size_t j = 0; j < in_m; ++j)
	{
		// down -> up
		b[j + 1] = up[2 * in_n + in_m + j];
	}
	// DOWN
	for (std::size_t j = 0; j < in_m; ++j)
	{
		// up -> down
		b[(n - 1) * m + j + 1] = down[2 * in_n + j];
	}
	for (std::size_t i = 1; i <= in_n; ++i) {
		for (std::size_t j = 1; j <= in_m; ++j) {
			b[i * m + j]  = a[i * m + j];
		}
	}
}

void print_iter(int iter)
{
	if (iter % 100 == 0)
	{
		std::cout << "IT: " << iter << "\n";
	}
}

void print_time()
{
	std::cout << "Calculation time: " << t << "\n";
}

