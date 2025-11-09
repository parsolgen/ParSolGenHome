#pragma once

#include <iostream>
#include <fstream>
#include <mpi.h>
#include <cassert>
#include <goopax>

double t = 0;
double total_exec_time = 0;

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

double goopax_time = 0.;

class MeasurementPointGoopax final
{
public:
	MeasurementPointGoopax()
	{
		time_ = MPI_Wtime();
	}
	~MeasurementPointGoopax()
	{
		goopax_time += (MPI_Wtime() - time_);
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

void calc_main_goopax(
	goopax::goopax_device dev,
	int in_n,
	int in_m,
	double a,
	double dt,
	double dx,
	double dy,
	const goopax::buffer<double>& prev,
	goopax::buffer<double>& curr)
{
	static constexpr const auto n_gpu_threads = 32;

	const auto n = in_n + 2;
	const auto m = in_m + 2;
	const auto dx2 = dx * dx;
	const auto dy2 = dy * dy;
	const auto nthreads = (n * m / n_gpu_threads + 1) * n_gpu_threads;

	MeasurementPointGoopax p1{};

	static goopax::kernel krnl(dev, [n, m, in_n, in_m, dx2, dy2, a, dt](
		const goopax::resource<double>& prev,
		goopax::resource<double>& curr)
	{
		const auto g_i = static_cast<goopax::gpu_int>(goopax::global_id());
		const auto i = g_i / n;
		const auto j = g_i % m;

		gpu_if(g_i < m)
		{
			curr[g_i] = prev[g_i];
			curr[(n - 1) * m + g_i] = prev[(n - 1) * m + g_i];
		}
		gpu_if(g_i < n)
		{
			curr[g_i * m] = prev[g_i * m];
			curr[g_i * m + m -1] = prev[g_i * m + m -1];
		}
		gpu_if(i >= 1 && i <= in_n)
		{
			gpu_if(j >= 1 && j <= in_m)
			{
				const auto ind = i * m + j;
				const auto ip = (i + 1) * m + j;
				const auto im = (i - 1) * m + j;
				const auto jp = i * m + j + 1;
				const auto jm = i * m + j - 1;
				curr[ind] = prev[ind] + a * dt *
			    	((prev[ip] - 2.0 * prev[ind] + prev[im]) / dx2 +
			    	(prev[jp] - 2.0 * prev[ind] + prev[jm]) / dy2);
			}
		}
	}, 0, nthreads);
	krnl(prev, curr).wait();
}

void print(double val)
{
	std::cout << "V = " << val << "\n";
}

void read_ghost(int in_n, int in_m, int ghost_edge_len, const double* a, double* ghost_edge)
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

void read_ghost_goopax(
	goopax::goopax_device dev,
	int in_n,
	int in_m,
	int ghost_edge_len,
	const goopax::buffer<double>& a,
	goopax::buffer<double>& ghost_edge)
{
	const auto n_gpu_threads = (in_n / 32 + 1) * 32;

	const auto n = in_n + 2;
	const auto m = in_m + 2;

	MeasurementPointGoopax p1{};

	static goopax::kernel krnl(dev, [n, m, in_n, in_m, ghost_edge_len](
		const goopax::resource<double>& a,
		goopax::resource<double>& ghost_edge)
	{
		const auto i = static_cast<goopax::gpu_int>(goopax::global_id());
		gpu_if(i >= 1 && i <= in_n)
		{
			//assert(i < ghost_edge_len);
			ghost_edge[i - 1] = a[i * m + 1];
		}
		// RIGHT
		gpu_if(i >= 1 && i <= in_n)
		{
			//assert(i < ghost_edge_len);
			ghost_edge[i + in_n - 1] = a[i * m + m - 2];
		}
		// UP
		gpu_if(i >= 1 && i <= in_m)
		{
			//assert(j < ghost_edge_len);
			ghost_edge[2 * in_n + i - 1] = a[m + i];
		}
		// DOWN
		gpu_if(i >= 1 && i <= in_m)
		{
			//assert(j < ghost_edge_len);
			ghost_edge[2 * in_n + in_m + i - 1] = a[(n - 2) * m + i];
		}
	}, 0, n_gpu_threads);
	krnl(a, ghost_edge).wait();
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

void write_ghost_goopax(
	goopax::goopax_device dev,
	int in_n,
	int in_m,
	int ghost_edge_len,
	const goopax::buffer<double>& left,
	const goopax::buffer<double>& right,
	const goopax::buffer<double>& up,
	const goopax::buffer<double>& down,
	const goopax::buffer<double>& a,
	goopax::buffer<double>& b)
{
	const auto n_gpu_threads = (in_n / 32 + 1) * 32;

	const auto n = in_n + 2;
	const auto m = in_m + 2;

	MeasurementPointGoopax p1{};

	static goopax::kernel krnl(dev, [n, m, in_n, in_m, ghost_edge_len](
		const goopax::resource<double>& left,
		const goopax::resource<double>& right,
		const goopax::resource<double>& up,
		const goopax::resource<double>& down,
		const goopax::resource<double>& a,
		goopax::resource<double>& b)
	{
		const auto i = static_cast<goopax::gpu_int>(goopax::global_id());
		// LEFT
		gpu_if (i < in_n)
		{
			// right -> left
			b[(i + 1) * m] = left[in_n + i];
		}
		// RIGHT
		gpu_if (i < in_n)
		{
			// left -> right
			b[(i + 1) * m + m - 1] = right[i];
		}
		// UP
		gpu_if (i < in_m)
		{
			// down -> up
			b[i + 1] = up[2 * in_n + in_m + i];
		}
		// DOWN
		gpu_if (i < in_m)
		{
			// up -> down
			b[(n - 1) * m + i + 1] = down[2 * in_n + i];
		}
		gpu_if (i >= 1 && i < in_n + 1) {
			goopax::gpu_for (1, in_m + 1, [m, i, &a, &b](goopax::gpu_int j) {
				b[i * m + j]  = a[i * m + j];
			});
		}
	}, 0, n_gpu_threads);
	krnl(left, right, up, down, a, b).wait();
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
	std::cout << "Total time taken: " << MPI_Wtime() - total_exec_time
		<< "\nCPU kernels took: " << t
		<< "\nGoopax kernels took: " << goopax_time << "\n";
}

void time_start()
{
	total_exec_time = MPI_Wtime();
}
