import init_el(in int ii, in int jj, in int n, in int m, in int nb, in int mb, out Array<double> a[n][m]);
import read_ghost(in int n, in int m, in int ghost_edge_len, in Array<double> a[n][m], out Array<double> ghost_edge[ghost_edge_len][1]);
import write_ghost(in int n,
	in int m,
	in int ghost_edge_len,
	in Array<double> left[ghost_edge_len][1],
	in Array<double> right[ghost_edge_len][1],
	in Array<double> up[ghost_edge_len][1],
	in Array<double> down[ghost_edge_len][1],
	in Array<double> A[n][m],
	out Array<double> B[n][m]);
import print(in double value);
import print_el(in int i, in int j, in int n, in int m, in Array<double> a[n][m]);
import calc_main(
	in int n,
	in int m,
	in double a,
	in double dt,
	in double dx,
	in double dy,
	in Array<double> A[n][m],
	out Array<double> B[n][m]);
import print_iter(in int iter);
import print_time();

export sub main()
{
	int nb = 10;
	int mb = 10;
	int block_n = 100;
	int block_m = 100;
	int ghost_edge_len = block_n * 2 + block_m * 2;
	int it = 0;
	double a = 100;
	double dx = 0.01;
	double dy = 0.01;
    double dx2 = dx * dx;
    double dy2 = dy * dy;
    double dt = dx2 * dy2 / (2.0 * a * (dx2 + dy2));

	Array<double[block_n + 2][block_m + 2]> A[nb][mb];
	Array<double[block_n + 2][block_m + 2]> B[nb][mb];
	Array<double[ghost_edge_len][1]> ghost_edges[nb][mb];

	for (int i = 0; i < nb; i++)
	{
		for (int j = 0; j < mb; ++j)
		{
			init_el(i, j, block_n, block_m, nb, mb, A[i][j]);
		}
	}

	while (it < 5000)
	{
		for (int i = 0; i < nb; i++)
		{
			for (int j = 0; j < mb; ++j)
			{
				calc_main(block_n, block_m, a, dt, dx, dy, A{it}[i][j], A{it+1}[i][j]);
				read_ghost(block_n, block_m, ghost_edge_len, A{it+1}[i][j], ghost_edges{it}[i][j]);
			}
		}

		write_ghost(
			block_n,
			block_m,
			ghost_edge_len,
			ghost_edges{it}[0][0],
			ghost_edges{it}[0][1],
			ghost_edges{it}[0][0],
			ghost_edges{it}[1][0],
			A{it + 1}[0][0],
			A{it + 2}[0][0]);

		write_ghost(
			block_n,
			block_m,
			ghost_edge_len,
			ghost_edges{it}[0][mb-2],
			ghost_edges{it}[0][mb-1],
			ghost_edges{it}[0][mb-1],
			ghost_edges{it}[1][mb-1],
			A{it + 1}[0][mb-1],
			A{it + 2}[0][mb-1]);

		write_ghost(
			block_n,
			block_m,
			ghost_edge_len,
			ghost_edges{it}[nb-1][0],
			ghost_edges{it}[nb-1][1],
			ghost_edges{it}[nb-2][0],
			ghost_edges{it}[nb-1][0],
			A{it + 1}[nb-1][0],
			A{it + 2}[nb-1][0]);

		write_ghost(
			block_n,
			block_m,
			ghost_edge_len,
			ghost_edges{it}[nb-1][mb-2],
			ghost_edges{it}[nb-1][mb-1],
			ghost_edges{it}[nb-2][mb-1],
			ghost_edges{it}[nb-1][mb-1],
			A{it + 1}[nb-1][mb-1],
			A{it + 2}[nb-1][mb-1]);

		// Left
		for (int i = 1; i < nb - 1; i++)
		{
			write_ghost(
				block_n,
				block_m,
				ghost_edge_len,
				ghost_edges{it}[i][0],
				ghost_edges{it}[i][1],
				ghost_edges{it}[i-1][0],
				ghost_edges{it}[i+1][0],
				A{it + 1}[i][0],
				A{it + 2}[i][0]);
		}
		// Right
		for (int i = 1; i < nb - 1; i++)
		{
			write_ghost(
				block_n,
				block_m,
				ghost_edge_len,
				ghost_edges{it}[i][mb-2],
				ghost_edges{it}[i][mb-1],
				ghost_edges{it}[i-1][mb-1],
				ghost_edges{it}[i+1][mb-1],
				A{it + 1}[i][mb-1],
				A{it + 2}[i][mb-1]);
		}
		// Up
		for (int j = 1; j < mb - 1; j++)
		{
			write_ghost(
				block_n,
				block_m,
				ghost_edge_len,
				ghost_edges{it}[0][j-1],
				ghost_edges{it}[0][j+1],
				ghost_edges{it}[0][j],
				ghost_edges{it}[1][j],
				A{it + 1}[0][j],
				A{it + 2}[0][j]);
		}
		// Down
		for (int j = 1; j < mb - 1; j++)
		{
			write_ghost(
				block_n,
				block_m,
				ghost_edge_len,
				ghost_edges{it}[nb-1][j-1],
				ghost_edges{it}[nb-1][j+1],
				ghost_edges{it}[nb-2][j],
				ghost_edges{it}[nb-1][j],
				A{it + 1}[nb-1][j],
				A{it + 2}[nb-1][j]);
		}
		// Center
		for (int i = 1; i < nb - 1; i++)
		{
			for (int j = 1; j < mb -1; j++)
			{
				write_ghost(
					block_n,
					block_m,
					ghost_edge_len,
					ghost_edges{it}[i][j-1],
					ghost_edges{it}[i][j+1],
					ghost_edges{it}[i - 1][j],
					ghost_edges{it}[i + 1][j],
					A{it + 1}[i][j],
					A{it + 2}[i][j]);
			}
		}
		print_iter(it);
		it = it + 2;
	}

	for (int i = 0; i < nb; i++)
	{
		for (int j = 0; j < mb; j++)
		{
			print_el(i, j, block_n, block_m, A{it}[i][j]);
		}
	}

	print_time();
}
