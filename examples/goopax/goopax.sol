import init_el(in int ii, in int jj, in int n, in int m, out Array<double> a[n][m]);
import print_el(in int ii, in int jj, in int n, in int m, in Array<double> a[n][m]);

import calc_main(in int n, in int m, in Array<double> prev[n][m], out Array<double> curr[n][m]) : @cpu, @goopax;

export sub main()
{
	int nb = 20;
	int mb = 20;
	int block_n = 10;
	int block_m = 10;

	Array<double[block_n][block_m]> A[nb][mb];

	for (int i = 0; i < nb; i++)
	{
		for (int j = 0; j < mb; ++j)
		{
			init_el(i, j, block_n, block_m, A[i][j]);
		}
	}
	for (int i = 0; i < nb; i++)
	{
		for (int j = 0; j < mb; ++j)
		{
			calc_main(block_n, block_m, A[i][j], A{1}[i][j]);
		}
	}
	for (int i = 0; i < nb; i++)
	{
		for (int j = 0; j < mb; ++j)
		{
			print_el(i, j, block_n, block_m, A{1}[i][j]);
		}
	}
}
