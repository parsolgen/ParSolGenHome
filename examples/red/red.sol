import init_el(in int ii, in int jj, in int n, in int m, out Array<double> a[n][m]);
import print(in double value);

export sub main()
{
	int nb = 10;
	int mb = 10;
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

	double m = reduce(%max, A[:][:]);
	print(m);
}
