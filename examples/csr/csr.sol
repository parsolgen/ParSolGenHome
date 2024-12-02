import init_params(in int nb, out int n, out int nnz, out int blockn);
import print(in int a);
import init_mat_arrays(in int n, in int nnz, out Array<int> ia[n+1], out Array<int> ja[nnz], out Array<double> a[nnz]);
import init_vec(in int i, in int blockn, out Array<double> vec[blockn]);
import init_zero_vec(in int i, in int blockn, out Array<double> vec[blockn]);
import init_mat_desc(in int nb);
import csr_mv(in int nb, in int n, in CSRBlock<double> A, in Array<double> x[n], in Array<double> ysrc[n], out Array<double> ydst[n]);
import print_vec(in int i, in int blockn, in Array<double> x[blockn]);
import vec_sum(in int blockn, in Array<double> x[blockn], in Array<double> x1[blockn], out Array<double> y[blockn]);
			
export sub main()
{
	int n = 0;
	int nnz = 0;
	int blockn = 0;
	int nb = 19;

	init_params(nb, n, nnz, blockn);
	init_mat_desc(nb);
	print(n);
	print(nnz);
	print(nb);
	print(blockn);

	Array<int> ia[n+1];
	Array<int> ja[nnz];
	Array<double> a[nnz];

	init_mat_arrays(n, nnz, ia, ja, a);
	BlockedCSRMatrix<double[][]> A = {n, n, nnz, nb, nb, ia, ja, a};

	SingleVersionArray<double[blockn]> x[nb];
	SingleVersionArray<double[blockn]> y[nb];

	for (int i = 0; i < nb; ++i)
	{
		init_vec(i, blockn, x[i]);
		init_zero_vec(i, blockn, y[i]);
	}
	for (int i = 0; i < nb; ++i) {
		int v = 0;
		for (int j = 0; j < nb; ++j) {
			if (!%zero(A[i][j])) {
				csr_mv(nb, blockn, A[i][j], x[j], y{v}[i], y{v+1}[i]);
				v = v + 1;
			}
		}
	}
	for (int i = 0; i < nb; ++i) {
		print_vec(i, blockn, y{%latest}[i]);
	}
}
