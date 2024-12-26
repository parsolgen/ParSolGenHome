import init_params(in int nb, out int n, out int nnz, out int blockn);
import print(in int a);
import print_d(in double a);
import init_mat_arrays(in int n, in int nnz, out Array<int> ia[n+1], out Array<int> ja[nnz], out Array<double> a[nnz]);
import init_vec(in int i, in int blockn, out Array<double> vec[blockn]);
import init_zero_vec(in int i, in int blockn, out Array<double> vec[blockn]);
import csr_mv(in int nb, in int n, in CSRBlock<double> A, in Array<double> x[n], in Array<double> ysrc[n], out Array<double> ydst[n]);
import print_vec(in int i, in int blockn, in Array<double> x[blockn]);
import vec_sum(in int blockn, in Array<double> x[blockn], in Array<double> x1[blockn], out Array<double> y[blockn]);
import init_mat_desc(in int nb);
import sqrt(in double a, out double b);
import csr_mv_t(in int nb, in int n, in CSRBlock<double> A, in Array<double> x[n], in Array<double> ysrc[n], out Array<double> ydst[n]);
import time_start();
import print_time();

sub main()
{
	int n = 0;
	int nnz = 0;
	int blockn = 0;
	int nb = 2;

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

	time_start();

	SingleVersionArray<double[blockn]> x[nb];
	SingleVersionArray<double[blockn]> B[nb];
	SingleVersionArray<double[blockn]> Rl[nb];
	SingleVersionArray<double[blockn]> Zl[nb];
	SingleVersionArray<double[blockn]> Pl[nb];
	SingleVersionArray<double[blockn]> Rr[nb];
	SingleVersionArray<double[blockn]> Zr[nb];
	SingleVersionArray<double[blockn]> Pr[nb];
	SingleVersionArray<double[blockn]> tmpr[nb];
	SingleVersionArray<double[blockn]> tmpl[nb];

	for (int i = 0; i < nb; ++i)
	{
		init_vec(i, blockn, B[i]);
		init_zero_vec(i, blockn, Rl[i]);
		init_zero_vec(i, blockn, Zl[i]);
		init_zero_vec(i, blockn, Pl[i]);
		init_zero_vec(i, blockn, Rr[i]);
		init_zero_vec(i, blockn, Zr[i]);
		init_zero_vec(i, blockn, Pr[i]);
		init_zero_vec(i, blockn, tmpr[i]);
		init_zero_vec(i, blockn, tmpl[i]);
	}

	Rr = B;
	Rl = Rr;
	Zr = Rr;
	Zl = Rl;
	Pr = Zr;
	Pl = Zl;

	double beta = (Zr, Rl);

	Zl = {};
	Zr = {};

	for (int i = 0; i < nb; ++i) {
		int v = 0;
		for (int j = 0; j < nb; ++j) {
			if (!%zero(A[i][j])) {
				csr_mv(nb, blockn, A[i][j], Pr[j], Zr{v}[i], Zr{v+1}[i]);
				v = v + 1;
			}
		}
	}

	for (int j = 0; j < nb; ++j) {
		int v = 0;
		for (int i = 0; i < nb; ++i) {
			if (!%zero(A[i][j])) {
				csr_mv_t(nb, blockn, A[i][j], Pl[i], Zl{v}[j], Zl{v+1}[j]);
				v = v + 1;
			}
		}
	}

	double dpi = (Rl, Zr{%latest});
	print_d(dpi);
	double a_c = beta / dpi;
	print_d(a_c);
	x{1} = x + a_c * Pr;
	double ma = a_c * -1.0;
	print_d(ma);
	Rr{1} = Rr + Zr{%latest} * ma;
	Rl{1} = Rl + Zl{%latest} * ma;
	double dp_t = (Rr{1}, Rr{1});
	double dp = 0;
	sqrt(dp_t, dp);
	print_d(dp);

	int it = 1;
	Zr = Rr{1};
	Zl = Rl{1};
	double betaold = beta;
	while (it < 100)
	{
		beta = (Zr, Rl);
		double b = beta / betaold;

		tmpr{it} = b * Pr{it};
		tmpl{it} = b * Pl{it};

		Pr{it + 1} = tmpr{it} + Zr{%latest};
		Pl{it + 1} = tmpl{it} + Zl{%latest};

		betaold = beta;

		Zr = {};
		Zl = {};

		for (int i = 0; i < nb; ++i) {
			int v = 0;
			for (int j = 0; j < nb; ++j) {
				if (!%zero(A[i][j])) {
					csr_mv(nb, blockn, A[i][j], Pr{it+1}[j], Zr{v}[i], Zr{v+1}[i]);
					v = v + 1;
				}
			}
		}

		for (int j = 0; j < nb; ++j) {
			int v = 0;
			for (int i = 0; i < nb; ++i) {
				if (!%zero(A[i][j])) {
					csr_mv_t(nb, blockn, A[i][j], Pl{it+1}[i], Zl{v}[j], Zl{v+1}[j]);
					v = v + 1;
				}
			}
		}

		double dpi = (Pl, Zr{%latest});
		double a_c = beta / dpi;

		x{it + 1} = x{it} + a_c * Pr{it + 1};
		double ma = a_c * -1.0;
		Rr{it + 1} = Rr{it} + Zr{%latest} * ma;
		Rl{it + 1} = Rl{it} + Zl{%latest} * ma;
		double dp_t = (Rr{it + 1}, Rr{it + 1});
		double dp = 0;
		sqrt(dp_t, dp);
		print_d(dp);

		it = it + 1;
	}

	print_time();
}