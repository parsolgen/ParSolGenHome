import init_vec(in int i, in int blockn, out Array<double> vec[blockn]);
import print_vec(in int i, in int blockn, in Array<double> x[blockn]);
import vec_sum(in int blockn,
	in Array<double> x[blockn],
	in Array<double> x1[blockn],
	out Array<double> y[blockn]);

export sub main()
{
	int blockn = 10;
	int nb = 5;

	Array<double[blockn]> a[nb];
	Array<double[blockn]> b[nb];
	Array<double[blockn]> c[nb];
	Array<double[blockn]> d[nb];

	for (int i = 0; i < nb; ++i)
	{
		init_vec(2, blockn, a[i]);
		init_vec(1, blockn, b[i]);
	}

    c = a + b;
    d = c * 2 - b;
        
	for (int i = 0; i < nb; ++i) {
		print_vec(i, blockn, d[i]);
	}
}
