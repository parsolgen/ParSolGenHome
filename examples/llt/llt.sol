import init_block(in int i, in int j, in int blockn, out Array<double> block[blockn][blockn]);
import potrf(in int n, in Array<double> a_in[n][n], out Array<double> a_out[n][n]);
import trsm(in int n, in Array<double> a_in[n][n], in Array<double> b_in[n][n], out Array<double> b[n][n]);
import syrk(in int n, in Array<double> a_in[n][n], in Array<double> c_in[n][n], out Array<double> c_out[n][n]);
import gemm(in int n, in Array<double> a[n][n], in Array<double> b[n][n], in Array<double> c_in[n][n], out Array<double> c_out[n][n]);
import hello();
import time(out double t);
import print_time(in double t);

export sub main()
{
        int n = 0;
        int blockn = 256;
        int nb = 100;
        int last_ver = nb-1;

        Matrix<double[blockn][blockn]> A[nb][nb];
        double t1 = 0;
        double t2 = 0;
        hello();
        
        for (int i = 0; i < nb; i++) {
                for (int j = 0; j < nb; ++j) {
                       init_block(i, j, blockn, A[i][j]); 
                }
        }
        time(t1);
        for (int i = 0; i < nb; i++) {
                potrf(blockn, A{i}[i][i], A{i+1}[i][i]);
        }
        for (int i = 0; i < nb; i++) {
                for (int j = i + 1; j < nb; j++) {
                        trsm(blockn, A{i+1}[i][i], A{i}[j][i], A{i+1}[j][i]);
                }
        }
        for (int v = 1; v < nb - 1; v++) {
                for (int i = v; i < nb; i++) {
                        syrk(blockn, A{v}[i][v-1], A{v-1}[i][i], A{v}[i][i]);
                        for (int i1 = i + 1; i1 < nb; i1++) {
                                gemm(blockn, A{v}[i1][v-1], A{v}[i][v-1], A{v-1}[i1][i], A{v}[i1][i]);
                        }
                }
        }
        syrk(blockn,
                A{last_ver}[last_ver][last_ver-1],
                A{last_ver-1}[last_ver][last_ver],
                A{last_ver}[last_ver][last_ver]);
        time(t2);
        print_time(t2-t1);
}
