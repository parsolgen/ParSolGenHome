#include <iostream>
#include <fstream>

#include <mpi.h>
#include <petsc.h>
#include <petscmat.h>

#include <algorithm>
#include <vector>

#include "mmio.h"

static char help[] = "Reads in a Symmetric matrix in MatrixMarket format. Writes\n\
it using the PETSc sparse format. It also adds a Vector set to random values to the\n\
output file. Input parameters are:\n\
  -fin <filename> : input file\n\
  -fout <filename> : output file\n\n";

struct CSRMAtrix final
{
    int* ia = nullptr;
    int* ja = nullptr;
    double* a = nullptr;
    int n = 0;
    int nnz = 0;

    CSRMAtrix() = delete;
    ~CSRMAtrix() {
        delete [] ia;
        delete [] ja;
        delete [] a;
    }
};

CSRMAtrix read_matrix_1(const char* file_name)
{
    int ret_code;
    MM_typecode matcode;
    FILE* f;

    int m, n, nnz;
    int i, *I, *J;
    double* val;
    if ((f = fopen(file_name, "r")) == NULL)
        exit(1);

    if (mm_read_banner(f, &matcode) != 0) {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && mm_is_sparse(matcode)) {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    if ((ret_code = mm_read_mtx_crd_size(f, &m, &n, &nnz)) != 0)
        exit(1);

    if (n != m) {
        std::cerr << "Matrix must be square\n";
        exit(1);
    }

    struct Coords final
    {
        int i;
        int j;
        double val;
    };

    std::vector<Coords> coords(nnz);

    val = new double[nnz];
    for (i = 0; i < nnz; i++) {
        int ii = 0, jj = 0;
        double val = 0.;
        auto r = fscanf(f, "%d %d %lg\n", &ii, &jj, &val);
        (void)r;
        ii--;
        jj--;
        coords[i] = Coords{ii, jj, val};
        /*if (I[i] != J[i] && val[i] == 0) {
            val[i] = 1.;
        }*/
    }
    std::sort(
        coords.begin(),
        coords.end(),
        [](const auto& lhs, const auto& rhs)
        {
            return lhs.i < rhs.i;
        });
    auto ia = new int[n + 1];
    auto ja = new int[nnz];
    for (int i = 0; i <= n; ++i) {
        ia[i] = 0;
    }
    for (int i = 0; i < nnz; i++) {
        ja[i] = coords[i].j;// J[i];
        ia[coords[i].i + 1]++;
        val[i] = coords[i].val;
    }
    for (int i = 0; i < n; i++) {
        ia[i + 1] += ia[i];
    }

    return CSRMAtrix{ia, ja, val, n, nnz};
}

int main(int argc, char **argv)
{
	Mat         A;
	Vec         b;
	char        filein[128],fileout[128];
	int ierr, size;

	PetscInitialize(&argc, &argv, (char *)0, help);

	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
	if (size > 1) {
		printf("too many processes\n");
		 exit(1);
	}

	/* Read in matrix and RHS */
	ierr = PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "-fin",filein,127,PETSC_NULL);CHKERRQ(ierr);

	auto matrix = read_matrix_1(filein);

	std::cout << "writing output file...\n";

	ierr = PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "-fout", fileout, 127, PETSC_NULL);CHKERRQ(ierr);

	std::ofstream os{fileout, std::ios::binary};
	const auto n = matrix.n;
	const auto nnz = matrix.nnz;
	os.write(reinterpret_cast<const char*>(&n), sizeof(int));
	os.write(reinterpret_cast<const char*>(&nnz), sizeof(int));
	std::cout << "nnz = " << nnz << "\n";
	os.write(reinterpret_cast<const char*>(matrix.ia), sizeof(int) * (n + 1));
	os.write(reinterpret_cast<const char*>(matrix.ja), sizeof(int) * (nnz));
	os.write(reinterpret_cast<const char*>(matrix.a), sizeof(double) * (nnz));

	ierr = PetscFinalize();CHKERRQ(ierr);

	return 0;
}