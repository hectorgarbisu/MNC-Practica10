#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <mkl.h>
#include "MDF.h"

void getparams(const int N, const double *x, double *params){

	double *A = params;
	double *B = &(params[N]);
	double *C = &(params[2*N]);
	double *D = &(params[3*N]);

	for (int i = 0; i < N; i++){
		A[i] = x[i] * x[i];
		B[i] = x[i];
		C[i] = -1.0;
		D[i] = 3.0*x[i] * x[i];
	}
}

double exactSolution(double x){
	return x + 1.0 / x + x*x;
}

int main(int argc, char *argv[]){
	MDF *mdf = new MDF();
	double a = 0.5;
	double b = 2.0;
	int Nmax = 150;
	for (int N = 1; N < Nmax; N++){
		//int N = 3;
		double h = (b - a) / (double)(N - 1);
		double *x = (double *)mkl_malloc(N*sizeof(double), 64);
		double *y = (double *)mkl_malloc(N*sizeof(double), 64);
		double *params = (double *)mkl_malloc(4*N*sizeof(double), 64);

		for (int i = 0; i < N; i++) x[i] = a + (double)i*h;
		getparams(N, x, params);

		int ret = mdf->solve(N, h, y, params, CONTOUR_TYPE1, exactSolution(a), exactSolution(b));
		double *yexact = (double *)mkl_malloc(N*sizeof(double), 64);
		for (int i = 0; i < N; i++) yexact[i] = exactSolution(x[i]);

		double error = 0.0;
		for (int i = 0; i < N; i++) error += fabs((yexact[i] - y[i]) / yexact[i]);
		error /= N;
		printf("\nError relativo promedio(%%) para N=%d: %g\n", N, 100.0*error);
		mkl_free(x);
		mkl_free(y);
		mkl_free(yexact);
	}
	delete mdf;
	std::getchar();
	return 0;

}