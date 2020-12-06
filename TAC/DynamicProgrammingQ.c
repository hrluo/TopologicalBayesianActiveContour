#include <math.h>
#include <stdlib.h>
#include "mex.h"
//Finding Minimum-Cost Path in a 2-D Matrix, w.r.t. equation (5.2.2) in J.D.Tucker's Thesis.
//Problem Statement : Given a cost matrix Cost[][] where Cost[i][j] denotes the Cost of visiting cell with coordinates (i,j), find a min-cost path to reach a cell (x,y) from cell (0,0) under the condition that you can only travel one step left/right or one step up/down. (We assume that all costs are positive integers)

// #define NNBRS	23

// const int Nbrs[NNBRS][2] = { 
// 	{ 1, 1 }, 
// 	{ 1, 2 },
// 	{ 2, 1 },
// 	{ 2, 3 },
// 	{ 3, 2 },
// 	{ 1, 3 },
// 	{ 3, 1 },
// 	{ 1, 4 },
// 	{ 3, 4 },
// 	{ 4, 3 },
// 	{ 4, 1 },
// 	{ 1, 5 },
// 	{ 2, 5 },
// 	{ 3, 5 },
// 	{ 4, 5 },
// 	{ 5, 4 },
// 	{ 5, 3 },
// 	{ 5, 2 },
// 	{ 5, 1 },
// 	{ 1, 6 },
// 	{ 5, 6 },
// 	{ 6, 5 },
// 	{ 6, 1 }
// };

#define NNBRS	63

const int Nbrs[NNBRS][2] = {
	{  1,  1 }, {  1,  2 }, {  1,  3 }, {  1,  4 }, {  1,  5 }, {  1,  6 }, {  1,  7 }, {  1,  8 }, {  1,  9 }, {  1, 10 },
	{  2,  1 }, {  2,  3 }, {  2,  5 }, {  2,  7 }, {  2,  9 }, {  3,  1 }, {  3,  2 }, {  3,  4 }, {  3,  5 }, {  3,  7 },
	{  3,  8 }, {  3, 10 }, {  4,  1 }, {  4,  3 }, {  4,  5 }, {  4,  7 }, {  4,  9 }, {  5,  1 }, {  5,  2 }, {  5,  3 },
	{  5,  4 }, {  5,  6 }, {  5,  7 }, {  5,  8 }, {  5,  9 }, {  6,  1 }, {  6,  5 }, {  6,  7 }, {  7,  1 }, {  7,  2 },
	{  7,  3 }, {  7,  4 }, {  7,  5 }, {  7,  6 }, {  7,  8 }, {  7,  9 }, {  7, 10 }, {  8,  1 }, {  8,  3 }, {  8,  5 },
	{  8,  7 }, {  8,  9 }, {  9,  1 }, {  9,  2 }, {  9,  4 }, {  9,  5 }, {  9,  7 }, {  9,  8 }, {  9, 10 }, { 10,  1 },
	{ 10,  3 }, { 10,  7 }, { 10,  9 }
};


int xycompare(const void *x1, const void *x2);
double CostFn2(const double *q1L, const double *q2L, int k, int l, int i, int j, int n, int scl);
void thomas(double *x, const double *a, const double *b, double *c, int n);
void spline(double *D, const double *y, int n);
void lookupspline(double *t, int *k, double dist, double len, int n);
double evalspline(double t, const double D[2], const double y[2]);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int i, j, k, l, n, M, N, Eidx, Fidx, Num, *Path, *xy, x, y, cnt;
	const int scl=1;
	const double *q1, *q2;
	double *q1L, *q2L, *yy, *D1, *D2, *tmp1, *tmp2, *E, Etmp, Emin, Ftmp, Fmin, t, a, b, lam = 0;

	if (nrhs != 4)
		mexErrMsgTxt("usage: [gam] = DynamicProgrammingQ(q1,q2,lam,Disp)");

	if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]))
		mexErrMsgTxt("Expected double precision arguments.");

	if (mxGetNumberOfDimensions(prhs[0]) != 2 || mxGetNumberOfDimensions(prhs[1]) != 2)
		mexErrMsgTxt("First two arguments expected to be two dimensional.");

	n = mxGetM(prhs[0]);//Number of rows in array, usually 2 in planar curve case.
	N = mxGetN(prhs[0]);//Number of cols in array, usually N in planar curve case.
    
	if (n != mxGetM(prhs[1]) || N != mxGetN(prhs[1]))
		mexErrMsgTxt("Dimension mismatch between first and second argument.");

	if (nlhs > 1)
		mexErrMsgTxt("Expected one return.");

	plhs[0] = mxCreateDoubleMatrix(1, N, mxREAL);

	yy = mxGetDoubles(plhs[0]);
	q1 = mxGetDoubles(prhs[0]);
	q2 = mxGetDoubles(prhs[1]);
	//Major loss of accuracy.
	//This a major loss of accuracy, MATLAB recommends mxGetDoubles instead of mxGetPr.
	M = scl*(N-1)+1;
	q1L = malloc(n*M*sizeof(double));
	q2L = malloc(n*M*sizeof(double));

	D1 = malloc(4*N*sizeof(double));
	tmp1 = D1 + N;
	D2 = D1 + 2*N;
	tmp2 = D2 + N;//tmp2=D1+3N, this is making sure that tmp1 and tmp2 are allocated to different spaces.

	//mexPrintf("Begin spline interp...\n");
	// compute spline interpolation 
	// for each dimension
	for (i = 0; i < n; ++i) {
		
		for (j = 0; j < N; ++j) {
			tmp1[j] = q1[n*j + i];
			tmp2[j] = q2[n*j + i];
		}

		//mexPrintf("Spline coeff for d=%d of q1\n",i);
		spline(D1, tmp1, N);
		//mexPrintf("Spline coeff for d=%d of q2\n",i);
		spline(D2, tmp2, N);

		// for each point in fine discretization
		for (j = 0; j < M; ++j) {
			//mexPrintf("Spline values at j=%d (out of M=%d)\n",j,M);
			lookupspline(&t, &k, j/(M-1.0), 1, N);
			q1L[n*j + i] = evalspline(t, D1+k, tmp1+k);
			q2L[n*j + i] = evalspline(t, D2+k, tmp2+k);
		}
	}

	free(D1);

	E = calloc(N*N, sizeof(double)); //E is the array that is of size N*N, which is the energy.
	Path = malloc(2*N*N*sizeof(int));//Path is the potential perturbation of the curve at N points.
	//Initialize the E matrix.

	for (i = 0; i < N; ++i) {
		E[N*i + 0] = 100000000000;
		E[N*0 + i] = 100000000000;
		Path[N*(N*0 + i) + 0] = -1;
		Path[N*(N*0 + 0) + i] = -1;
		Path[N*(N*1 + i) + 0] = -1;
		Path[N*(N*1 + 0) + i] = -1;
	}
	E[N*0 + 0] = 0;

	//mexPrintf("Begin DP...\n");
	for (j = 1; j < N; ++j) {
		for (i = 1; i < N; ++i) {

			Emin = 100000000000;
			Eidx = 0;

			for (Num = 0; Num < NNBRS; ++Num) {
				k = i - Nbrs[Num][0];
				l = j - Nbrs[Num][1];

				if (k >= 0 && l >= 0) {
					Etmp = E[N*l + k] + CostFn2(q1L,q2L,k,l,i,j,n,scl);
					if (Num == 0 || Etmp < Emin) {
						Emin = Etmp;
						Eidx = Num;
					}
				}
			}

			E[N*j + i] = Emin;
			Path[N*(N*0 + j) + i] = i - Nbrs[Eidx][0];
			Path[N*(N*1 + j) + i] = j - Nbrs[Eidx][1];
		}
	}

	free(E);
	free(q2L);

	xy = malloc(2*N*sizeof(int));
	xy[2*0 + 0] = N-1;
	xy[2*0 + 1] = N-1;

	cnt = 1;
	while (x = xy[2*(cnt-1) + 0], x > 0) {
		y = xy[2*(cnt-1) + 1];

		xy[2*cnt + 1] = Path[N*(N*0 + x) + y];
		xy[2*cnt + 0] = Path[N*(N*1 + x) + y];
		++cnt;
	}

	free(Path);

	qsort(xy, cnt, 2*sizeof(int), xycompare);

	for (i = 0; i < N; ++i) {

		Fmin = 100000000000;
		Fidx = 0;

		for (j = 0; j < cnt; ++j) {
			x = xy[2*j + 0];
			Ftmp = fabs(i - x);
			if (j == 0 || Ftmp < Fmin) {
				Fmin = Ftmp;
				Fidx = j;
			}
		}

		x = xy[2*Fidx + 0];
		y = xy[2*Fidx + 1];

		if (x == i) {
//			yy[i] = (y+1);
			yy[i] = y;
		}
		else {
			if (x > i) {
				a = x - i;
				b = i - xy[2*(Fidx-1) + 0];
				yy[i] = (a*xy[2*(Fidx-1) + 1] + b*y)/(a+b);
			}
			else {
				a = i - x;
				b = xy[2*(Fidx+1) + 0] - i;
				yy[i] = (a*xy[2*(Fidx+1) + 1] + b*y)/(a+b);
			}
		}

//		yy[i] /= N;
		yy[i] = (yy[i]-yy[0])/(N-1);
	}

	free(xy);
    free(q1L);
}

int xycompare(const void *x1, const void *x2) {
	return (*(int *)x1 > *(int *)x2) - (*(int *)x1 < *(int *)x2);
}

double CostFn2(const double *q1L, const double *q2L, int k, int l, int i, int j, int n, int scl) {
	double m = (j-l)/(double)(i-k), sqrtm = sqrt(m), E = 0, y, tmp, ip, fp;
	int x, idx, d, iL=i*scl, kL=k*scl, lL=l*scl;

	for (x = kL; x <= iL; ++x) {
		y = (x-kL)*m + lL;
		fp = modf(y, &ip);
		idx = (int)(ip + (fp >= 0.5));

		for (d = 0; d < n; ++d) {
			tmp = q1L[n*x + d] - sqrtm*q2L[n*idx + d];
			E += tmp*tmp;
		}
	}

	return E;
}

void thomas(double *x, const double *a, const double *b, double *c, int n) {
	//This function is the tridiagonal matrix algorithm, also known as the Thomas algorithm (named after Llewellyn Thomas), is a simplified form of Gaussian elimination that can be used to solve tridiagonal systems of equations. Here it is used for cubic spline.
	double tmp;
	int i;

	c[0] /= b[0];
	x[0] /= b[0];

	for (i = 1; i < n; ++i) {
		tmp = 1/(b[i] - c[i-1] * a[i]);
		c[i] *= tmp;
		x[i] = (x[i] - x[i-1] * a[i])*tmp;
	}

	for (i = n-2; i >= 0; --i) {
		x[i] -= c[i]*x[i+1];
	}
}

// input:  y is array to interpolate, n is array length
// output: D will be array of spline data
void spline(double *D, const double *y, int n) {
	int i;
	double *a, *b, *c;

	a = malloc(3*n*sizeof(double));
	b = a + n;
	c = b + n;

	if (n < 4) {
		a[0] = 0;
		b[0] = 2;
		c[0] = 1;
		D[0] = 3*(y[1]-y[0]);

		a[n-1] = 1;
		b[n-1] = 2;
		c[n-1] = 0;
		D[n-1] = 3*(y[n-1]-y[n-2]);
	}
	else {
		a[0] = 0;
		b[0] = 2;
		c[0] = 4;
		D[0] = -5*y[0] + 4*y[1] + y[2];

		a[n-1] = 4;
		b[n-1] = 2;
		c[n-1] = 0;
		D[n-1] = 5*y[n-1] - 4*y[n-2] - y[n-3];
	}

	for (i = 1; i < n-1; ++i) {
		a[i] = 1;
		b[i] = 4;
		c[i] = 1;
		D[i] = 3*(y[i+1]-y[i-1]);
	}

	thomas(D, a, b, c, n);

	free(a);
}

void lookupspline(double *t, int *k, double dist, double len, int n) {
	*t = (n-1)*dist/len;
	*k = (int)floor(*t);

	*k = (*k > 0)*(*k);
	*k += (*k > n-2)*(n-2-*k);

	*t -= *k;
}

double evalspline(double t, const double D[2], const double y[2]) {
	double c[4];

	c[0] = y[0];
	c[1] = D[0];
	c[2] = 3*(y[1]-y[0])-2*D[0]-D[1];
	c[3] = 2*(y[0]-y[1])+D[0]+D[1];

	return t*(t*(t*c[3] + c[2]) + c[1]) + c[0];
}
