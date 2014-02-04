/**********************************************************************
betabinomial C implementation
***********************************************************************/

#include <stdio.h>
#include <math.h>

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


double lbeta(double alpha, double beta){
	return  lgamma(alpha) + lgamma(beta) - lgamma(alpha + beta);
}

// Beta-binomial probability
double dbb(int x, int n, double mu, double disp, int logp) {
	double y = mu * disp;
	double p = lbeta(x + y, n - x - y + disp) - lbeta(y, disp - y) + lgamma(n+1) - lgamma(x+1) -lgamma(n-x+1);
	if(! logp)
		p = exp(p);
	return p;
}

// Cumulative beta-binomial probability
double pbb(int x, int n, double mu, double disp, int logp) {
	int i, log=0;
	double p = 0;
	for (i=0; i<=x; i++){
		p += dbb(i, n, mu, disp, logp);
	}
	return p;
}

// Vectorized beta-binomial probability
int dbetabinom(double* p, int* lp, int *x, int* lx, int *n, int* ln, double *mu, int* lmu, double *disp, int* ldisp, int *logp) {
	int i, ix, in, imu, idisp;
	ix=in=imu=idisp=0;
	for(i=0; i<*lp;i++){
		p[i] = dbb(x[ix],n[in],mu[imu],disp[idisp],*logp);
		ix = (++ix==*lx) ? 0 : ix;
		in = (++in==*ln) ? 0 : in;
		imu = (++imu==*lmu) ? 0 : imu;
		idisp = (++idisp==*ldisp) ? 0 : idisp;
	}
	return 0;
}

// Vectorized cumulative beta-binomial probability
int pbetabinom(double* p, int*lp, int *x, int* lx, int *n, int* ln, double *mu, int* lmu, double *disp, int* ldisp, int *logp) {
	int i, ix, in, imu, idisp;
	ix=in=imu=idisp=0;
	for(i=0; i<*lp;i++){
		p[i] = pbb(x[ix],n[in],mu[imu],disp[idisp],*logp);
		ix = (++ix==*lx) ? 0 : ix;
		in = (++in==*ln) ? 0 : in;
		imu = (++imu==*lmu) ? 0 : imu;
		idisp = (++idisp==*ldisp) ? 0 : idisp;
	}
	return 0;
}

R_CMethodDef cMethods2[] = {
   {"dbetabinom", (DL_FUNC) &dbetabinom, 11},
   {"pbetabinom", (DL_FUNC) &pbetabinom, 11}
};

void R_init_betabinom(DllInfo *info) {
   R_registerRoutines(info, cMethods2, NULL, NULL, NULL);
}
