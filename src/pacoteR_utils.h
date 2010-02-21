#include <R.h>
#include <Rdefines.h>


double rpriors(int, double, double);
double dpriors(double, int, double, double, int);
double dinvwish(double *, double, double *);
void rinvwish(double, double *, double *);
void rmtvnorm(int, double *, double *, double *);

SEXP call_lsoda(SEXP, SEXP, SEXP, SEXP, SEXP rtol,
		SEXP, SEXP, SEXP, SEXP, SEXP,
		SEXP, SEXP, SEXP);
