#include <Rcpp.h>
using namespace Rcpp;

double StoEval(NumericMatrix X0, int q, int crit);
double MD2(NumericMatrix X0, int q);
double CL2(NumericMatrix X0, int q);
