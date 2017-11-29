#include "doe_Eval.h"
#include "doe_utility.h"
#include <string.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
double StoEval(NumericMatrix X0, int q, int crit=0)
{
  switch(crit)
  {
  case 0:
    return MD2(X0,q);
  case 1:
    return CL2(X0,q);
  default:
    return MD2(X0,q);
  }
}

double MD2(NumericMatrix X0, int q)
{
  int i,j,k,nsamp=X0.nrow(),nv=X0.ncol();
  double result=1,part1=0,part2=0,mul=1;
  NumericMatrix X(clone(X0));
  for(i=0;i<nsamp;i++) for(j=0;j<nv;j++) X(i,j) = (X(i,j)-0.5)/q;
  for(i=0;i<nv;i++) result *= 19/12.0;
  for(i=0;i<nsamp;i++)
  {
    mul=1;
    for(j=0;j<nv;j++) mul *= (5.0/3-0.25*ABS(X(i,j)-0.5)-0.25*ABS(X(i,j)-0.5)*ABS(X(i,j)-0.5));
    part1 += mul;
  }
  part1 *= (-2.0)/nsamp;
  for(i=0;i<nsamp;i++)
  {
    for(k=0;k<nsamp;k++)
    {
      mul=1;
      for(j=0;j<nv;j++) mul *=(15.0/8 - 0.25*ABS(X(i,j)-0.5) - 0.25*ABS(X(k,j)-0.5) -
          0.75*ABS(X(i,j)-X(k,j)) + 0.5*(X(i,j)-X(k,j))*(X(i,j)-X(k,j)));
      part2 += mul;
    }
  }
  part2 /= (nsamp*nsamp);
  result = result + part1 + part2;
  return result;
}

double CL2(NumericMatrix X0, int q)
{
  int i,j,k,nsamp=X0.nrow(),nv=X0.ncol();
  double result=1,part1=0,part2=0,mul=1;
  NumericMatrix X(clone(X0));
  for(i=0;i<nsamp;i++) for(j=0;j<nv;j++) X(i,j) = (X(i,j)-0.5)/q;
  for(i=0;i<nv;i++) result *= 13/12.0;
  for(i=0;i<nsamp;i++)
  {
    mul = 1.0;
    for(k=0; k<nv; k++) mul *= (1.0+0.5*ABS(X(i,k)-0.5)-0.5*(X(i,k)-0.5)*(X(i,k)-0.5));
    part1 += mul;
  }
  part1 *= (-2.0/nsamp);
  for(i=0;i<nsamp;i++) for(j=0;j<nsamp;j++)
  {
    mul = 1.0;
    for(k=0; k<nv; k++)mul *= (1.0+0.5*ABS(X(i,k)-0.5)+0.5*ABS(X(j,k)-0.5)-0.5*ABS(X(i,k)-X(j,k)));
    part2 += mul;
  }
  part2 *= (1.0/(nsamp*nsamp));
  result = result + part1 + part2;
  return result;
}
