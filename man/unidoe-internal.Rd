\name{UniDOE-internal}
\title{Internal UniDOE functions}
\alias{StoAUDC}
\alias{StoUDC}
\alias{StoLP}
\alias{StoEval}
\description{Internal UniDOE functions}
\usage{
StoUDC(n,s,q,init,initX,crit,maxiter,hits_ratio)
StoAUDC(XP, n, s, q, init, initX, crit, maxiter, hits_ratio)
StoLP(X0,q,crit,maxiter,hits_ratio)
StoEval(X0,q,crit)
}
\author{
  Aijun Zhang and Haoyu Li, Shijie Quan
}
\details{These are not intended for use by users. \code{StoUDC}Uniform Design Construction. \code{StoAUDC}Uniform Design Construction via row-wise augmentation. \code{StoLP}Level permutation for design construction. \code{StoEval}Evaluate discrepancy value for design.}
\keyword{internal}
