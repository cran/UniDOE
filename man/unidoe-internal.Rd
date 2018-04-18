\name{Internal functions}
\title{UniDOE internal functions}
\alias{SATA_AUD}
\alias{SATA_UD}
\alias{SATA_LP}
\alias{CritEval}
\alias{panel.bar}
\alias{panel.scatter}
\alias{panel.heatmap}
\description{Internal core functions for constructing uniform design of experiments}
\usage{
SATA_UD(n,s,q,init,initX,crit,maxiter,hits_ratio)
SATA_AUD(XP, n, s, q, init, initX, crit, maxiter, hits_ratio)
SATA_LP(X0,q,crit,maxiter,hits_ratio)
CritEval(X0,q,crit)
panel.bar(x, ...)
panel.scatter(x, y, ...)
panel.heatmap(x, y, ...)
}
\author{
  Aijun Zhang, Haoyu Li, Shijie Quan
}
\details{These functions are internal functions, not directly recalled by R users.}
\keyword{internal}
