\name{UniTracePlot}
\alias{UniTracePlot}
\title{Draw discrepancy value trace in optimization process.}
\usage{
  UniTracePlot(output, skip=0)
}
\description{
  This function takes an output from GenUD or GenAUD function and a integer value skip to draw trace plot.
}
\arguments{
  \item{output}{an output from GenUD or GenAUD function.}

  \item{skip}{an integer value. Number of values in the head to skip for the given sequence.}
}
\value{
  A trace plot
}

\author{
  Aijun Zhang, Haoyu Li, Shijie Quan
}
\references{
  Zhang, A. and Li, H. (2017). UniDOE: An R package for constructing uniform design of experiments via stochastic and adaptive threshold accepting algorithm. Technical Report.
}


\examples{
  ##e.g.
  n=12 #(must be multiples of q)
  s=3
  q=4
  Dlst = GenUD(n=n,s=s,q=q,crit="MD2",maxiter=100)
  UniTracePlot(output=Dlst)
}
