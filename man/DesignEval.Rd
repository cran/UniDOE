\name{DesignEval}
\alias{DesignEval}
\title{Evaluate design in terms of criteria}
\usage{
DesignEval(X0,crit)
}
\description{
  This function takes matrix X0,q and crit to output the criterion value.
}
\arguments{
\item{X0}{an integer matrix object}

\item{crit}{an integer object, criterion to choose:

              "MD2" -- MD2

              "CD2" -- CD2

              "maximin" -- maximin

              Default: "MD2"}
}
\value{
 criterion value.
}

\author{
   Aijun Zhang, Haoyu Li, Shijie Quan
}
\references{
 Zhang, A. and Li, H. (2017). UniDOE: An R package for constructing uniform design of experiments via stochastic and adaptive threshold accepting algorithm. Technical Report.
}

\examples{

  x0 = matrix(c(1,1,1,2,2,2,3,3,3),nrow=3,byrow=TRUE)
  crit="MD2"
  value = DesignEval(x0,crit)

}
