\name{DesignQuery}
\alias{DesignQuery}
\title{Evaluate design in terms of criteria}
\usage{
DesignQuery(n,s,q,crit, ShowCrit)
}
\description{
  This function takes size of desired design,criterion crit.
  If the required design exists in database, then return the design, else return NULL.
}
\arguments{
\item{n}{an integer R object. Run of Experiment}

\item{s}{an integer R object. Factor of Experiment.}

\item{q}{an integer R object. Level of Experiment.}

\item{crit}{an integer object, criterion to choose:

              "MD2" -- MD2

              "CD2" -- CD2

              "maximin" -- maximin

              Default: "CD2"}
\item{ShowCrit}{If TRUE, print CD2,MD2,maximin value of required design.

                Default: TRUE}

}

\value{
 Desired design
}

\author{
   Aijun Zhang, Haoyu Li, Shijie Quan
}
\references{
 Zhang, A. and Li, H. (2017). UniDOE: An R package for constructing uniform design of experiments via stochastic and adaptive threshold accepting algorithm. Technical Report.
}

\examples{
  n = 9
  s = 3
  q = 3
  crit="MD2"
  D = DesignQuery(n,s,q,crit)
  D
}
