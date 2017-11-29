\name{LP}
\alias{LP}
\title{Improve design with level permutation}
\usage{
  LP( X0, crit, maxiter, hits_ratio, vis)
}
\description{
  This function aims at improving current design using exhaustive level permutation.
}
\arguments{
\item{X0}{Current design }

\item{crit}{an character R object. Type of criterion to use.

  "maximin" -- maximin Discrepancy ;

  "CL2"  --Centered L2 Discrepancy ;

  "MD2"  --Mixture L2  Discrepancy ;}

\item{maxiter}{a positive integer  R  object}


\item{hits_ratio}{an float R object. Default value is 0.1, which is the ratio to accept changes of design in inner for loop. Details can be checked in (Zhang, A. and Li, H. (2017). UniDOE: an R package for uniform design construction via stochastic optimization.)}

\item{vis}{an boolean R object. If true, plot the criterion value sequence.}

}
\value{
  A list that contains Initial design matrix(initial_design),optimal design matrix(final_design), initial criterion value(initial_criterion), final criterion value(criterion_value) and criterion list(criterion_lists) in update process.
}

\author{
   Aijun Zhang, Haoyu Li
}
\references{
 Zhang, A. and Li, H. (2017). UniDOE: An R package for constructing uniform design of experiments via stochastic and adaptive threshold accepting algorithm. Technical Report.
}

\examples{

  n=27#(must be multiples of q)
  s=13
  q=3
  crit = "MD2"#(Mixture L2 criteria)
  list1=UDC(n,s,q,crit=crit,maxiter=100)
  LP_list = LP(list1$final_design,vis=TRUE,maxiter=100)
}
