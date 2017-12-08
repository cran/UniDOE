\name{GenLP}
\alias{GenLP}
\title{Generating uniform designs by level permutation}
\usage{
  GenLP(X0, crit, maxiter, hits_ratio, vis)
}
\description{
Function for generating uniform designs by level permutation. It can be also used to improve a user-specified design with randomized level permutation.
}
\arguments{
\item{X0}{Current design }

\item{crit}{an character R object. Type of criterion to use.

  "maximin" -- maximin Discrepancy ;

  "CL2"  --Centered L2 Discrepancy ;

  "MD2"  --Mixture L2  Discrepancy ;}

\item{maxiter}{a positive integer  R  object}


\item{hits_ratio}{an float R object. Default value is 0.1, which is the ratio to accept changes of design in inner for loop. Details can be checked in (Zhang, A. and Li, H. (2017). UniDOE: An R package for constructing uniform design of experiments via stochastic and adaptive threshold accepting algorithm. Technical Report.)}

\item{vis}{an boolean R object. If true, plot the criterion value sequence.}

}
\value{
  A list that contains Initial design matrix(initial_design),optimal design matrix(final_design), initial criterion value(initial_criterion), final criterion value(criterion_value) and criterion list(criterion_lists) in update process.
}

\author{
   Aijun Zhang, Haoyu Li, Shijie Quan
}
\references{
 Zhang, A. and Li, H. (2017). UniDOE: An R package for constructing uniform design of experiments via stochastic and adaptive threshold accepting algorithm. Technical Report.
}

\examples{

  n=27#(must be multiples of q)
  s=13
  q=3
  crit = "MD2" #(Mixture L2 criteria)
  res1 = GenUD(n,s,q,crit=crit,maxiter=100)
  res2 = GenLP(res1$final_design,vis=TRUE, maxiter=100) # To improve ...
}
