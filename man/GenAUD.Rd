\name{GenAUD}
\alias{GenAUD}
\title{Generating Augmented Uniform Design of Experiments}
\usage{
GenAUD(X0,n,crit,maxiter,hits_ratio,vis)
}
\description{
    This function takes n,s,q; a unchanged initial design and other arguments to output a list (described below).}
\arguments{
\item{X0}{an integer matrix R object}

\item{n}{an integer R object. Run of Experiment}

\item{crit}{an character R object. Type of criterion to use.

      "maximin" -- maximin Discrepancy;

      "CD2"     --Centered L2 Discrepancy;

      "WD2" -- Wrap-around L2 Discrepancy;

      "MD2"     --Mixture L2  Discrepancy;}

\item{maxiter}{a positive integer  R  object. Maximum iteration number in outer while loop of SOAT algorithm.}

\item{hits_ratio}{an float R object. Default value is 0.1, which is the ratio to accept changes of design in inner for loop. Details can be checked in (Zhang, A. and Li, H. (2017). UniDOE: An R package for constructing uniform design of experiments via stochastic and adaptive threshold accepting algorithm. Technical Report.)}

\item{vis}{an boolean variable. If true, plot the trace of criterion values.}
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
#Example 1.
#Set a fixed initial matrix:
n=12 #(must be multiples of q)
mat0 = matrix(c(1,1,1,2,2,2,3,3,3),ncol=3,byrow=TRUE)# nb. of columns=s
crit = "MD2" #(Mixture L2 criteria)
res = GenAUD(X0=mat0,n,crit=crit,maxiter=100)

# Example 2.
# Set a fixed initial matrix with visualization:
n=9 #(must be multiples of q)
mat0 = matrix(c(1,1,1,2,2,2,3,3,3), ncol = 3, byrow = TRUE)
crit = "MD2" #(Mixture L2 criteria)
list1=GenAUD(X0=mat0,n,crit=crit,vis=TRUE,maxiter=100)
}
