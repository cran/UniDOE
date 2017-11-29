\name{UDC}
\alias{UDC}
\title{Uniform Design Construction}
\usage{
UDC(n,s,q,init,initX,crit,maxiter,hits_ratio,vis)
}
\description{
  This function takes n,s,q and other arguments to output a list(described below).
}
\arguments{
\item{n}{an integer R object. Run of Experiment}

\item{s}{an integer R object. Factor of Experiment.}

\item{q}{an integer R object. Level of Experiment.}

\item{init}{an string vector R object: "rand"(default), "orth" or "input". Criterion used}

\item{initX}{an user-defined integer matrix R object, which is fixed in search. Needed when init="orth" or init="input"}

\item{crit}{an character R object. Type of criterion to use.

         "maximin" -- maximin Discrepancy ;

         "CL2"  --Centered L2 Discrepancy ;

         "MD2"  --Mixture L2  Discrepancy ;}

\item{maxiter}{a positive integer  R  object. Maximum iteration number in outer while loop of SOAT algorithm.}


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
##e.g.1
n=12#(must be multiples of q)
s=3
q=4
crit = "MD2"#(Mixture L2 criteria)
list1=UDC(n,s,q,crit=crit,maxiter=100)

##e.g.2
n=10 #(must be multiples of q)
s=3
q=5
init = "rand"
crit = "MD2" #(Mixture L2 criteria)
vis=TRUE
list1=UDC(n,s,q,init=init,crit=crit,maxiter=100,vis=vis)

##e.g.3
#If init = "orth",
#algorithm will fix columns of initX and expand total number of columns to s for final design.
n=6
s=5
q=3
init = "orth"
# initX should be orthogonal design with nrow=n&level=q
initX = matrix(c(1,1,2,2,3,3,3,3,1,1,2,2),ncol=2)
list = UDC(n,s,q,init=init,initX = initX,maxiter=100)

##e.g.4
#If init="input", algorithm will search for better a better design with same size as initX.
init = "input"
# initX should be orthogonal design with nrow=n&level=q
initX = matrix(c(1,1,2,2,3,3,3,3,1,1,2,2),ncol=2)
list = UDC(init=init,initX = initX,maxiter=100)

}
