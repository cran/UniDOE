\name{GenAUD_MS}
\alias{GenAUD_MS}
\title{Generating sequential Uniform Design of Experiments using diffrent initial designs}
\usage{
GenAUD_MS(X0, n, crit, maxiter, nshoot, vis=FALSE)
}
\description{
  This function takes n,s,q and other arguments to output a list(described below).
}
\arguments{
\item{X0}{an integer matrix R object. Fixed design to be used.}

\item{n}{an integer R object. Number of rows to be added to design}

\item{crit}{an character R object. Type of criterion to use.

         "maximin" -- maximin Discrepancy ;

         "CD2"  -- Centered L2 Discrepancy ;

         "WD2" -- Wrap-around L2 Discrepancy;

         "MD2"  -- Mixture L2 Discrepancy ;}

\item{maxiter}{a positive integer  R  object. Maximum iteration number in outer while loop of SATA algorithm in each shoot.}

\item{nshoot}{Total counts to try different initial designs.}

\item{vis}{an boolean R object. If true, plot the criterion value sequence for all shoots.}
}

\value{
  Best design over all shoots.
}

\author{
   Aijun Zhang, Haoyu Li, Shijie Quan
}

\references{
 Zhang, A. and Li, H. (2017). UniDOE: An R package for constructing uniform design of experiments via stochastic and adaptive threshold accepting algorithm. Technical Report.
}


\examples{
D1 = DesignQuery(n = 18, s = 7, q = 3, crit = "MD2")
D2 = GenAUD_MS(X0=D1+1, n=18, crit="MD2",
               maxiter = 100, nshoot = 5,
               vis = TRUE)
}
