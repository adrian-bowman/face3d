\name{princurve.face3d}

\alias{princurve.face3d}

\title{Find a principal curve through a set of points}

\description{
  This function implements the idea of a principal curve, which is a smooth path through a set of observed points, in two or three dimensions.
}

\usage{
princurve.face3d(x, df = 12, maxiter = 20,
                          monitor = FALSE, threshold = 1e-3)
}

\arguments{
\item{x}{a matrix with two or three columns containing the co-ordinates of the points in two or three dimensions to be used to contruct the principal curve.}
\item{df}{the degrees of freedom which specify the smoothness required in the relationship between each co-ordinate and the arc length of the curve.}
\item{maxiter}{The maximum permitted number of iterations of the principal curve algorithm.}  
\item{threshold}{The value of the relative absolute difference in sum-of-squares below which the iterations will be halted.}  
}

\details{
Reference to Hastie & Tibshirani (?) 1989, JASA.
}

\value{
A list with the following components:
\item{s}{a matrix of principal curve fitted values corresponding to the observed values \code{x}.}
\item{lambda}{the arc length values which correspond to the rows of \code{s}.}
\item{ord}{the permutation required to put the entries of \code{s} and \code{lambda} into sequential order along the curve.}
\item{dist}{the minimised sum-of-squares of the observations from the curve.}
}

\references{
Hastie ...
}

\examples{
\dontrun{
}}
