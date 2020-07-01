\name{smquad.face3d}

\alias{smquad.face3d}

\title{Locating ridges and valleys by quadratic fitting.}

\description{
  This function creates a smooth ridge path through a quadratic assumption for the ridge.
}

\usage{
smquad.face3d(x, pdist, y, lambda, df = 5, method = "df",
                        weights, eval.points, ngrid,
                        nseg, bdeg = 3, pord = 2, display = "lines", fixed) 
}

\arguments{
\item{x}{abcissae of data}
\item{pdist}{response}
\item{y}{values of pdist}
\item{lambda}{smoothness parameter}
\item{df}{degrees of freedom}
\item{method}{method of smoothing}
\item{weights}{weights attached to observations}
\item{eval.points}{points at which the estimate is evaluated.}
\item{ngrid}{the number of regular grid points in \code{eval.points}.}
\item{nseg}{number of intervals for B-splines}
\item{bdeg}{degreee of B-splines}
\item{pord}{order of difference penalty}
\item{display}{a character variable which identified the form of display.}
\item{fixed}{a matrix containing fixed points.}
}

\details{
This function creates a smooth ridge path through a quadratic assumption for the ridge.
}

\value{
  None
}

