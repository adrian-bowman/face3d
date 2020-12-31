\name{sm.map}

\alias{sm.map}

\title{Finding the appropriate smoothing parameter in 2 dimensions for a curve through a 2d points space}

\description{
This function finds the appropriate level of smoothing through a 2d point space.
}

\usage{
sm.map(x, y, weights, fixed, nseg = 25, deg = 3, pdiff = 2,
             ig.a = 1e-3, ig.b = 1e-3, lambdas = exp(seq(-10, 10, length = 50)),
             prior = function(lambda) 1,
             ngrid = 401, eval.points = seq(min(x), max(x), length = ngrid))
}

\arguments{
\item{x}{covariate data.}
\item{y}{response data.}
\item{weights}{the weights attached to each observation.}
\item{fixed}{a matrix of points which should be held fixed.}
\item{nseg}{the number of segments used in the p-splines.}
\item{deg}{the degree of the b-splines.}
\item{pdiff}{the order of differencing used for the penalty.}
\item{ig.a, ig.b}{the parameters of the inverse Gaussian prior for the variance.}
\item{lambdas}{a vector of values for the penalty parameter \code{lambda} over which the solution will be sought.}
\item{prior}{a function which defines a prior over the penalty parameter \code{lambda}.}
\item{ngrid}{the number of grid points at which the estimate will be evaluated.}
\item{eval.points}{a vector of points at which the estimate will be evaluated.}
}

\details{
This function finds the appropriate level of smoothing through a 2d point space.
}

\value{
  None
}

% \examples{
% \dontrun{
  
% }}
