\name{closestcurve.face3d}

\alias{closestcurve.face3d}

\title{Find the closest distance from a curve in 3 dimensions to all points on the surface containing it.}

\description{
  This function constructs finds the closest distance from a curve in 3 dimensions to all points on the surface containing it.
}
\usage{
closestcurve.face3d(shape, curve) 
}
\arguments{
\item{shape}{A face3d object or a subset thereof.}
\item{curve}{A n by 3 matrix where the rows are the individual curve points and the columns are the dimensions.}
}

\details{
The closest distance from every point on the shape surface to the curve is given, along with the indice of the closest curve point.  }

\value{
  A list with a component \code{closest.distance} providing the distance between every surface point and its closest curve point.  In addition a component \code{closest.curvept} providing the indice of the closest curve point to each surface point.
}

\examples{
\dontrun{
  ccdist  <- closestcurve.face3d(shape, curve) 
  lambda  <- cumsum(c(0, sqrt(apply((diff(curve))^2, 1, sum))))
  xcoord  <- lambda[ccdist$closest.curvept] # distance along the curve
  ycoord  <- ccdist$closest.distance # perpindicular distance from curve 
 }}
