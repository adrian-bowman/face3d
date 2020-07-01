\name{interp.barycentric}

\alias{interp.barycentric}

\title{A barycentric interpolation of 2d points to 3d points.}

\description{
  This function interpolates a set of two dimensionsal points into three dimensional points. It is used specifically for projecting a 2 dimensional set of points at 3d curve.
}

\usage{
interp.barycentric(X, f, Xnew)
}

\arguments{
\item{X}{A n by 2 matrix where column 1 is the x coordinates and column 2 is the y coordinates (the full space).}
\item{f}{The value of the three dimensional coordinates the full set of two dimensional coordinates were projected from.}
\item{Xnew}{A new set of points (the curve) in the same two dimensional space.}
}

\details{
The 2 dimensional curve points are returned as three dimensional curve points with associated values attached
}

\value{
 A list containing the component \code{Xnew} and \code{fnew} where \code{fnew} is the coordinates of the curve in three dimensions and \code{Xnew} is the associated points in 2d.
}

\examples{
\dontrun{
   interp.x      <- interp.barycentric(X, f = c(lmk1[1], shape$coords[ , 1], lmk2[1]), Xnew)$fnew
   interp.y      <- interp.barycentric(X, f = c(lmk1[2], shape$coords[ , 2], lmk2[2]), Xnew)$fnew
   interp.z      <- interp.barycentric(X, f = c(lmk1[3], shape$coords[ , 3], lmk2[3]), Xnew)$fnew
   new.curve   <- cbind(interp.x, interp.y, interp.z)

}}
