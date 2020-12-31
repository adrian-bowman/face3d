\name{gcurvature.face3d}

\alias{gcurvature.face3d}

\title{The geodesic curvature of a three-dimensional curve.}

\description{
The curvature of a geodesic curve, derived from a three-dimensional surface, is calculated.
}

\usage{
gcurvature.face3d(curve, df, n = 200, graphics = FALSE)
}

\arguments{
\item{curve}{a matrix with three columns containing the co-ordinate positions of points along the curve.}
\item{df}{an integer defining the degrees of freedom used when smoothing the curve points to estimate curvature.}
\item{n}{an integer determining the number of points to be created in the resampled version of the curve.}
\item{n}{a logical values which determines whether the curvature is plotted in the standard graphics window..}
}

\details{
The points in \code{curve} are assumed to lie on the lines defining the triangulation of the surface.  The resampled curve is then created by the \code{approx} function, using linear interpolation so that the resampled points also lie on the surface.  A smoothing technique is used to estimate the derivatives required to calculate curvature.
}

\value{
The function returns a list with the following components.
\item{resampled.curve}{a matrix with three columns giving the co-ordinate positions of the points on the resampled curve.}
\item{gcurvature}{a vector of estimated curvatures at each point on the resampled curve.}
\item{arclength}{a vector of arc length locations of each point on the resampled curve.}
\item{pos.localmax}{the 3-d positions of points which sit at local maxima of the curvature.}
\item{ind.localmax}{a vector of the indices of the local maxima of the curvature.}
\item{pos.max}{the 3-d position of the point on the curve where the curvature is largest.}
\item{ind.max}{the index of the 3-d position of the point on the curve where the curvature is largest.}
\item{d2.x, d2.y, d2.z}{the second derivatives with respect to x, y and z along the curve.}
}
