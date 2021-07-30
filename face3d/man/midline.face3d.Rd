\name{midline.face3d}

\alias{midline.face3d}

\title{Locating the line of symmetry in a shape}

\description{
This function identified the midline of a shape by minimising local asymmetry.
}

\usage{
midline.face3d(shape, x1, x2, nv = 51, nh = 51, dh = 50, lambda = 1, monitor = FALSE)
}

\arguments{
\item{shape}{a \code{face3d} object.}
\item{x1, x2}{vectors of length 3 which define two landmarks through which the midline is expected to run.  These landmarks can be approximate.}
\item{nv, nh}{integers defining the number of grid points for the search area, in vertical and horizontal directions respectively.}
\item{dh}{a numeric value defining the spatial extent of the horizontal transects to be searched on each side of the initial midline.}
\item{lambda}{a nunmerica value which defines the weight attached to the roughness penalty used when constructing the midline curve.}
\item{monitor}{a logical value which controls whether intermediate plots are created, to monitor the estimation process.}
}

\details{
An initial estimate of the midline is created by applying the \code{planepath} function, using \code{x1} and \code{x2} as end-points.  At each vertical position on this initial curve, narrow horizontal strips centred at different horizontal locations are used to measure asymmetry and the location with minimum asymmetry is found.  The final curve is produced with the use of a roughness penalty.
}

\value{
A matrix with three columns whose rows give a regularly spaced set of points along the estimated midline curve.
}

\examples{
\dontrun{
x1   <- c(18, -7, 85)
x1   <- closest.face3d(x1, face)$point
x2   <- c(17.47, -75, 96.70)
x2   <- closest.face3d(x2, face)$point
path <- midline.face3d(face, x1, x2, lambda = 1, monitor = TRUE)
plot(face, new = FALSE)
spheres3d(rbind(x1, x2), col = "red")
}}
