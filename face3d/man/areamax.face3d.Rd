\name{areamax}

\alias{areamax}

\title{Locate the maximum local areal integral of a variable}

\description{
The integral of a variable over a radial area is computed at all vertices within a \code{face3d} object and the location which maximises this is identified.
}

\usage{
areamax(shape, values, radius)
}

\arguments{
\item{shape}{a \code{face3d} object.}
\item{values}{a vector of values at the locations in \code{shape$vertices}.}
\item{radius}{a numerical value specifying the radius around each vertex over which the integral of \code{values} is computed.}
}

\details{
The use of integration over a local area helps to ameliorate the problem of individual points of local spikes with very high values dominating.

At each location within \code{shape$vertices} the areal integral is approximated by the sum of \code{values} weighted by the areas associated with each element of \code{shape$vertices} within a distance of \code{radius}; see \code{\link{areas}} for the calculation of the areas.
}

\value{
A list with three components:
\item{point}{a vector of length three giving the location of the maximum areal integral.}
\item{value}{the value of the corresponding areal integral.}
\item{id}{an integer giving the row index of \code{point} in the \code{shape$vertices} matrix.}
}

\examples{
\dontrun{
nose <- subset(template_male, edist(template_male$landmarks["pn", ], template_male$vertices) < 30)
nose <- curvature(nose, distance = 10)
gc   <- nose$kappa1 * nose$kappa2
pt   <- areamax(nose, gc, 10)$point
plot(nose, col = gc)
spheres3d(pt, col = "green")
}}
