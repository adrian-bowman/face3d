\name{projectline.face3d}

\alias{projectline.face3d}

\title{Project the co-ordinates of a shape onto the line between two points.}

\description{
This function projects the co-ordinates of a \code{shape} object onto the line between the points \code{x1} and \code{x2}.  If a \code{normal} vectr is supplied then the distance have an appropriate sign attached.
}

\usage{
projectline.face3d(shape, x1, x2, normal)
}

\arguments{
\item{shape}{a \code{face3d} object.}
\item{x1, x2}{vectors of length three defining points which lie on the projection line.}
\item{normal}{an optional vector of length three which defines a normal.  This will be used to attached a sign to the distances associated with each of the co-ordinates in the \code{face3d} object.}
}

\details{
This function finds the perpendicular projection of the points in \code{shape$coords} onto the line containing \code{x1} and \code{x2}.  It returns the distance of each projected point along this line and the distance of the projected point from the corresponding co-ordinate of \code{shape$coords}.  These distances are returned in the \code{x} and \code{y} components of the result of the function.  The origin is at \code{x1} so that points which project there have \code{x = 0}.

If a \code{normal} vector is supplied then a new axis is created from the crossproduct of \code{x2 - x1} and \code{normal}.  The returned disstances \code{y} then have a sign attached, determined by the sign of the projection of each element of \code{shape$coords}, re-centred at \code{x1}, onto this new axis.
}

\value{
A list with two components:
\item{x}{the distance along the projection line.  The origin is at \code{x1} so that points which project there have \code{x = 0}.}
\item{y}{the perpendicular distance from each projection point to each element of \code{shape$coords}.}
}

\examples{
\dontrun{
   # To be written.
}}
