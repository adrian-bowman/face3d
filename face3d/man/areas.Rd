\name{areas.face3d}

\alias{areas.face3d}

\title{Calculate the areas of the triangles in an object}

\description{
This function calculates the area of each triangle in the triangulation of a \code{face3d} object.  The total surface area of the object is also returned.  Areas are also associated with each vertex.
}

\usage{
areas(shape)
}

\arguments{
\item{shape}{a \code{face3d} object.}
}

\details{
The function calculates the area of each triangle which makes up the \code{face3d} object.  Care needs to be taken when the function is applied to a raw image, as the surface area may be affected by complex, and potentially noisy, surface shapes at the perimeter.

There is some provision to associate area with curves in a model mesh but this is not yet consistent.  In this case, the \code{coords} of shape needs to have rownames.
}

\value{
A list with three components:
\item{area}{the total area of \code{shape}.}
\item{triangles}{a vector containing the area of each triangle.  This vector is in the order which corresponds to the triangles defined by \code{shape$triples}.}
\item{points}{a vector containing the area around each vertex.  This vector is in the order which corresponds to the points defined by \code{shape$coords}.  The areas are computed as one third of the sum of the areas of the triangles which have a particular point as a vertex.}
}

\examples{
\dontrun{
nose <- subset(template_male, edist(template_male$landmarks["pn", ], template_male$vertices) < 30)
ars  <- areas(nose)
}}
