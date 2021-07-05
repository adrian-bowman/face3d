\name{distance.face3d}

\alias{distance.face3d}

\title{The distances between points or shape vertices}

\description{
Euclidean distances between 3d co-ordinates, or more complex distances between \code{face3d} objects, are computed.  The Details section below provides a description of the distances computed between different types of input object.
}

\usage{
distance.face3d(x1, x2)}

\arguments{
\item{x1, x2}{each of \code{x1} and {x2} may be either a vector of 3d co-ordinates, a three-column matrix whose rows are 3d co-ordinates, or a \code{face3d} object.}
}

\details{
When either \code{x1} or \code{x2} is a vector or matrix of 3d co-ordinates, the Euclidean distances are calculated.  When one of \code{x1} or \code{x2} is a \code{face3d} object then the distances are calculated to the \code{vertices} of this object. In either case, these distances are computed by calling the \code{rdist} function from \pkg{fields} package.

When \code{x1} and {x2} are both \code{face3d} objects whose \code{vertices} have geometrical correspondence, several types of distance are calculated from the \code{vertices} of \code{x1} to the \code{vertices} of \code{x2}. These are documented below in the Values section.
}

\value{
When at least one of \code{x1} and \code{x2} is a matrix of 3d co-ordinates, a matrix of distances between the relevant points is returned. If \code{x1} or \code{x2} is a vector corresponding to a single 3d point then a vector of distances is returned.

When \code{x1} and \code{x2} are both \code{face3d} objects whose \code{vertices} have geometrical correspondence, a list with components which contain different distances.  In each case, the distances are signed, with \code{x1} as the reference.  Linear distances are from \code{x1} to \code{x2}.  The areal distance gives the ratio of the areas of the triangles in \code{x2} to the corresponding triangular areas in \code{x1}, on a log scale.
\item{x,y,z}{the distances between the corresponding points of the two shapes in the x, y and z co-ordinate axes.}
\item{xyz}{the Euclidean distances between the corresponding points of the two shapes.}
\item{normal}{the distances between the corresponding points of the two shapes, projected along the normal directions of \code{x1}.}
\item{log.area.ratio}{the log ratio of the areas of the corresponding triangles of the two shapes.}
}

\seealso{
The function \code{closest.face3d} computes the distances of points to a \code{face3d} object by also considering distances to edges and triangular faces.

The function \code{project.face3d} projects points onto a \code{face3d} object, by moving along the directions specified in the corresponding \code{normals} vectors
}
