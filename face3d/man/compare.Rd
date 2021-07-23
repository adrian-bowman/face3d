\name{compare.face3d}

\alias{compare.face3d}

\title{Euclidean and other distances between the vertices of two \code{face3d} objects}

\description{
Euclidean distances between the vertices of two \code{face3d} obje cts are computed.  The Details section below provides a description of the distances produced.
}

\usage{
compare(x1, x2)
}

\arguments{
\item{x1, x2}{each of \code{x1} and {x2} should be a \code{face3d} object.}
}

\details{
\code{x1} and {x2} must be \code{face3d} objects whose \code{vertices} have geometrical correspondence. Several types of distance are calculated from the \code{vertices} of \code{x1} to the \code{vertices} of \code{x2}. These are documented below in the Values section.  A particular use of these distances is in quantifying and plotting the differences between the model or template representations of two objects.
}

\value{
The function returns a list with components which contain different distances.  In each case, the distances take \code{x1} as the reference object.  Linear distances are signed and correspond to the distances from \code{x1} to \code{x2}.  The areal distance gives the ratio of the areas of the triangles in \code{x2} to the corresponding triangular areas in \code{x1}, on a log scale.
\item{x,y,z}{the distances between the corresponding points of the two shapes in the x, y and z co-ordinate axes.}
\item{xyz}{the Euclidean distances between the corresponding points of the two shapes.}
\item{normal}{the distances between the corresponding points of the two shapes, projected along the normal directions of \code{x1}.}
\item{log.area.ratio}{the log ratio of the areas of the corresponding triangles of the two shapes.}
}

\seealso{
The function \code{closest.face3d} computes the distances of points to a \code{face3d} object by also considering distances to edges and triangular faces.

The function \code{project.face3d} projects points onto a \code{face3d} object, by moving along the directions specified in the corresponding \code{normals} vectors
}
