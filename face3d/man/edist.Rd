\name{edist}

\alias{edist}

\title{Euclidean distances between 3D co-ordinates}

\description{
Euclidean distances between 3d co-ordinates are computed.  The Details section below provides a description of the distances computed between different types of input object.
}

\usage{
edist(x1, x2, type = "pairwise", index = FALSE)
}

\arguments{
\item{x1, x2}{each of \code{x1} and {x2} may be either a vector of 3d co-ordinates or a three-column matrix whose rows are 3d co-ordinates.  One of \code{x1} or \code{x2} may be a \code{face3d} object.  \code{x2} may be omitted.}
\item{type}{a character string which determines which type of distance is computed.}
\item{index}{a logiocal argument which, in the case of \code{type = "min"}, determines whether an index is retured to identify which point in \code{x2} is closest to each point of \coide{x1}.}
}

\details{
When either \code{x1} or \code{x2} is a vector or matrix of 3d co-ordinates, Euclidean distances are calculated.  When one of \code{x1} or \code{x2} is a \code{face3d} object then the distances are calculated from or to the \code{vertices} of this object. In either case, these distances are computed by calling the \code{rdist} function from \pkg{fields} package.
}

\value{
When at least one of \code{x1} and \code{x2} is a matrix of 3d co-ordinates, a matrix of distances between the relevant points is returned. If \code{x1} or \code{x2} is a vector corresponding to a single 3d point then a vector of distances is returned.

When \code{type = "min"}, the minimum distance from each set of co-ordinates in \code{x1} to those in \code{x2} is computed.  If \code{index = FALSE} then the vector of minimum distanvces is returned.  If \code{index = TRUE} then a list with two components is returned.  The minimum distances are in the component \code{distance},  Indices which identify the point in \code{x2} which is closest to each point in \code{x1} are returned in the component \code{indices}.

When \code{type = "minum"}, the minimum distance from each set of co-ordinates in \code{x1} to those in \code{x2} is computed and these are summed to provide a single overall distance measure.

When \code{x2} is missing, the lengths of the vector, matrix or \code{face3d} object vertices are returned.
}

\seealso{
The function \code{closest.face3d} computes the distances of points to a \code{face3d} object by also considering distances to edges and triangular faces.

The function \code{project.face3d} projects points onto a \code{face3d} object, by moving along the directions specified in the corresponding \code{normals} vectors.

The function \code{compare} calculates a variety of distances between two \code{face3d} objects with the same number of vertices.
}
