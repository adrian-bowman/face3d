\name{closest.face3d}

\alias{closest.face3d}

\title{Find the closest point in a face3d object}

\description{
  This function finds the closest point to \code{x} in the face3d object \code{shape}.
}

\usage{
closest.face3d(x, shape, nearest = 30)
}

\arguments{
\item{x}{a vector of length 3 giving a co-ordinate position in three-dimensional space.}
\item{shape}{a \code{face3d} object.}
\item{nearest}{a positive integer controlling the number of nearest neighbours which will be used to define the neighbourhood within which the closest point is sought.}
}

\details{
A subset of \code{shape} is defined, using the \code{subset.face3d} function, through the  \code{nearest} points within \code{shape$coords} which are closest to \code{x}.  This neighbourhood is then searched to identify the closest point to \code{x}.  This closest point may be a vertex of \code{shape}, it may lie on an edge or it may be an interior point of a tringle within the mesh.
}

\value{
A list with two components: 
\item{point}{a vector of length 3 which defines the closest point.}
\item{id}{an integer vector which identifies the points, as rows of \code{shape$coords}, which encompass the closest points; where \code{point} is a vertex this will be of length 1, where it lies on an edge this will be of length 2 to define the end points, and where it is an interior point of a triangle this will be of length 3 giving the indices of the triangle vertices.}
\item{distance}{the distance of the closest point from \code{x}.}
\item{subset}{a \code{face3d} object containing the subset of \code{shape} which was used to search for the closest point.}
}
