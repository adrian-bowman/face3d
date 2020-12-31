\name{distance.face3d}

\alias{distance.face3d}

\title{Comparing two shapes which have point-to-point correspondece}

\description{
The calculation of several types of distance between two shapes, where the rows of \code{shape1$coords} and \code{shape2$coords} have geometrical correspondence.
}

\usage{
distance.face3d(shape1, shape2)}

\arguments{
\item{shape1, shape2}{two \code{face3d} objects.}
}

\details{
Several different distances between \code{shape1} and \code{shape2} are computed.  These are documented below in the Values section.
}

\value{
A list with components which contain different distances.  In each case, the distances are signed, with \code{shape1} as the reference.  Linear distances are from \code{shape1} to \code{shape2}.  The areal distance gives the ratio of the areas of the triangles in \code{shape2} to the corresponding triangular areas in \code{shape1}, on a log scale.
\item{x,y,z}{the distances between the corresponding points of the two shapes in the x, y and z co-ordinate axes.}
\item{xyz}{the Euclidean distances between the corresponding points of the two shapes.}
\item{normal}{the distances between the corresponding points of the two shapes, projected along the normal directions of \code{shape1}.}
\item{log.area.ratio}{the log ratio of the areas of the corresponding triangles of the two shapes.}
}
