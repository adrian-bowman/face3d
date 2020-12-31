\name{subset.face3d}

\alias{subset.face3d}

\title{Create a subset of a face3d object}

\description{
  This function creates a subset of a face3d object.
}

\usage{
\method{subset}{face3d}(x, ind, remove.singles = TRUE, \dots)
}

\arguments{
\item{x}{a \code{face3d} object.}
\item{ind}{a set of indices, or a logical vector with length equal to the number of rows of \code{shape$coords}, which defines the subset.}
\item{remove.singles}{a logical variable which determines whether individual points, not connected to the triangular mesh, are removed.  When \code{remove.singles} is \code{TRUE} and isolated points are present, the subset object will have fewer points thanthose specified in the input argument \code{ind}.  The indices defining the final subset of the original object \code{x} can be identified as \code{as.numeric(rownames(object))}, where \code{object} denotes the subset object.}
\item{...}{additional arguments.  These are not used.  They are present simply to provide compatability with the generic \code{subset} function.}
}

\details{
The function returns a new face3d object which has \code{coords} and \code{colour} reduced to match the required subset and with \code{triples} constructed to index the triangular faces of the new object.  The rownames of \code{coords} give the indices of the points from \code{x}.  Isolated points which are not part of any triangles are removed. 
}

\value{
A face3d object in the form of a list with two or three components.  The component \code{coords} is a matrix with three columns containing the three-dimensional co-ordinates of the points.  The component \code{triples} is a vector of integers containg the indices of \code{coords} which should be passed to the \pkg{rgl} function \code{triangles} to render the surface.  If the third component \code{colour} is present in the original object the the corresponding components of the new object contains the colours of the points corresponding to the rows of \code{coords}.
}

\examples{
\dontrun{
data(face)
display.face3d(face)

face1 <- subset.face3d(face, face$coords[ , 1] > 0)
clear3d()
display.face3d(face1, add = TRUE)

clear3d()
display.face3d(face1, subset = face$coords[ , 1] > 0, add = TRUE)

display.face3d(face)
fn    <- select3d()
ind   <- which(fn(face$coords))
face2 <- subset.face3d(face, ind)
clear3d()
display.face3d(face2, add = TRUE)
}}
