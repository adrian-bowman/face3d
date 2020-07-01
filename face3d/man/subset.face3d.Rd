\name{subset.face3d}

\alias{subset.face3d}

\title{Create a subset of a face3d object}

\description{
  This function creates a subset of a face3d object.
}

\usage{
\method{subset}{face3d}(shape, subset, remove.singles = TRUE, retain.indices = FALSE, \dots)
}

\arguments{
\item{shape}{a \code{face3d} object.}
\item{subset}{a set of indices, or a logical vector with length equal to the number of rows of \code{shape$coords}, which defines the subset.}
\item{remove.singles}{a logical variable which determines whether individual points, not connected to the triangular mesh, are removed.  When \code{remove.singles} is \code{TRUE} and isolated points are present, the subset object will have fewer points thanthose specified in the input argument \code{ind}.  The indices defining the final subset of the original object \code{shape} can be identified as \code{as.numeric(rownames(object))}, where \code{object} denotes the subset object.}
\item{retain.indices}{a logical variable which determines whether the indices defining the subset of the input \code{shape} are placed in a component \code{subset} of the returned face3d object.  See Details below for further information.}
\item{...}{additional arguments.  These are not used.  They are present simply to provide compatability with the generic \code{subset} function.}
}

\details{
The function returns a new face3d object which has \code{coords} and \code{colour} reduced to match the required subset and with \code{triples} constructed to index the triangular faces of the new object.  The rownames of \code{coords} are transferred to the new object.  Isolated points which are not part of any triangles are removed.

When \code{remove.singles = TRUE}, the subset returned may be different from the subset defined by the argument \code{subset}.  The indices of the returned subset can be returned by setting \code{retain.indices = TRUE}.  The indices are then placed in the component \code{subset} of the returned face3d object.

In addition to the components \code{coords}, \code{triples} and \code{colour}, any vector, matrix, array or list components of \code{shape} which have exactly one dimension matching \code{nrow(shape$coords)} are transferred to the returned face3d object after applying the subset operation defined by \code{ind}.  Other information is transferred unchanged.
}

\value{
A face3d object in the form of a list with two or three components.  The component \code{coords} is a matrix with three columns containing the three-dimensional co-ordinates of the points.  The component \code{triples} is a vector of integers containg the indices of \code{coords} which should be passed to the \pkg{rgl} function \code{triangles} to render the surface.  If the third component \code{colour} is present in the original object the the corresponding components of the new object contains the colours of the points corresponding to the rows of \code{coords}.

When \code{retain.indices = TRUE}, the returned object contains a component \code{subset} which gives the row indices of \code{shape$coords} which correspond to the subset created by \code{ind}.  This may not be identical to the definition of the subset defined by \code{ind} when \code{remove.singles = TRUE}.  If \code{shape} contains a \code{subset} component this will not be transferred to the returned object.

As described in Details above, information present in \code{shape} which has exactly one dimension matching \code{nrow(shape$coords)} is also subsetted.  Other information in \code{shape} is transferred without change.
}

\examples{
\dontrun{
sbst <- subset.face3d(template_male, template_male$shape.index > 0)
plot(sbst, col = "shape index")
plot(sbst, display = "normals", add = TRUE)

plot(template_male)
fn   <- rgl::select3d()
sbst <- subset.face3d(template_male, fn(template_male$coords))
plot(sbst)
}}
