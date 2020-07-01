\name{as.face3d}

\alias{as.face3d}

\title{Create a face3d object}

\description{
This function reates a \code{face3d} object from appropriate input.
}

\usage{
as.face3d(object, model.mesh = FALSE)
}

\arguments{
\item{object}{this may be a list whose components include \code{coords}, containing a matrix  of 3-d point positions, and \code{triples}, a vector of indices of this points where each successive triple identifies a triangle on the surface of the object.  Alternatively, a matrix of points corresponding to a facial model mesh can be passed when \code{model.mesh} is set to \code{TRUE}.}
}

\details{
When \code{object} is a list with components \code{coords} and \code{triples} (and possibly other components), the list is returned with the \code{class} assigned as \code{face3d}.

When \code{object} is a model mesh, the rownames of \code{object} are used to construct a triangulation.  A \code{face3d} object is then created, using \code{object} as the \code{coords} component.
}

\value{
A \code{face3d} object.
}

\examples{
\dontrun{
shape <- as.face3d(list(face$coords, face$triples)
}}
