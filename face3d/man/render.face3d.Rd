\name{render.face3d}

\alias{render.face3d}

\title{Creates a \code{face3d} object from a set of mesh landmarks that have been slid and prepared.}

\description{
  This function creates a \code{face3d} object from a set of mesh landmarks that have been slid and prepared.
}

\usage{
render.face3d(rend.object, subset = NULL)}

\arguments{
\item{rend.object}{A matrix of a set of mesh landmarks to be rendered.}
\item{subset}{a subset of of points.  This argument is not currently used.}

}

\details{
A \code{face3d} object is returned of the input landmark points \code{rend.object}.
}

\value{
  An object called \code{rendered.object} of class \code{"face3d"}.
}

\examples{
\dontrun{
   data(face)
face <- render.face3d(mesh.pts)
}}
