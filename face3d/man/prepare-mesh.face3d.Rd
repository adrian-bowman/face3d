\name{prepare_mesh.face3d}

\alias{prepare_mesh.face3d}

\title{Preparing a full mesh for rendering.}

\description{
  This function takes a full facial mesh and prepares it for rendering. Can be done directly after the sliding stage.
}

\usage{
prepare_mesh.face3d(face)}

\arguments{
\item{face}{A \code{face3d} object with a full facial mesh that has already been slided.}
}
\details{
A \code{face3d} object containing a prepared mesh for rendering.
}

\value{
  None
}

\examples{
\dontrun{
   data(face)
face <- prepare-mesh.face3d(face)
}}
