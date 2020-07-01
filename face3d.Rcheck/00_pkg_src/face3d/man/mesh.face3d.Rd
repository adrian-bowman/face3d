\name{mesh.face3d}

\alias{mesh.face3d}

\title{Find a planar path between anatomicaly pre-specified set of curve pairs along the surface of a face3d object}

\description{
This function finds a planar path between anatomicaly pre-specified set of curve pairs along the surface of a face3d object.
}

\usage{
mesh.face3d(face, meshes, template, graphics=FALSE, monitor=FALSE)
}

\arguments{
\item{face}{A \code{face3d} object.}
\item{meshes}{the specified mesh patch to find.}
\item{template}{a template used in the mesh construction.}
\item{graphics}{a logical value which controls whether diagnostic graphics are created as the function executes.}
\item{monitor}{a logical value which controls whether progress is monitored as the function executes.}
}

\details{
A plane path (curve) between anatomicaly pre-specified set of 20 curve pairs as follows: mid-face R/L (cheek-lip R/L - cheek-nose R/L), upper mid face R/L (cheek-nose R/L - cheek-eye R/L), upper lip R/L (upper lip R/L - mid-line lip R/L), lower lip R/L (lower lip R/L - mid-line lip R/L), philtrum R/L (upper lip R/L - nasal base R/L), lower face R/L (mandible R/L - lower lip R/L), nose R/L (nasal bridge R/L - nasal root R/L), upper face R/L 1 (brow ridge R/L - cheek-eye R/L), upper face R/L 2 (brow ridge R/L - upper eye socket R/L), upper face R/L 3 (brow ridge R/L- nasal root R/L).
}

\value{
Returning \code{shape} as \code{face3d} object with a set of meshes as a list concatenated to a matrix \code{shape$meshes}.
}

\examples{
\dontrun{
   data(face)
   face <- mesh.face3d(face, "mid-face right")
}}
