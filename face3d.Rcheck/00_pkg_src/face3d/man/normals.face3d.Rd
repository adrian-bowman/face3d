\name{normals.face3d}

\alias{normals.face3d}

\title{Compute the normal vectors at each vertex in a face3d object.}

\description{
  This function computes the normal vectors at each vertex in a face3d object.
}

\usage{
normals.face3d(shape)
}

\arguments{
\item{shape}{the face3d object.}
}

\value{
  A new version of the face3d object, with the normals added, is returned.
}

\examples{
\dontrun{
   data(face)
   normals <- normals.face3d(face)$normals
}}
