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
  The function returns a new version of the \code{shape} with the additional components: \item{normals}{a matrix whose rows contain the normal at each vertex;}
\item{axes}{a list of 3x3 matrices, one for each vertex.  The columns of each matrix give the normal vector and two arbitrary axes in the tangent plane at each vertex.}
}

\examples{
\dontrun{
   normals <- normals.face3d(template_male)
}}
