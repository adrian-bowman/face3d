\name{slidingmeshes.face3d}

\alias{slidingmeshes.face3d}

\title{Slides the meshes created using \code{mesh.face3d} to minimize bending energy between the found mesh patches and the template.}

\description{
  This function slides the mesh patches created using \code{mesh.face3d} to minimize bending energy between found mesh patches and the template.
}

\usage{
slidingmeshes.face3d(shape1,shape2, n = 100)
}

\arguments{
\item{shape1}{The template.}
\item{shape2}{A \code{face3d} objects containing found meshes from \code{mesh.face3d}.}
\item{n}{The multiplicative number of points for resampling.}

}

\details{
A matrix of slid and resampled mesh points is returned.
}

\value{
A matrix of landmark mesh points \code{curves2.optim} that have been slid and resampled.
}

\examples{
\dontrun{
   data(face)
   meshes.new  <- slidingmeshes.face3d(template, face)
}}
