\name{render_color.face3d}

\alias{render_color.face3d}

\title{Creates a new \code{face3d} object rendered by topographic color differences from an original set of landmarks to a target set of landmarks. }

\description{
  This function creates a new \code{face3d} object rendered by topographic color differences from an original set of landmarks to a target set of landmarks.}

\usage{
 render_color.face3d(face, reference, target, landmark.names,
                     nbrk = c(-30, seq(-10,-1,1), seq(1,10,1), 30))}

\arguments{
\item{face}{The original \code{face3d} object of the reference landmarks.}
\item{reference}{A n by 3 mesh matrix which is the the main reference set of landmarks that will be coloured.}
\item{target}{A n by 3 mesh matrix with is the target set of landmarks.}
\item{landmark.names}{The rownames of the reference and target landmark sets}
\item{nbrk}{Defines the topographic colouring breaks.}
}

\details{
The rendered reference \code{face3d} object and the target are produced which are coloured by differences}

\value{
  A list containing the two elements \code{rend.reference} and \code{rend.target}.
}

\examples{
\dontrun{
   data(face)
}}
