\name{face}

\alias{face}

\title{A face shape}

\description{
The object consists of a list with three components.  The component \code{coords} is a matrix with three columns containing the three-dimensional co-ordinates of the points.  The component \code{triples} is a vector of integers containg the indices of \code{coords} which should be passed to the \pkg{rgl} function \code{triangles} to render the surface.  The component \code{colour} contains the colours of the points corresponding to the rows of \code{coords}.
}

\examples{
\dontrun{
   data(face)
   display.face3d(face)
}}

\keyword{iplot}
\keyword{dynamic}

