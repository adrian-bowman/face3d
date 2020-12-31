\name{approximatelandmarks.face3d}

\alias{approximatelandmarks.face3d}

\title{Estimate approximate landmark positions on a human face}

\description{
Approximate landmark positions are estimated on a human face.
}

\usage{
approximatelandmarks.face3d(face, landmark.names = c("pn", "enL", "enR", "se"),
                                   threshold, trim = 3 * threshold, distance = 5 * threshold,
                                   monitor = 1, overwrite = FALSE)

}

\arguments{
\item{face}{a \code{face3d} object contain a 3D image of a human face.}
\item{threshold}{a numeric value which sets the upper bound on the distance from any vertex of \code{face$coords} to the closest of the sampled points.}
}

\details{
If the special value \code{"none"} is given in the \code{landmark.names} argument then the returned value contains the shape index and curvatures computed from the value of \code{distance}.

The setting \code{monitor = 3} allows the stages of the process to be inspected by rotating the images in the \pkg{rgl} window.  This can indicate where the problem lies and suggest alterations to the defaukt setting of the distance parameters.
}

\value{
A \code{face3d} object which has estimated positions of the nominated landmarks in the names rows of the \code{landmark} component.

Note that any existing information on curvature (\code{shape.index}, \code{kappa1}, \code{kappa2}) will be overwritten by the operations of the function.
}

\examples{
\dontrun{

}}
