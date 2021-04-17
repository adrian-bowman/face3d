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
If the special value \code{"none"} is given in the \code{landmark.names} argument then the returned value contains the shape index and curvatures computed from the value of \code{distance}.  THis can be useful in troubleshooting.

The setting \code{monitor = 3} allows the stages of the process to be inspected by rotating the images in the \pkg{rgl} window.  This can indicate where the problem lies and suggest alterations to the defaukt setting of the distance parameters.

gn is estimated by subsetting the low resolution to include only areas where the Gaussian curvature is above the 80th percentile of the positive values.  Points which are within a distance 20 of pn are also excluded.  At each point in the remaining object, the Gaussian curvature, the cosine of the angle between the vector from pn to the point and the vector from se to pn, and an indicator of whether the distance of the point from pn lies between 0.5 and 2 times the distance between se and pn, are all multiplied together to create an index.  This is then integrated over each distanct part of teh object.  The part with the highest integral is selected and gn is estimated as the local model using a distance of 15.
}

\value{
A \code{face3d} object which has estimated positions of the nominated landmarks in the names rows of the \code{landmark} component.

Note that any existing information on curvature (\code{shape.index}, \code{kappa1}, \code{kappa2}) will be overwritten by the operations of the function.
}

\examples{
\dontrun{

}}
