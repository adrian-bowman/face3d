\name{findlandmarks.face3d}

\alias{findlandmarks.face3d}

\title{Automatically locate landmarks on a facial surface.}

\description{
This function finds a small number of anatomical landmarks on a face, in an automatic manner.
}

\usage{
findlandmarks.face3d(face, lmks = c("pn", "enL", "enR", "se"),
                     orient = TRUE, monitor = FALSE, graphics = FALSE)
}

\arguments{
\item{face}{A \code{face3d} object representing a facial surface.}
\item{lmks}{The anatomical landmarks to be located.  The four listed in the default setting are the only ones which are currently available.}
\item{orient}{a logical value which controls whether the face is oriented, using the function \code{orient.face3d}, before the landmarks are sought.}
\item{monitor}{a logical value which controls whether text is printed to monitor the progress of the algorithm.}
\item{graphics}{a logical value which determines whether graphics are produced to onitor the progress of the algorithm.}
}

\details{
The landmark pn is located by starting with the closest point to the viewer, among those identified (\code{nearest}) in the output of the \code{orient.face3d} function.  The curvature values of the neighbouring points are then used to find the weighted average of the mesh points among those with suitable large shape index.  The pn landmark is then taken to be the closest mesh point to this mean value.

The enR/L landmarks are sought within the two largest connected sets of points with small shape index.  The landmarks are taken to be the points whose minimum curvature is maximum within these connected sets.

The se landmark is located from a curve between the enL/R landmarks using the \code{planepath.face3d} function based on minimising distance.  The landmark is then taken to be the point of maximum curvature along this path, identified through the \code{gcurvature.face3d} function.
}

\value{
THe function returns the \code{face3d} object with the requested landmarks placed inside the \code{lmks} component.
}

\examples{
\dontrun{
   face <- findlandmarks.face3d(face)
}}
