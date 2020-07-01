\name{smoothpath.face3d}

\alias{smoothpath.face3d}

\title{Create a smooth path through 2d point space using different methods of choice.}

\description{
  This function creates a smooth path through 2d point space using different methods of choice.
}

\usage{
smoothpath.face3d(shape, shape.smooth, lmk1, lmk2, df = 5, distance = 5,
                         curve.name = NA, penalty = 1,
                         reference.path = NA, rotation.range = rotation.range,
                         si.target,
                         boundary = c(0.2, 0.5), npts, graphics = FALSE,
                         perp.dist.bound = 10, pcrv.path = FALSE, start.pc = NA,
                         smooth.method = NULL, average.alpha = average.alpha)
}

\arguments{
\item{shape}{The subsetted shape from the \code{face3d} object of interest.}
\item{shape.smooth}{The optional more subsetted shape from the subsetted shape of the \code{face3d} object of interest. }
\item{lmk1}{A vector anatomical landmark that the curve must start at.}
\item{lmk2}{A vector anatomical landmark that the curve must finish at.}
\item{df}{A numerical value determining the degrees of freedom for smoothing.}
\item{distance}{A numerical value determining the distance for the shape index to look over.}
\item{curve.name}{The curve name from the list given in \code{facecurves.face3d} determined by \code{lmk1} and \code{lmk2}.}
\item{penalty}{Description required.}
\item{reference.path}{A different 3 dimensional matrix reference curve path than that created by \code{planepath.face3d} betweeen \code{lmk1} and \code{lmk2}.}
\item{rotation.range}{A value between 0 and \code{pi/2} for \code{planepath.face3d} to test over.}
\item{si.target}{A numerical value determining the curvature value of interest chosen to be \code{0.5} or \code{-0.5}.}
\item{boundary}{a numerical vector that determines the percentage boundary expansion past \code{lmk1} and \code{lmk2} that determines the subsetted shape surface.}
\item{npts}{ A numerical value determining the amount the curve will be resampled to when necessary}
\item{graphics}{a logical value which controls whether diagnostic graphics are produced as the function executes.}
\item{perp.dist.bound}{A numerical value determining the y-axis cut off.}
\item{pcrv.path}{a logical value which controls whether a principal curve or a planepath is used as a reference path.}
\item{start.pc}{an starting path used to initialise the principal curve algorithm.}
\item{smooth.method}{a text string which controls the type of smoothing used.  The options are quadratic (\code{"quad"}) and principal curve (\code{"pcurve"}) smoothing.}
\item{average.alpha}{Description required.}
}

\details{
A smoothed path in three dimensions and its associated curvature values will be returned.
}

\value{
  A list containing the component \code{path} and \code{values}.
}


