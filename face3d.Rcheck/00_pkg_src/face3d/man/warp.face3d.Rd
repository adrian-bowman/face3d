\name{warp.face3d}

\alias{warp.face3d}

\title{Construct a function which reproduces specified information at particular 3D locations exactly and which allows interpolation at further 3D locations.}

\description{
This function creates a transformation which reproduces specified information at particular 3D locations exactly and which allows interpolation at further 3D locations.  This can be used to warp a high resolution 3D template onto a lower resolution image, using landmarks or specified points of shape correspondence.  It can also be used to interpolate scalar information, such as curvature values, across an image, using information available at a subset of locations.
}

\usage{
warp.face3d(from, to, carry)
}

\arguments{
\item{from}{a three-column matrix supplying the 3D locations of the points to be mapped exactly.}
\item{to}{a matrix or vector defining the information to be mapped exactly.}
\item{carry}{a matrix of 3-d co-ordinates, or a list of such matrices, or a \code{face3d} object.  In each case, the results of applying the warp function to this additional shape information are returned.  The Value section below explains how the different forms of this argument are handled.}
}

\details{
The warp function is constructed as transformation of 3-d space which maps the points in \code{from} to the points in \code{to}. The transformation is based on a thin-plate splines whose coefficients then allow the transformation to be applied to any other 3-d locations.

Note that when \code{carry} is a \code{face3d} object and \code{to} is a matrix with 3 columns it is assumed that the warping function is operating on 3D co-ordinates.  If this is not the case, and \code{to} refers to different information, then \code{carry} should be set to a matrix of co-ordinates, not a \code{face3d} object.
}

\value{
There are three cases, depending on the form of the input in the \code{carry} argument of the function.
\itemize{
\item
When \code{carry} is a single matrix or vector, the returned object is a single transformed matrix or vector.
\item
When \code{carry} is a list of matrices or vectors, the returned object is list of corresponding transformed matrices or vectors.
\item
When \code{carry} is a \code{face3d} object, the transformation is applied to \code{coords} and any of the components \code{lmks}, \code{curves}, \code{mesh} which are present.  A new \code{face3d} object containing these transformed components, and the \code{triples} component which defines the triangulation, is then returned.  Note that no other components of the original \code{face3d} object are carried over to the new \code{face3d} object.
}
}
