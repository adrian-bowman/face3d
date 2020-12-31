\name{rotate.face3d}

\alias{rotate.face3d}

\title{Rotate the matrix of landmarks or face3d object to pre-specified planes based on carefully chosen landmarks.}

\description{
This function the matrix of landmarks or \code{face3d} object to pre-specified planes - coronal (frontal), sagittal sinister (lateral sinister), sagittal dexter (lateral dexter), and transversal (vertical/basal).
}

\usage{
rotate.face3d(shape, id.lndms, rotation)}

\arguments{
\item{shape}{a matrix of landmark coordinates or \code{face3d} object in the form of a list.}
\item{id.lndms}{ids of two mid-sagittal landmarks (direction upwards) and two paired landmarks (direction right-left); the default is \code{c("sn","n","exR","exL")}.}
\item{rotation}{\code{"coronal"}, \code{"sagittal sinister"}, \code{"sagittal dexter"} or \code{"transversal"}.}
}

\details{
Rotated shape as a matrix (if shape is a matrix) or a list (if shape is an \code{face3d} object).
}

\value{
If \code{shape} is a \code{face3d} object, then a \code{face3d} is returned, with the components \code{coords}, and if present \code{lmks}, \code{curves} and \code{mesh}, rotated; all other components are left unchanged.   If \code{shape} is a matrix, a list is returned with components \code{ROT.shape} containing the rotated shape, and \code{angles} containing the angles of rotation.
}
