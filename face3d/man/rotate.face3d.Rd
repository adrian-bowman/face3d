\name{rotate.face3d}

\alias{rotate.face3d}

\title{Rotate a face3c object or a matrix.}

\description{
This function rotates a \code{face3d} object or a matrix of 3D coordinates.  The rotation can be specified by and angle, rotation axis and centre, or by planes specified by landmarks: coronal (frontal), sagittal sinister (lateral sinister), sagittal dexter (lateral dexter), and transversal (vertical/basal).
}

\usage{
rotate.face3d(shape, angle, raxis, centre = rep(0, 3), center,
              landmarks = c("sn", "se", "exR", "exL"), rotation = "coronal",
              exclude = NULL)
}

\arguments{
\item{shape}{a \code{face3d} object or matrix of 3D coordinates to be rotated.}
\item{angle}{a numeric value specifying the angle of rotation in radians.}
\item{raxis}{a numeric vector of length 3 defining the axis of rotation.}
\item{centre, center}{a numeric vector of length 3 defining the origin around which rotation will take place.  \code{center} is a synonym for \code{centre}.}
\item{landmark.names}{a character vector of length 4, with the first two entries giving mid-sagittal landmarks which define the vertical direction (bottom to top) and the second two entries giving paired landmarks which define the horizontal direction (right to left); the default is \code{c("sn","se","exR","exL")}.  The landmarks are identified from the named rows of the three-column matrix supplied in \code{landmarks} or, if this is not supplied, \code{shape$landmarks}.  The function will fail if neither of these is present.}
\item{rotation}{\code{"coronal"}, \code{"sagittal sinister"}, \code{"sagittal dexter"} or \code{"transversal"}.}
\item{exclude}{a character vector giving the names of components of \code{shape}, when it is a \code{face3d} object, which should not be rotated.}
}

\details{
When \code{shape} is a \code{face3d} object, the rotation is applied by translating the relevant component from \code{centre} to the origin, applying the rotation and the translating back to \code{centre}.  The \code{triangles} component is always excluded.
}

\value{
If \code{shape} is a \code{face3d} object, then a \code{face3d} is returned, with the appropriate components rotated.   If \code{shape} is a matrix, a list is returned with components \code{ROT.shape} containing the rotated shape, and \code{angles} containing the angles of rotation.
}
