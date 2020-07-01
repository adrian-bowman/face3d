\name{orient.face3d}

\alias{orient.face3d}

\title{Find the orientation of a face into a frontal view}

\description{
  This function finds the rotation which orients a face into a frontal view.
}

\usage{
orient.face3d(face, ngrid = c(16, 4), angle.max = c(pi/4, pi/8), nlevels = 1,
              graphics = FALSE)
}

\arguments{
\item{face}{a facial image in \code{face3d} form.}
\item{ngrid}{the number of angles which will be searched in each direction to identify the frontal position.}
\item{angle.max}{a vector of length 2 which controls the range of rotations which are searched
	for an optimal position.  The first component controls the range for the x and y rotations
	at the first stage.  The second component controls the range of y rotations within the y and z
	rotations at the second stage.}
\item{nlevels}{a positive integer which controls the number of levels over which the search is
	carried out.  At the moment, this is restricted to \code{nlevels = 1}, which corresponds to
	a search over x and y rotations, or \code{nlevels > 1}, when a further search over y and z
	rotations is performed.}
\item{graphics}{a logical value which controls whether grpahics are produced to monitor the
	progress of the search.}
}

\details{
This is an old desciption, which needs to be updated.

This function uses some simple information from the normal vectors at each point on the image
to orient the face into a frontal view.  If the normal vectors are not already available, these
can be created with the \code{index.face3d} function.  Angular rotations around the y and z axes
are considered, in order to maximise the asymmetry of the projection of the absolute value of the
x component of the normal vectors onto the x-y plane.  This asymmetry is measured by binning these
absolute values into a regular grid, placing a reference line at the position on the x-axis
corresponding to the mean x-position of the absolute values and then measuring the distances
between the absolute values between each bin and its reflection.  The rotations considered are
defined by a regular grid of \code{n.angle} values between \code{-max.angle} and \code{max.angle}
around the y and z axes.  Finally, rotations around the x-axis are considered, again using a
regular grid of \code{n.angle} values between \code{-max.angle} and \code{max.angle}.  The rotation
which maximise the sum of the z components of the normal vectors is selected.

This function can also be useful in identifying facial features such as the nose, eyes and mouth,
through the x, y and z components of the normal vectors when the face is oriented into a frontal view.

New material.

The output of the function is a set of paramters which define the rotation required, plus a function which can implement this rotation.  The face itself is not rotated.  There may be various components, other than \code{coords}, of the \code{face3d} object to which the rotation needs to be applied so implementation is left to the user.
}

\value{
A list with the following components:
\item{centroid}{a vector of length 3 which gives the origin around which rotation takes place.}
\item{rotation.angles}{a vector of length 3 which the angles of rotation around the x,y,z axes respectively.}
\item{orient}{a function, taking arguments \code{(x, centroid, rotation.angles)}, which will apply the rotation to the object \code{x}, which is a matrix with three columns giving the co-ordinates of the object to be rotated.}
\item{nearest}{an integer vector containing the indices of points, from horizontal strips defined by intervals in the y-axis, which are most prominent in the z-axis.  This information is used by other functions in estimating the positions of key landmarks.}

The \code{lmks} component of the \code{face} object are also rotated.  So to are the \code{curves} component if this is present.

In additional, a vector of length 3 is added to the object as the  component \code{rotation.angles} to
record the rotation angles around the y, z and x axes respectively.  These angles can be useful in
rotating the face back to its original position, if required.  However, note that in that case the 
rotations must be applied in reverse order, around the x-axis, then z-axis, then y-axis.
}

\examples{
\dontrun{
# Rotate a face away from frontal position
data(face)
face1         <- face
face1$coords  <- rotate3d(face1$coords,  pi/8, 0, 1, 0)
face1$coords  <- rotate3d(face1$coords,  pi/8, 1, 0, 0)
face1$coords  <- rotate3d(face1$coords,  pi/8, 0, 0, 1)
face1$normals <- rotate3d(face1$normals, pi/8, 0, 1, 0)
face1$normals <- rotate3d(face1$normals, pi/8, 0, 0, 1)
face1$normals <- rotate3d(face1$normals, pi/8, 1, 0, 0)
plot(face1)

# Identify the rotation required to produce a frontal view
face1 <- orient.face3d(face1, verbose = 1)
plot(face1)

# Now rotate back to the original position
rot           <- face1$rotation.angles
face1$coords  <- rotate3d(face1$coords,  -rot[3], 1, 0, 0)
face1$coords  <- rotate3d(face1$coords,  -rot[2], 0, 0, 1)
face1$coords  <- rotate3d(face1$coords,  -rot[1], 0, 1, 0)
face1$normals <- rotate3d(face1$normals, -rot[3], 1, 0, 0)
face1$normals <- rotate3d(face1$normals, -rot[2], 0, 0, 1)
face1$normals <- rotate3d(face1$normals, -rot[1], 0, 1, 0)
plot(face1)
}
}
