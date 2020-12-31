\name{editlandmarks.face3d}

\alias{editlandmarks.face3d}

\title{Interactive editing of the positions of landmarks}

\description{
  This function display the shape in an \pkg{rgl} window and provides an interactive control panel which allows landmark positions to be adjusted.
}

\usage{
editlandmarks.face3d(shape, lmk.names = rownames(shape$lmks), panel = TRUE,
                     lmk.name = "none", directions = FALSE)
}

\arguments{
\item{shape}{the face3d object.}
\item{lmk.names}{a character vector giving the names of the landmarks to be edited.  If new landmarks are required then their names should be added to the \code{lmk.names} vector.}
\item{panel}{a logical value which determines whether a control panel is created for interactive use.}
\item{lmk.name}{a text string identifying the landmark to be used for current editing.}
\item{directions}{a logical value controlling whether axis directions are shown, to indicate the current movement of the identified landmark.}
}

\details{
The function launches an \pkg{rgl} window to display the \code{shape} and a control panel (created using the \pkg{rpanel} package) which allows interactive specification and movement of landmarks.  Landmarks are specified by name through a combo box.  A 'Zoom in/out' button provides the option to zoom in on the area immediately surrounding the specified landmark.

A 'Show path' checkbox allows a path across the shape to be created and displayed, using the \code{planepath.face3d} function.  With large shapes, it can be useful to activate this checkbox only when zoom is activated.  This is because the calculation of a path across a very large object can be a little slow.

The direction of the path across the object is controlled by an 'angle' doublebutton while the movement of the specified landmark along this path is controlled by a distance doublebutton.  A 'Perpendicular direction' button provides a convenient mechanism for switching the current direction of travel by 90 degrees.

Radiobuttons are used to control the type of display (points, triangular mesh or rendered surface) and the colour (natural texture or shape index).

If a specified landmark is not present in the \code{shape} object passed to the function then a new landmark is created.  Its initial position is approximately defined by the closest point which lies along the line of view to the origin of rotation of the shape.  This means that new landmarks can be placed in approximately the correct initial position by rotating the shape to put the desired location approximately in the centre of the field of view before selecting the landmark name.  The zoom facility can also be useful here.

The 'Save shape' button must be pressed if the edited landmarks are to be inserted into the shape object returned by the function.
}

\value{
The entire \code{shape} object is returned by the function.  If the "Save shape" button is pressed before the control panel is closed and the function call ends, then the current state of the landmarks is stored in \code{shape$lmks}.  Note that it is the state of the landmarks at the time of the button press which is stored.  The button should therefore be pressed immediately before the control panel is closed, if the most recent landmarks are to be returned.
}

\examples{
\dontrun{
   data(face)
   editlandmarks.face3d(face)
}}
