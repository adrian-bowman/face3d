\name{planepath}

\alias{planepath}

\title{Find a planar path through the mesh of a face3d object}

\description{
This function finds a planar path through the mesh of a face3d object.  The plane is defined by two conditions.  The first is that it passes through a specified point, the second is that it contains either a second point or a specified direction.
}

\usage{
planepath(shape, x1, x2, pts1, pts2, direction, normal, bothways = FALSE,
                         distance = 5, boundary = c(0.2, 0.5),
                         rotation = "optimise", rotation.range = pi/2,
                         ngrid = 50, bridge.gaps = FALSE, si.target, directions = FALSE,
                         monitor = FALSE, monitor.prompt = monitor)
}

\arguments{
\item{shape}{a face3d object.}
\item{x1}{a point which the cutting plane must contain.}
\item{x2}{an optional second point which the cutting plane must contain.}
\item{pts1, pts2}{optional vectors which identifies the rows of \code{shape$coords} which are associated with \code{x1} and, if present, \code{x2}.  See the help file for \code{closest.face3d} for the length and meaning of this.}
\item{direction}{a vector of length 3 defining a direction from \code{x1} which the plane will contain if \code{x2} is not specified.}
\item{normal}{a vector of length 3 defining a direction from \code{x1} which the plane will contain, or be rotated away from.  If this is not specified then the normal will be computed from the local surface at \code{x1}, and \code{x2} if supplied.}
\item{bothways}{a logical value which determines whether the path is drawn in both directions when only \code{x1} is specified.}
\item{distance}{the distance used in the call to the function \code{curvatures} if \code{si.target} is not missing, so the shape index is required, but this is not already present.}
\item{boundary}{this can be set to a vector of length 2 which can be used to restrict the search area within the object.  When \code{x1} and \code{x2} are both present, a numeric vector of length 2, interpreted as multipliers of the distance between \code{x1} and \code{x2}, is expected.  When only \code{x1} is present, a single numeric value is expected.  See the Details section below for how these values are used and for the default behaviour.}
\item{rotation}{if a numerical value, or a vector of numerical values, is supplied then this is taken to be an angle(s) by which the cutting plane containing \code{x1} and \code{x2}, or from \code{x1} in the specified \code{direction}, is rotated away from the normal at \code{x1}.  If the character values \code{"optimize"} or \code{"optimise"} is supplied then a range of angles between (but not including) \code{-rotation.range} and \code{rotation.range} are used and the path with smallest length is returned.  The character values \code{"maxmise"} or \code{"maxmize"} is supplied then the path with largest length is returned.}
\item{rotation.range}{the distance covered, in both positive and negative directions, by the range of rotation values when optimising over rotation.}
\item{ngrid}{the number of values used to define the range of rotation angle.}
\item{bridge.gaps}{a logical argument which controls whether the identified path will cross gaps in the mesh or not.}
\item{si.target}{a numerical value which, if present, determines whether shape index is used to guide the planepath.  Currently only the sign of \code{si.target} is used.}
\item{directions}{a logical argument which controls whether the directions of principal curvature are stored, if shape index or curvature information is computed.}
\item{monitor}{an integer value which determines whether graphics are displayed as the algorithm progresses.  A value greater than 0 initiates feedback, with the values \code{1} or \code{2} control the level of feedback supplied.}

}

\details{
The plane is defined by two conditions.  The first is that it passes through a specified point, the second is that it contains either a second point or a specified direction.

If the object is large, with many vertices, then it can be very inefficient to search for paths in the entire object.  A better strategy is to subset the object appropriately and the \code{boundary} argument allows this to be controlled.  There are two cases to consider.
\enumerate{
\item
If \code{x2} is specified, two constraints are applied.  The first includes only vertices whose projections onto the line through \code{x1} and \code{x2} lies between these two points, or whose distance from \code{x1} or \code{x2} is no more than the distance between \code{x1} and \code{x2} multiplied by \code{boundary[1]}. The second constraint includes only vertices whose distance from the cutting plane is no more than the distance between \code{x1} to \code{x2} multiplied by \code{boundary[2]}.  The default setting is \code{boundary = c(0.2, 0.5)}.  If any element of boundary is set to \code{NA} then the whole shape will be used.
\item
If \code{x2} is not specified but \code{direction} is, then vertices within a distance \code{boundary[1]} of \code{x1} are included.  In this case the default setting is to include the whole shape.
}
In either case, if the boundary constraints lead to a subset with fewer than \code{3} vertices then a warning is issued and the whole shape is used.

The arguments \code{pts1, pts2} allow the positions of \code{x1} and \code{x2} within the mesh to be specified.  Each of these arguments is an integer vector which identifies the points, as rows of \code{shape$coords}, which encompass the point of interest.  Where this point is a vertex then the vector will be of length 1; where the point lies on an edge the vector will be of length 2 to define the end points; where the point is an interior point of a triangle the vector will be of length 3 giving the indices of the triangle vertices.  The returned values of \code{pts1} and \code{pts2} have the same interpretation.
}

\value{
A list with the following components:
\item{path}{a matrix of ordered points at which the plane cuts the mesh.  If \code{x1} and \code{x2} are both specified then these points are also added in the top and bottom rows.  The ordering is from \code{x1} to \code{x2}.}
\item{arclength}{a vector which gives the cumulative length of the identified path along its successive points.}
\item{criterion}{the minimum value of the criterion used to identify the path.}
\item{shape}{the subset of the input shape used to identify the path.}
\item{pts1, pts2}{the index of the rows in \code{shape$coords} which surround \code{x1} and \code{x2}.  Note that this refers to the returned \code{shape} which is a subset of the \code{shape} which passed as an argument to the \code{planepath} function.}
\item{angle}{the angle away from the normal associated with the identified path.}
\item{normal}{the normal vector which defines the plane which corresponds to \code{rotation = 0}.}
When the \code{directions} argument is set to \code{TRUE}, the following additional components are added to \code{shape}:
\item{directions}{the two principal directions and the normal direction at the points along the path.}
\item{kappa1, kappa2, shape.index}{the value of \code{kappa1},\code{kappa2} ad \code{shape.index} at the points along the path.}
\item{triangles}{when \code{x1} and \code{x2} are both specified, the triangles which lie along the path are returned.  Specifically, \code{triangles} is a matrix with 3 columns, with each row giving the indices of the points in \code{shape} which form the vertices of a particular triangle.  Note that these indices refer to points in the shape which was originally passed to the function, not in the subset shape returned by the function.}
If, for any reason, a path cannot be identified, a message is printed to that effect and the returned value is \code{NULL}.
}

\examples{
\dontrun{
data(face)
nose <- subset.face3d(face, face$coords[ , 3] > 85 & 
           face$coords[ , 2] > -60 & face$coords[ , 2] < -15)
plot(nose, type = "mesh")
x1 <- nose$coords[2737, ]
x2 <- nose$coords[146, ]
spheres3d(t(matrix(x1)), col = "red", radius = 0.2)
spheres3d(t(matrix(x2)), col = "red", radius = 0.2)
p <- planepath(nose, x1, x2)
lines3d(p$path, col = "green")
}
}
