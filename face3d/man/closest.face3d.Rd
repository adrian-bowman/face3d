\name{closest.face3d}

\alias{closest.face3d}

\title{Find the closest points in a face3d object}

\description{
This function finds the closest point to each row of \code{x} in the face3d object \code{shape}.
}

\usage{
closest.face3d(x, shape, nearest = 30, tol = 1e-8)
}

\arguments{
\item{x}{a vector of length 3 or a matrix with 3 columns, giving the co-ordinate position(s) of the points in three-dimensional space whose closest points on \code{shape} are to be located.}
\item{shape}{a \code{face3d} object.}
\item{nearest}{a positive integer controlling the number of nearest neighbours which will be used to define the neighbourhood within which the closest point is sought.}
\item{tol}{the tolerance used in deciding whether a point lies on an edge.}
}

\details{
A subset of \code{shape} is defined, using the \code{subset.face3d} function, through the  \code{nearest} points within \code{shape$coords} which are closest to \code{x}.  This neighbourhood is then searched to identify the closest point to \code{x}.  This closest point may be a vertex of \code{shape}, it may lie on an edge or it may be an interior point of a triangle within the mesh.

The returned list has a component \code{ids} which contains the indices of \code{shape} which identify the vertices which surround the closest points.

Rows of \code{x} which contain missing values are removed. In this case the components of the object returned by the function refer to the remaining rows of \code{x}.
}

\value{
A list with three components: 
\item{points}{a matrix with three columns giving the co-ordinates of the closest points. If \code{x} is a vector then this component is also returned as a vector.}
\item{ids}{a matrix with three columns where each row gives the indices of the rows of \code{shape$coords} which encompass the closest point.  Where the corresponding row of \code{point} is a vertex this is a single index, where it lies on an edge there are two indices to define the end points, and where it is an interior point of a triangle there will be three indices for the triangle vertices.  Each row of \code{ids} is padded with \code{NA} as required, to allow it to be returned as a matrix.  If the input argument \code{x} corresponds to a single point then \code{NA} padding is not used.}
\item{distances}{the distance of each closest \code{point} from the corresponding row of \code{x}.}
}

\examples{
\dontrun{
# Projection onto vertices, edges and faces
shape   <- as.face3d(list(coords = matrix(c(0, 0, 0, 1, 0, 0, 0.5, 1, 0, 0.5, 0.5, 1),
                              ncol = 3, byrow = TRUE),
                        triples = c(1, 2, 4, 2, 3, 4, 1, 4, 3, 3, 2, 1)))
pts     <- matrix(c(1.2, -0.2, 0, 0.5, 1.2, 0.5, 0.3, 0.5, 0, 0, 0, 0, 0.5, 0, 0, 0.5, 0, 0.5),
                  ncol = 3, byrow = TRUE)
closest <- t(apply(pts, 1, function(x) closest.face3d(x, shape)$point))
plot(shape)
spheres3d(pts,     radius = 0.02, col = "red")
spheres3d(closest, radius = 0.03, col = "green")
segments3d(matrix(c(t(cbind(pts, closest))), ncol = 3, byrow = TRUE))

# Projection onto the nose of template_male
x    <- template_male$lmks["pn", ] + c(0.2, 0.1, 1)
sbst <- subset(template_male, rdist(t(x), template_male$coords) < 10)
cx   <- closest.face3d(x, sbst)
plot(sbst, display = "mesh")
spheres3d(x, col = "red", radius = 0.1)
spheres3d(cx$point, col = "blue", radius = 0.1)
spheres3d(sbst$coords[cx$id, ], col = "green", radius = 0.1)
}}
