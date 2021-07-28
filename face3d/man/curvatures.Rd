\name{curvatures}

\alias{curvatures}

\title{Compute the curvatures across a face3d object}

\description{
  This function constructs the principal curvatures and directions at each vertex of a \code{face3d} object.
}

\usage{
curvatures(shape, distance = 10, sample.spacing,
             type = "euclidean", subset = 1:nrow(shape$coords),
             extension = "NA", overwrite = FALSE, directions = FALSE, monitor = 0)
}

\arguments{
\item{shape}{a face3d object.}
\item{distance}{the size of the radial neighbourhood which will be used to calculate the curvatures at each location.}
\item{subset}{a vector which defines a subset of \code{shape$coords} on which the curvatures will be calculated.  This may be in the form of a set of indices, or a logical vector, defining a set of vertices.  It may also be a number between 0 and 1, which will be used to define a random proportion of vertices.  If \code{subset} is missing, all co-ordinates will be used.}
\item{interpolate}{a logical value which controls whether the curvatures evaluated on a subset are interpolated to all vertices, using the \code{warp.face3d} function.  if the \code{subset} argument is not used then the \code{interpolate} argument is ignored.}
\item{extension}{a logical value or vector of indices which define a slightly larger area of the object than \code{subset}.  This allows the patches around particular points to extend beyond \code{subset}, hence improving the estimation of curvature at the edges.  When \code{extension} is logical and \code{TRUE}, all points within \code{distance} of any point in \code{shape$coords} will be included in the curvature calculations.  When \code{extension} is an integer vector, the extension is provided by the vertices of \code{shape} which have these indices.  If this argument is missing, curvature is calculated using only points defined by \code{subset}, except in the case where the number of points in \code{sbst} is less than 10, when \code{extension} is set to \code{TRUE}.}
\item{overwrite}{a logical value determining whether existing values of curvature will be overwritten by new ones.  It will clearly be faster not to repeat the caluclation of indices if the same value of \code{distance} is used.}
\item{directions}{a logical value which controls whether the principal directions are placed in the \code{face3d} object returned by the function.}
\item{monitor}{an integer value which controls the level of reporting.  If \code{monitor} is > 0 and some vertices have insufficient (< 8) neighbours to compute the local curvature then a message reporting the number of these vertices is printed.}
}

\details{
For each point on the mesh which defines the face3d object, the neighbours are identified.  This is done either by identyfing connected points and then subsequent connected points (\code{type = "mesh"}), to a degree specified by the \code{extent} argument, or using all points within a distance (\code{type = "euclidean"}) determined by the \code{distance} argument.  A `normal' surface direction at each point is identified by averaging the normals for each of the (1st) neighbouring triangles.  A local co-ordinate system is then contructed, using the normals as the `'response' directions, and a cubic regression surface fitted to the neighbours.  The principal curvatures are then computed, from which the shape indices can be found.

A common use of subset is through expressions such as \code{subset = (dst < 70)}, where \code{dst} denotes the distances of each vertex in the \code{face3d} object rom some point of interest.  In this case, the \code{extension} argument can be conveniently set as \code{extension = (dst < 70 + distance)}.
}

\value{
An augmented version of the input \code{shape} with the additional components:
\item{kappa1, kappa2, shape.index}{vectors containing the the maximal and minimal curvature and the shape index at each vertex.}
\item{normals}{a matrix whose rows give the co-ordinate direction of the estimated surface normal at each vertex, if this was not already present.}
\item{distance}{the distance value used to compute the curvatures.}
}

\examples{
\dontrun{
nose <- subset(template_male, edist(template_male$landmarks["pn", ], template_male$vertices) < 30)
nose <- curvatures(nose, distance = 10)
plot(nose, col = "shape index")
}
}
