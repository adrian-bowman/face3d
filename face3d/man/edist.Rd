\name{edist}

\alias{edist}

\title{Euclidean-based distances between 3D co-ordinates}

\description{
The separation between two sets of 3d co-ordinates are computed, based on Euclidean distance.  The Details section below provides a description of the different distances available for different types of input object.
}

\usage{
edist(x1, x2, type = "pairwise", index = FALSE, x2curve = FALSE)
}

\arguments{
\item{x1, x2}{each of \code{x1} and {x2} may be either a vector of 3d co-ordinates or a three-column matrix whose rows are 3d co-ordinates.  One or both of \code{x1} or \code{x2} may be a \code{face3d} object.  \code{x2} may be omitted.}
\item{type}{a character string which determines which type of distance is computed. The available options are \code{"pairwise"}, \code{"min"} and \code{"minsum"}. The distances computed in each case are described in the Values section below.}
\item{index}{a logical argument which, in the case of \code{type = "min"}, determines whether an index is returned to identify which point in \code{x2} is closest to each point of \code{x1}.}
\item{x2curve}{a logical value which identifies \code{x2} as a curve or set of curves.  In this case distances will be signed when \code{type = "min"}.}
}

\details{
When one of \code{x1} or \code{x2} is a \code{face3d} object then the distances are calculated from or to the \code{vertices} component of this object.

When \code{x1} is a \code{face3d} object and \code{x2curve = TRUE}, it is assumed that \code{x2} is a matrix giving equally-spaced locations along a curve or set of curves. If \code{x2} has no rownames then it is assumed to be a single curve.  When rownames are present, these must follow a particular structure.  Any digits appearing at the end of the rownames are assumed to index points within a particular curve and these digits are removed.  The remaining character strings are used to identify individual curves.  The 3d co-ordinates for each curve are assumed to identify sequential, equally-spaced points along the curve.

The pairwise Euclidean distances which provide the starting point for all forms of distance available in the function are computed by calling the \code{rdist} function from the \pkg{fields} package.
}

\value{
When \code{type = "pairwise"}, a matrix of dimension \code{nrow(x1)} x \code{nrow(x2)} distances between each pair of points is returned. (Recall that a \code{face3d} object is replaced by its matrix of \code{vertices}.) If \code{x1} or \code{x2} is a vector corresponding to a single 3d point then a vector of distances is returned.

When \code{type = "min"}, the minimum distance from each set of co-ordinates in \code{x1} to those in \code{x2} is computed.  If \code{index = FALSE} then the vector of minimum distances is returned.  If \code{index = TRUE} then a list with two components is returned.  The minimum distances are in the component \code{distance},  Indices which identify the point in \code{x2} which is closest to each point in \code{x1} are returned in the component \code{indices}.

A further special case occurs when \code{type = "min"}, \code{x1} is a \code{face3d} object and \code{x2curve = TRUE}.  In this case, for each point in \code{x2}, the crossproduct of the normal vector for the closest point in \code{x1} and the vector direction along the curve \code{x2} creates an axis which allows each distance to be signed, indicating that it lies 'above' or 'below' the curve.  The distances returned correspond to the signed distance of each point in \code{x1} to its projection onto the vector direction along the curve.

When \code{type = "minsum"}, the minimum distance from each set of co-ordinates in \code{x1} to those in \code{x2} is computed and these are summed to provide a single overall distance measure.

When \code{x2} is missing, the lengths of the vector, matrix or \code{face3d} object vertices are returned.
}

\seealso{
The function \code{project.face3d} projects points onto a \code{face3d} object, by moving along the directions specified in the corresponding \code{normals} vectors.

The function \code{compare} calculates a variety of distances between two \code{face3d} objects with the same number of vertices.
}
