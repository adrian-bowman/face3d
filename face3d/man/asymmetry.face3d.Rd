\name{asymmetry.face3d}

\alias{asymmetry.face3d}

\title{Construct asymmetry scores for shapes}

\description{
  This function constructs asymmetry scores for shapes.  Sub-scores for regions can also be computed.}

\usage{
asymmetry.face3d <- function(X, sides = c("left", "right"), reference.size,
                             region, match.subset, controls)
}

\arguments{
\item{X}{a matrix or three-way array whose rows provide the three-dimensional co-ordinates of the vertices of landmarks, curves, or surface meshes.  In the ase of an array, the third dimension indexes a set of different shapes.  The rownames of X are used to identify the two matching parts of the shape to which each vertex belongs.}
\item{sides}{a character vector of length 2 which gives the identifiers of the two parts of the shape to which symmetry applies.  The points are paired by matching the rownames after the removal of the \code{sides} string.}
\item{reference.size}{a scaling factor which allows the asymmetry score to be placed on a scale which is comparable with a reference set.}
\item{region}{a list whose components give the indices of the vertices associated with locql regions of the shape.  Asymmetry scores will then be computed for eaqch region.  The componenets of the list should be named.  A \code{"global"} component will be created automatically, to ensure thqt the global asymmetry is computed and to all a symmetric shape to be constructed.}
\item{match.subset}{a set of indices which defines the vertices to be used in Procrustes matching when the asymmetry score is computed.  This applies only to the global case.  The default uses all available vertices.}
\item{controls}{a matrix containing asymmetry scores (columns) for a reference set of individuals (rows), to be used to construct a plot which compares the scores on the new shapes with a control set.  The object created by a call to \code{asymmetry.face3d} on this control set can also been passed in this argument for convenience.  The \code{asymmetry} scores and the \code{reference.size} will be picked up from this list, as will \code{region} if it is not specified in the \code{region} argtument of the call to the function.}
}

\details{
A global asymmetry score is calculated as well as scores for sub-regions when these are specified.  When \code{X} is a matrix specifying an individual shape, a symmetrised version of \code{X} is also returned.}

\value{
A list with components:
\item{asymmetry}{a matrix providing the asymmetry scores for each shape; these scores measure the distance from the observed shape to its symmetrised version;}
\item{region}{a list whose comnponents give the vertices used to define each region of the shape for which an asymmetry score has been computed;}
\item{size}{a matrix giving the size of each global and regional component, for each shape;}
\item{reference.size}{the scaling factor used to compare the asymmetry scores to those of a reference set; when no reference set or \code{reference.size} has been specified in the call to the function, the mean sizes across the shapes in \code{X} will be returned;}
\item{symmetric}{a symmetrised version of \code{X}, when this contains a single shape;}
\item{plot}{a plot of the asymmetry scores compared to a reference set, when this set has been supplied; the plot is in \pkg{ggplot2} form which allows subsequent amendment of details.}
}

\examples{
\dontrun{
   asym <- asymmetry.face3d(template.male$mesh)
}}
