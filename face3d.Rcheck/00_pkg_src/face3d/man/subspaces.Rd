\name{subspaces.face3d}

\alias{subspaces.face3d}

\title{Separate shape space into affine and non-affine subspaces.}

\description{
  Each individual shape is regressed against the mean shape, with the fitted values giving the affine co-ordinates and the residuals plus mean shape providing the non-affine co-ordinates.
}

\usage{
subspaces.face3d(x, weights)
}

\arguments{
\item{x}{either a three-way array giving the co-ordinates of aligned shapes, or a list object containing aligned shapes in the component \code{rotated} and a set of weights in the component \code{weights}.  The latter option makes it convenient simply to pass the output of the \code{gpa.face3d} function.}
\item{weights}{a vector of weights attached to the points which represent each shape.  When \code{x} is a list object, the component \code{weights} will be used, irrespeective of anything passed through this argument.}
}

\details{
The weights created by the \code{gpa.face3d} function are the square roots of the areas of teh model surface which surrounds each point location.  This enables principal components and other analysese to be performed in functional, ratherthan point-based, form.
}

\value{
A list with two or three components:
\item{affine}{the co-ordinates of each shape in affine space;}
\item{nonaffine}{the co-ordinates of each shape in non-affine space;}
\item{rsq}{a two-way array giving the R-squared values which quntify the proportions of variation captured by the affine space, for each shape (row) in each dimension (column).}
}

\examples{
\dontrun{
}}
