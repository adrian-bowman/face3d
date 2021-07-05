\name{summary.face3d}

\alias{summary.face3d}

\title{Provide a simple summary of a face3d object}

\description{
  This function provides a simple summary of a \code{face3d} object.  The \code{class} of the object is checked to be \code{face3d}.  Checks for duplicated co-ordinates and degenerate triangles can also be performed.  Printing can be suppressed by adding the argument \code{print = FALSE}.
}

\usage{
\method{summary}{face3d}(shape, checks = FALSE, fix = FALSE, ...)
}

\arguments{
\item{shape}{the face3d object.}
\item{checks}{a logical value which controls whether various checks are applied to the face3d object.}
\item{fix}{a logical value which controls whether various checks are applied to the face3d object.}
\item{...}{other arguments which may be passed.  Currently, the only argument used is \code{print} which is interpreted as a logical value which determines whether a summary table and messages about isolated, duplicated and collinear points are printed.  Messages about inappropriate values are always printed.}
}

\details{
If \code{fix = TRUE} then the argument \code{checks} is also forced to be \code{TRUE} so that any problems needing to be fixed are identified. 
}

\value{
When \code{fix = FALSE}, the returned object is a list with components:
\item{n.vertices}{the number of vertices.}
\item{n.faces}{the number of faces.}
\item{ranges}{the ranges of x, y, and z.}

If, in addition, \code{checks == TRUE}, additional components are:
\item{duplicated}{the indices of any duplicated vertices.}
\item{collinear}{a three-column matrix containing the indices of the points in any collinear triangles.}
\item{isolated}{the indices of any isolated points which are not part of teh triangulation.}

If \code{fix = TRUE} then the checks described above are carried out, any duplicated \code{vertices} and collinear \code{triangles} are removed and the returned object is the adjusted \code{shape}.
}

\examples{
\dontrun{
   data(face)
   summary.face3d(face)
   summary.face3d(face, checks = TRUE)
   summary.face3d(face, fix = TRUE)
}}
