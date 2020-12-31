\name{summary.face3d}

\alias{summary.face3d}

\title{Provide a simple summary of a face3d object}

\description{
  This function provides a simple summary of a face3d object.  The \code{class} of the object is checked to be \code{face3d}.  Checks for duplicated co-ordinates and degenerate triangles can also be performed.  Printing can be suppressed by adding the argument \code{print = FALSE}.
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
When \code{fix = FALSE}, the returned object as a list with components for the number of vertices \code{n.vertices}, the number of faces \code{n.faces} and the ranges of x, y, and z \code{ranges}. If, in addition, \code{checks == TRUE}, checks for duplicated co-ordinates (\code{object$coords}) and collinear triangles (indices in \code{object$triples}) are carried out.  If any are found, the indices of any duplicated co-ordinates are returned in component \code{duplicated} and the indices of the points in any collinear triangles are returned in component \code{collinear}.

If \code{fix = TRUE} then the checks described above are carried out, any duplicated \code{coords} and collinear triangles are removed and the returned object is the adjusted \code{shape}.
}

\examples{
\dontrun{
   data(face)
   summary.face3d(face)
   summary.face3d(face, checks = TRUE)
   summary.face3d(face, fix = TRUE)
}}
