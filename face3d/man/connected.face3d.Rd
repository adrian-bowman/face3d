\name{connected.face3d}

\alias{connected.face3d}

\title{Identify the connected parts of a shape object}

\description{
  This function identifies the connected parts of the \code{face3d} object \code{shape}.
}

\usage{
connected.face3d(shape, order.by.area = FALSE)
}

\arguments{
\item{shape}{a \code{face3d} object.}
\item{order.by.area}{a logical value which determines whether the different parts of the object will be ordered by surface area (largest to smallest).  With large objects surface area can be slow to compute (with \code{area.face3d}) so the default is set ot \code{FALSE}.}
}

\details{
The information in \code{shape$triangles} is examined to identify the connected parts of the object.  Within each connected part it is possible to move between any two vertices along the edges of the triangulation defined by \code{shape$triangles}.
}

\value{
An integer vector of indicators, of length \code{nrow(shape$coords)}, is returned to identify the connected parts of \code{shape}.  Values of \code{0} indicate isolated points.  The connected parts are in decreasing order of the number of vertices or, if \code{order.by.area} is set to \code{TRUE}, surface area.
}

\examples{
\dontrun{
   data(face)
   connected.face3d(face)
}}
