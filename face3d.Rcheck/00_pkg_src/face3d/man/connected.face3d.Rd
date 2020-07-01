\name{connected.face3d}

\alias{connected.face3d}

\title{Identify the connected parts of a shape object}

\description{
  This function identifies the connected parts of the \code{face3d} object \code{shape}.
}

\usage{
connected.face3d(shape)
}

\arguments{
\item{shape}{a \code{face3d} object.}
}

\details{
The information in \code{shape$triples} is examined to identify the connected parts of the object.  Within each connected part it is possible to move between any two vertices along the edges of the triangulation defined by \code{shape$triples}.
}

\value{
An integer vector of indicators, of length \code{nrow(shape$coords)}, is returned to identify the connected parts of \code{shape}.  Values of \code{0} indicate isolated points.
}

\examples{
\dontrun{
   data(face)
   connected.face3d(face)
}}
