\name{sample.face3d}

\alias{sample.face3d}
\alias{sample.face3d}

\title{Create a sample of vertices with minimum spacing from a shape or set of vertices}

\description{
Create a sample of vertices with minimum spacing from a shape or set of vertices
}

\usage{
sample.face3d(x, spacing)
}

\arguments{
\item{x}{a \code{face3d} objecy or a three-column matrix containg the co-ordinates of a set of vertices.}
\item{spacing}{the maximum distance from any vertex of \code{x} to the closest sampled vertex.}
}

\details{
}

\value{
A vector of indices for the sampled vertices.
}

\examples{
\dontrun{
siclr <- sicolour.face3d(template_female$shape.index)
}}
