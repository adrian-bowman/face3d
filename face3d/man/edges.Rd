\name{edges}

\alias{edges}

\title{Locating the edges of a shape}

\description{
This function located the points which lie on the edges of a shape.
}

\usage{
edges(shape)
}

\arguments{
\item{shape}{a \code{face3d} object.}
}

\details{
The function identifies the pairs of points which define an edge which appears only once in the triangulation of the shape.  Edges which lie in the interior part of the surface are present in adjacent triangles and so appear more than once.

If a \code{lines} component is present in \code{shape}, this is ignored.

Note that single edges can sometimes also appear in areas of the surface where the trianglution is confused, due to imperfections in the captured image.  For example, strands of hair can sometimes cause this kind of problem in facial images.  There can also be confusion at the perimeter of the surface.  In these cases, multiple edges are produced.
}

\value{
A list where each element is a vector of integers which identify the elements of \code{face$coords} which form a particular edge. The identified points form a sequence of adjacent locations, with the last index repeating the starting point, so that the edges can be drawn by the \code{rgl::lines3d} function, as in the example below.

The elements of the returned list is ordered by the number of points in each edge.
}

\examples{
\dontrun{
edgs <- edges(template.male)
sapply(edgs, length)
plot(template.male)
lapply(edgs, function(x) rgl::lines3d(template_male$coords[x, ], lwd = 2, col = "blue"))
}}
