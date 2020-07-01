\name{rp.shapeindex}

\alias{rp.shapeindex}

\title{An interactive demonstration of shape indices}

\description{
  This function produces an interactive demonstration of shape indices.
}

\usage{
rp.shapeindex(panel = TRUE, s = -1)
}

\arguments{
\item{panel}{a logical value which creates an interactive control panel when \code{TRUE} and a static plot when \code{FALSE}.}
\item{s}{the value of shape index used for a static plot.}
}

\details{
The shape index describes the local characteristics of a surface through a quadratic representation.  The two principal curvatures are used to construct this index, which ranges from -1 to 1.  This function shows the surface shapes which correspond to different shape indices.
}

\value{
Nothing is returned
}

\examples{
\dontrun{
rp.shapeindex()
rp.shapeindex(panel = FALSE, s = 0)
}
}
