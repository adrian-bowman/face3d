\name{sicolour.face3d}

\alias{sicolour.face3d}
\alias{sicolor.face3d}

\title{Return the colour associated with values of shape index}

\description{
There is a standard colour palette associated with values of shape index.  This function returns the colours associated with the values of the shape indices in \code{x}.
}

\usage{
sicolour.face3d(x)
}

\arguments{
\item{x}{a vector of numeric values.}
}

\details{
The shape index lies between \code{-1} and \code{1} so a check is carried out on the elements of \code{x}.  If any values outside this range are identified then an error is reported.
}

\value{
A vector of colour values.
}

\examples{
\dontrun{
siclr <- sicolour.face3d(template_female$shape.index)
}}
