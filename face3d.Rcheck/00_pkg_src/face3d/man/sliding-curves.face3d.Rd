\name{slidingcurves.face3d}

\alias{slidingcurves.face3d}

\title{Slides the curves created using \code{facecurves.face3d} to minimize bending energy between the found curves and the template.}

\description{
  This function slides the curves created using \code{facecurves.face3d} to minimize bending energy between found and subsequently cut curves and the template.
}

\usage{
slidingcurves.face3d(shape1,shape2, n = 100)
}

\arguments{
\item{shape1}{The template.}
\item{shape2}{A \code{face3d} objects containing found curves from \code{facecurves.face3d} that have been cut using \code{cutface.face3d}.}
\item{n}{The multiplicative number of points for resampling.}
}

\details{
A matrix of slid and resampled curve points is returned}

\value{
 A matrix of landmark curve points \code{curves2.optim} that have been slid and resampled.
}

\examples{
\dontrun{
   data(face)
   curves.new  <- slidingcurves.face3d(template, face)
}}
