\name{cutface.face3d}

\alias{cutface.face3d}

\title{Removing the ear effects and adjusting the mandible of face3d object curves.}

\description{
  This function removes the ear effects and adjusts the mandible of the face3d object curves.
}

\usage{
 cutface.face3d(face, graphics = FALSE)
}
\arguments{
\item{face}{A face3d object with \code{landmarks} and \code{curves}}
\item{graphics}{a logical value which controls whether diagnostic graphics are created as the function executes.}
}

\details{
A face3d object with the ear effects and adjusted mandible is returned}

\value{
   An object with the class \code{"face3d"} is returned.
}

\examples{
\dontrun{
   data(face)
   face <- cutface.face3d(face)
}}
