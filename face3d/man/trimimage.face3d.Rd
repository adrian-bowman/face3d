\name{trimimage.face3d}

\alias{trimimage.face3d}

\title{Trim an image using landmarks and an underlying model mesh as guides.}

\description{
This function uses supplied landmarks to align the shape horizontally and vertically, and then trims the parts of the image which do not lie close to the model mesh.
}

\usage{
   trimimage.face3d(shape, lmks.vertical = c("sn", "n"),
                    lmks.horizontal = c("exR", "exL"))
}

\arguments{
\item{shape}{a \code{face3d} object to be trimmed.  This must contain both a \code{landmarks} and a \code{mesh} component.}
\item{lmks.vertical}{the names of two landmarks which will be used to rotate the shape to a vertical orientation (first representing bottom and second top).}
\item{lmks.horizontal}{the names of two landmarks which will be used to rotate the shape to a horizontal orientation (first representing right and second left).}
}

\value{
  A \code{face3d} object which is a trimmed version of the input object.
}
