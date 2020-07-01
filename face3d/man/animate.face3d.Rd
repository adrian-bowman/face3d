\name{animate2.face3d}

\alias{animate2.face3d}

\title{Create a sequence of images to display the interpolated path between two shapes.}

\description{
This function creates a sequence of images that interpolate between two shapes.  The images can be displayed as an animation in an \code{rgl} window and/or written to a folder as the raw material for the construction of a movie.
}

\usage{
animate2.face3d(shapes, ngrid = 4, display = TRUE, new = FALSE, image.folder, ...)
}

\arguments{
\item{shapes}{a list of \code{face3d} objects to be interpolated in successive pairs.}
\item{ngrid}{the number of interpolated steps between each successive pair of \code{face3d} objects.  This has the default value 4 for viewing on-screen.  Larger values will provide smoother animations which will be helpful in the creation of a movie from stored files.}
\item{display}{a logical value which controls whether images are displayed in the \code{rgl} window.}
\item{new}{a logical value which controls whether a new \code{rgl} windown is created or the current one is used.} 
\item{image.folder}{the name of a folder in the present working directory which will be used to store .png files with the images of the interpolated steps.  If this argument is not supplied then no images will be stored.}
\item{...}{other arguments which will be passed on to the \code{plot.face3d} function.}
}

\details{
Note that the \code{image.folder} should be present in the current working directory.  If this folder does not exist the function will stop and report that.

Files stored in the \code{image.folder} have names of the form \code{temp-i.png} where \code{i} indexes the sequence.

Stored images may be converted into a movie through the use of tools such as \code{ImageMagick}.
}

\value{
  None.
}

\examples{
\dontrun{
   data(face)
}}
