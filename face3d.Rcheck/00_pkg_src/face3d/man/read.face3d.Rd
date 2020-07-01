\name{read.face3d}

\alias{read.face3d}

\title{Read obj, pts, dilm, tps, landmarkAscii, jpg, and tif files}

\description{
This function reads seven different types of files which contain information on three-dimensional images. The seven file types are \code{.obj}, which contains the information defining three-dimensional image, \code{.pts}, \code{.dilm}, \code{.tps}, and \code{.landmarkAscii} which contain data on landmarks and possibly also curves, and \code{.jpg} and \code{.tif} which contains the information defining two-dimensional image. In the case of \code{.obj} files, a \code{face3d} object is created, while from \code{.pts}, \code{.dilm}, \code{.tps}, and \code{.landmarkAscii} files a list containing landmarks, and curves if these are present, is created. From \code{.jpg}, RGB imagematrix object out of 3 rasters (red, green, and blue) or one raster (grey) is created. From \code{.tif}, RGB pixmap object out of 3 rasters (red, green, and blue) or one raster (grey) is created.
}

\usage{
read.face3d(filename, colour = TRUE, jpgfile.addition = "-small",
            convert.gray = FALSE, quality = 100)
}

\arguments{
\item{filename}{the name of the .obj, .pts, .dilm, .tps, .landmarkAscii, .jpg, or .tif file to be read.}
\item{colour}{a logical value which determines whether colour information for each point in the three-dimensional image is read from an associated \code{.jpg} file, if this is present.}
\item{jpgfile.addition}{the addition to the part of \code{objfile} before the \code{.obj} file extension, to form the name of the \code{.jpg} file containing the photographic image from which colour information can be extracted.}
\item{convert.gray}{a logical value which determines whether RGB or gray raster is created.}
\item{quality}{optional scaling factor to change the two-dimensional image size, should be a value between 0 and 100 (in percentages).}

}

\details{
The type of information to be read is identified by the extension of \code{filename}, either \code{.obj}, \code{.pts}, \code{.dilm}, \code{.tps}, \code{.landmarkAscii}, \code{.jpg} or \code{.tif}.

In the case of \code{.obj} files, information on the co-ordinates and indices is extracted.  If the \pkg{jpeg} package is available then the colours of the co-ordinates are also extracted from the file with the name \code{objfile} where the \code{.obj} extension is replaced by \code{paste(jpgfile.addition, "jpg", sep = ".")}.

In the case of \code{.pts} files, information on landmarks, and if present curves, is extracted.
}

\value{
In the case of \code{.obj} files, the function returns a \code{face3d} object in the form of a list with two or three components.  The component \code{coords} is a matrix with three columns containing the three-dimensional co-ordinates of the points.  The component \code{triples} is a vector of integers containg the indices of \code{coords} which should be passed to the \pkg{rgl} function \code{triangles} to render the surface.  If the \pkg{ReadImages} package is available, the third component \code{colour} contains the colours of the points corresponding to the rows of \code{coords}.

In the case of \code{.pts} files, the function returns a list with a component \code{lmks} containing landmarks and, if curves are present, a second component \code{curves} containing curves as sets of locations along each curve length.

In the case of \code{.dilm}, \code{.tps}, and \code{.landmarkAscii} files, the function returns a list with a component \code{lmks} containing landmarks.

}
