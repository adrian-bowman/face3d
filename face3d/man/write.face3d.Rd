\name{write.face3d}

\alias{write.face3d}

\title{Write obj, ply, dmp, jpg, and tif files}

\description{
This function writes five different types of files which contain information on three-dimensional images. The five file types are \code{.obj}, with simple information on vertices and faces, \code{.ply}, which contains the information defining three-dimensional image with colour scaled 0 to 1, \code{.jpg} and \code{.tif} which contains the information defining two-dimensional image. In the case of \code{.dmp} files, a \code{face3d} object is exported, while from \code{.jpg}, RGB imagematrix object out of 3 rasters (red, green, and blue) or one raster (grey) is exported. From \code{.tif}, RGB pixmap object out of 3 rasters (red, green, and blue) or one raster (grey) is exported.  
}

\usage{
write.face3d(filename, shape, jpgfile, tiffile, quality = 100) 
}

\arguments{
\item{filename}{the name of the .ply, .dpm, .jpg, or .tif file to be write.}
\item{shape}{a \code{face3d} object in the form of a list with two or three components created by \code{read.face3d} function.}
\item{jpgfile}{RGB imagematrix object created by \code{read.face3d} function.}
\item{tiffile}{RGB pixmap object created by \code{read.face3d} function.}
\item{quality}{optional scaling factor to change the two-dimensional image size, should be a value between 0 and 100 (in percentages).}
}

\details{
The type of information to be write is identified by the extension of \code{filename}, either \code{.ply}, \code{.dmp}, \code{.jpg} or \code{.tif}.

In the case of \code{.obj} a very simplified form of file is exported, with information only on vertices and faces.  This can be helpful as a means of creating representations which only have triangular faces.  The \code{read.face3d} function will create a \code{face3d} object in R which is based only on triangles, even where the input file has more complex faces.  The \code{write.face3d} function can then be used to create a \code{.obj} representation if the simplified triangular form is required.

In the case of \code{.ply} and \code{.dmp} files, information on the co-ordinates, indices, and colour is exported.

In the case of \code{.jpg} and \code{.tif} files, the figure is ploted in two-dimensional R-window and then exported.
}

\value{
In the case of \code{.ply} file, \code{.ply} file is exported.

In the case of \code{.dmp} file, \code{.dmp} file is exported.

In the case of \code{.jpg} file, \code{.jpg} figure is exported.

In the case of \code{.tif} file, \code{.tif} figure is exported.

}
