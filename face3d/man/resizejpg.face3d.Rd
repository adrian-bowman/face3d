\name{resizejpg.face3d}

\alias{resizejpg.face3d}

\title{Resize the jpg file associated with a camera capture}

\description{
This function uses the \code{convert} facility from the ImageMagick suite of tools to resize a \code{.jpg} image.  This is usually necessary in order to read the jpg file as part of the creation of an \code{R} face3d object which has colour associated with each three-dimensional surface point.
}

\usage{
   resizejpg.face3d(filename, recursive = TRUE, xsize = 1024, ysize = 1024,
                     file.extension = ".jpg", file.addition = "-small",
                     overwrite = FALSE, verbose = 1)
}

\arguments{
\item{filename}{the name of the .jpg file to be reduced, or a vector of filenames, or a folder name.  If the filname does not end in \code{file.extension} then it is assumed to be a folder name within which \code{list.files} will be used to identify all files in that folder.  If \code{filename} is missing, files are sought in the current directory.  File searching is performed recursively if \code{recursive == TRUE}.}
\item{recursive}{a logical value which determines whether file searching is performed recursively, when \code{filename} is missing.}
\item{xsize}{an integer which determines the horizontal number of pixels in the reduced image.}
\item{ysize}{an integer which determines the vertical number of pixels in the reduced image.}
\item{file.extension}{the text which defines the file extension of the files sought when \code{filename} is missing.}
\item{file.addition}{the text to be added before the file extension, to provide the name of the reduced file.}
\item{overwrite}{a logical value which determines whether existing files with names which end in \code{file.addition} followed by \code{file.extension} are over-written.}
\item{verbose}{an integer which determines the level of messages produced by the function.  The function will operate without any messages if \code{verbose = 0}.  When \code{verbose} is set to \code{1} the names of the newly created files will be printed.  When \code{verbose} is set to \code{2} a message will be printed for all files, including those which already exist and are not recreated.}
}

\value{
  Nothing is returned.
}
