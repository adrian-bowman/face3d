\name{face3d-package}
\alias{face3d-package}
\alias{face3d}

\docType{package}

\title{Tools for the analysis of three-dimensional surface images}

\description{
The \pkg{face3d} package provides tools to analyse data in the form of three-dimensional facial surfaces, such as those arising from laser scanning or stereo-photogrammetry.  Some functions can handle surfaces of any type while other functions are tailored to human faces.  Some tools are for organising and manipulating the data through operations such as reading, structuring, subsetting and displaying.  Other tools allow the identification of specific features, particularly for facial images.  Tools for statistical modelling are also provided.
}

\details{
\tabular{ll}{
Package: \tab face3d\cr
Type: \tab Package\cr
Version: \tab 0.1-1\cr
Date: \tab 2012-12-06\cr
License: \tab GNU\cr
}

\code{\link{asymmetry.face3d}}: compute the asymmetry of a face3d object\cr
\code{\link{closest.face3d}}: find the closest vertex of a face3d object to a point of interest\cr
\code{\link{connected.face3d}}: identify the connected parts of a face3d object\cr
\code{\link{distances.face3d}}: compute distnaces of various types between shapes\cr
\code{\link{facecurves.face3d}}: identify facial curves on a face3d object\cr
\code{\link{plot.face3d}}: display of face3d objects\cr
\code{\link{procrustes.face3d}}: generalised Procrustes registration\cr
\code{\link{index.face3d}}: compute shape and other indices\cr
\code{\link{interpolate.face3d}}: interpolation of a scalar quantity across a face3d object\cr
\code{\link{orient.face3d}}: orient a face into a frontal view\cr
\code{\link{planepath.face3d}}: identify a planar facial curve on a face3d object\cr
\code{\link{read.face3d}}: create face3d and other objects from a file\cr
\code{\link{resample.face3d}}: create face3d and other objects from a file\cr
\code{\link{resizejpg.face3d}}: resize jpg files for colour reading\cr
\code{\link{rotate.face3d}}: rotate a face3d object into a standardised position\cr
\code{\link{rp.shapeindex}}: view an interactive demonstration of shape indices\cr
\code{\link{subset.face3d}}: subset a face3d object\cr
\code{\link{summary.face3d}}: give a summary of a face3d object\cr
\code{\link{warp.face3d}}: warp one face3d object onto another\cr
\code{\link{write.face3d}}: write a face3d and other objects to a file\cr
}

\author{
Adrian Bowman, Stanislav Katina and Liberty Vittert.

Maintainer: Adrian Bowman <adrian.bowman@glasgow.ac.uk>
}

\keyword{package}
