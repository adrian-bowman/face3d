\name{plot.face3d}

\alias{plot.face3d}
\alias{display.face3d}

\title{Display the shape in an rgl window}

\description{
  This function display the shape in an \pkg{rgl} window, with some simple controls on the nature of the display.
}

\usage{
plot.face3d(shape, display = "surface", colour = "texture", palette,
               breaks, scaling = "quantile", range.colour, key = FALSE, cex.axis = 1,
               new = FALSE, add = FALSE, subset, vector.scale = 1, windowRect = c(300, 300, 800, 800), ...)
\method{plot}{face3d}(x, ...)
}

\arguments{
\item{shape, x}{the face3d object.}
\item{display}{a character vector which determines whether the shape is displayed as \code{"points"}, \code{"spheres"}, \code{"mesh"} or \code{"surface"}.  With \code{spheres}, the additional argument \code{radius} can be passed to control size.  The setting \code{"normals"} adds the normals to the current plot while the settings \code{"directions 1"} and \code{"directions 2"} add the principal directions.  The argument \code{vector.scale} can be used to adjust the length of each vector for a suitable visual display.  Note that \code{display} is a vector, so that multiple features can be created in the plot with a single function call.}
\item{colour}{If this argument is set to "texture" (the default) then the natural texture of the shape (in \code{shape$colour}) will be used, if available.  If a single colour is supplied, this will be replicated across all the points in \code{shape} and the display will use the \pkg{rgl} control \code{lit = TRUE}.  
	
If a vector, whose length matches the number of points or number of triangles (\code{shape$triples)/3}) in \code{shape}, is suppled then colours will be created using the \code{cut} function.  If \code{colour.range} is supplied then this defines the range of colours for plotting, otherwise the range of \code{colour} will be used.

If \code{colour} is set to \code{"shape index"} then the shape indices in the matrix \code{shape$shape.index}, if present, are used to define the colour.  If \code{colour} is set to \code{"normal-x"}, \code{"normal-y"} or \code{"normal-z"}, then the absolute values of the coordinates of the normals in the particular direction (or its difference) from the matrix \code{shape$normals}, if present, are used to define the colour.}

\item{palette}{This defines the colour palette which is used when the \code{colour} argument is set to a numeric vector.  If \code{palette} is set to \code{"diverge"} then the \code{diverge_hcl} palette from the \code{colorspace} package will be used, if this is available, otherwise topo-like colours are used.  This latter palette can be requested directly by setting \code{palette = "diverge_topo"}.  If \code{col.palette} is set to a vector of colours then this colour palette will be used.  If \code{palette} is not specified and \code{colour} is a numeric variable whose range spans \code{0} then \code{palette} will be set to \code{"diverge"}.}
\item{breaks}{a vector of numerical values which will be used to define the allocation of colours when a numeric variable is passed to the \code{colour} argument.  If \code{breaks} is set to a single value then this will be taken as the number of break points which are to be computed from the range of the colour variable.  In this latter case the default number of breaks is \code{50}.}
\item{scaling}{a character value which determines whether colour breaks are equally spaced a standard scale (\code{"linear"}), equally spaced on the interior section but with final categories (indicated by white hatching) which accommodate outliers (\code{"trimmed"}), or determined by quantiles (\code{"quantile"}).  If the \code{range.colour} argument is specified then \code{range.colour} is set to \code{"linear"}, otherwise the default is \code{"trimmed"}.  The details of these options are described below.}
\item{range.colour}{a numeric vector of length 2 which defines the end points of the scale when \code{colour} is set to a numeric variable.  The default is to use the range of the colour variable.  Setting \code{range.colour} to the range of a collection of colour variables can be a useful way of retaining a common colour key across a collection of separate plots.}
\item{key}{a logical value which determines whether a colour key is displayed in the R plot window.}
\item{cex.axis}{when \code{key} is \code{TRUE}, this numeric value allows control over the size of the printed axis labels.}
\item{new}{a logical value which determines whether a new \pkg{rgl} window is created or the current window is cleared before the new object is displayed.}
\item{add}{a logical value which determines whether the object is added to the current scene or the current scene is cleared and a new one is created.}
\item{subset}{a logical vector, whose length matches \code{nrow(shape$coords)}, which is used in a call to \code{subset} to subset the shape before it is displayed.}
\item{vector.scale}{a numeric value which is used to scale the normal vectors hen these are displayed through the setting \code{display = "normals"}.  This allows the display to be adjusted to fit the scale of the co-ordinates of the image.}
\item{windowRect}{an integer vector of length 4 which is passed to the \pkg{rgl} function \code{open3d} to controls the c(left, top, right, bottom) pixel co-ordinates of the \pkg{rgl} window.}
\item{theta, phi, zoom}{numerical values which are passed to the corresponding arguments of \code{view3d} function.  The default values are 0, 0, and 0.7 respectively.}
\item{...}{further arguments which are passed to \pkg{rgl} functions.  This allows detailed control of \code{points3d} and \code{triangles3d}, when \code{display} is set to \code{"points"} or \code{"surface"} respectively.}
}

\details{
If \code{scaling} is set to \code{"trimmed"} then the \code{boxplot} function is used to identify the extremes of the lower and upper whiskers.  The scale is then set to linear within this range, with any values of the \code{colour} variable which are more extreme allocated to the lower or upper end of the colour range.
}

\value{
  Nothing is returned.
}

\examples{
\dontrun{
   plot(template.male)
}}
