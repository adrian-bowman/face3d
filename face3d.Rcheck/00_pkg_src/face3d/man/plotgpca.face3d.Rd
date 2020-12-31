\name{plotgpca.face3d}
\alias{plotgpca.face3d}

\title{Plots of (generalised) principal components for shape data.}

\description{
The output of the \code{gpca.face3d} function is plotted in a variety of ways to show the results of a principal components analysis.
}

\usage{
plotgpca.face3d(gpca, code, high, type = c("scores", "ci"), npc = 10, magnification = 2)
}

\arguments{
\item{gpca}{a list object created by the \code{gpca.face3d} function.}
\item{code}{an optional character variable or factor which splits the data into groups.}
\item{high}{an character value which defines the group which will be assigned the highest possible mean scores.  This is achieved by changing the signs of the scores as required.}
\item{type}{a character variable which controls the kinds of plots which are produced.  \code{"scores"} creates boxplots of the scores associated with each component.  When there are only two unique values in \code{code}, \code{"ci"} draws confidence intervals, using Bonferroni correction for the number of components considered, for the difference between the means of each component.  Again when there are only two unique values in \code{code}, \code{"shapes"} produces an \code{rgl} plot of the two shapes associated with low and high values on the score scale; specifically, these are determined by \code{magnification} standard deviations away from the mean, both positively and negatively. }
\item{npc}{the number of principal components to be considered.}
\item{magnification}{the scaling factor for the standard deviation of the scores when \code{type = "shapes"}.}
}

\value{
  the \code{ggplot} objects, when \code{type} is set to \code{"scores"} and/or \code{"ci"}.
}

\details{
Plotting is performed with the aid of the \pkg{ggplot2} and \pkg{rgl} packages.
}
