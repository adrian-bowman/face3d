\name{closestcontrol.face3d}

\alias{closestcontrol.face3d}

\title{Find the shape which is as close as possible to the cases while remaining consistent with a control set}

\description{
  This function finds the shape which is as close as possible to the cases while remaining consistent with a control set.
}

\usage{
closestcontrol.face3d(cases, gpa, npc, pc.proportion = 0.8,
                      display = c("none", "unusual")
}

\arguments{
\item{cases}{a three-dimensional array giving the three-dimensional co-ordinates (array dimension 2) of surface locations (array dimension 1) for a set (array dimension 3) of cases.}
\item{gpa}{a three-dimensional array, in the same format as the \code{cases}, with the co-oreinates of a set of cases.  Note that these are assumed to have been Procrustes registered.}
\item{npc}{a positive integer controlling the number of principal components used in the calculations.  If this argument is missing, \code{npc} will be set to the smallest value for which the cumulative proportion of variance explained in the controls is at least \code{pc.proportion}.}
\item{pc.proportion}{the threshold for the cumulative proportion of variance explained in the controls, used to set the value of \code{npc} if this is not specified explicitly in the call to the function.}
\item{display}{a character value or vector which controls the plots produced.  If \code{"none"} is an element of \code{display} then no plots will be displayed, although \pkg{ggplot2} plots will be returned by the function.  If \code{"unusual"} is an element of \code{display} then only unusual cases will be highlighted in the plots, while if \code{"all"} is an element then all cases will be shown.  The setting \code{"unusual"} is activated when none or both of \code{"all"} and \code{"unusual"} is present.}
}

\details{
Something.
}

\value{
A list with the following components: 
\item{closest.control}{.}
\item{plot.pc}{.}
\item{plot.residuals}{.}
\item{others}{.}
}
