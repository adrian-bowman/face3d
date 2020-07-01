\name{procrustes.face3d}

\alias{procrustes.face3d}

\title{Generalized Procrustes registration of (semi)landmarks}

\description{
This function takes a three-dimensional array of point configuations and registers these by matching the configurations together as closely as possible.  Generalised Procrustes registration is used.
}

\usage{
procrustes.face3d(lmks, scale = TRUE, threshold = 0.03, rotation = "coronal",
                  rotation.lmks = c("sn", "n", "exR", "exL"), subset)}

\arguments{
\item{lmks}{a three-dimensional array of (semi)landmark coordinates.  The first dimension indexes the points in each configuration, the second dimension refers to the three spatial co-ordinates and the third dimension indexes the different configurations.}
\item{scale}{a logical value which controls whether registration includes the scaling of the configurations.}
\item{threshold}{a value which controls the stopping rule for the iterations.}
\item{rotation}{a character value which specifies a type of rotation for the mean of the registeted configurations.  Valid values are "coronal", "sagittal sinister", "sagittal dexter" and "transversal", as used in \code{rotate.face3d}.}
\item{rotation.lmks}{the indices of the landmarks used for rotation of the mean configuration.  See \code{rotate.face3d}.}
\item{subset}{a logical vector or a set of indices which defines a subset of the \code{lmks} which will be used for registration..}
}

\details{
An array of registered coordinates.
}

\value{
The function returns a list with the following components.
\item{rotated}{an array of registered coordinates;}
\item{mean}{a mean configuration;}
\item{centroid.sizes}{the centroid sizes of the configurations.}
}
