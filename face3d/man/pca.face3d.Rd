\name{pca.face3d}
\alias{pca.face3d}

\title{Interactive plots of three-dimensional shape data, including principal components.}

\description{
This function computes a PCA of aligned shape data.
}

\usage{
pca.face3d(x, group, weights = rep(1, nrow(psc)), alpha = 0)
}

\arguments{
\item{x}{either a matrix of aligned shapes or a list object, created by the \code{gpa.face3d} function. In the latter case the weights are set to \code{x$weights}.}
\item{group}{a factor with identifies each individual shape to one of two groups.}
\item{weights}{a vector of weights used to perform functional PCA when the objects are model meshes.  The weights are the square roots of surface areas associated with each point on the model mesh.  When \code{x} is the output of the \code{gpa.face3d} function, \code{weights} is set to \code{x$weights}.}
\item{alpha}{an index which allows 'generalised PCA' to be produced.}
}

\value{
   a list object which contains information about the PCA.
}

