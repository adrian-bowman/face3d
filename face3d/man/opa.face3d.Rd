\name{opa.face3d}

\alias{opa.face3d}

\title{Register one configuration of 3D points onto another}

\description{
Procrustes registration is used to register one configuration of 3D points onto another.  If the points constitue a discrete representation of a surface, then weighted Procrustes analysis is used to create a functional registration.
}

\usage{
opa.face3d(from, to, carry = from, scale = TRUE, weights, model.mesh = FALSE,
                       return.parameters = FALSE)
}

\arguments{
\item{from}{the point set to be registered.}
\item{to}{the point set onto which \code{from} is to be registered.}
\item{carry}{the point set, or \code{face3d} object, to which the registration will be applied.  This allows the registration to be defined by a subset of points but then applied to a larger set.}
\item{scale}{a logical value which determines whether a scale factor is applied to \code{from} as part of the registration process.}
\item{weights}{a set of weights which will be used for each point when the registration is computed.  When the points constitute a discrete representation of a surface an area weighting is appropriate and previously computed weights can be passed to this argument.}
\item{model.mesh}{a logical value which determines whether the points constitute a discrete representation of a surface.  If the \code{weights} argument is missing and \code{model.mesh} is set to \code{TRUE} then the \code{as.face3d} and \code{area.face3d} functions are used to compute weights.  If the \code{weights} argument is not missing then the \code{model.mesh} argument is ignored.}
\item{return.parameters}{a logical value which controls whether the parameters of the registration are returned.}
\item{exclude}{a character vector giving the names of components of \code{carry}, when it is a \code{face3d} object, which should not be transformed.}
}

\details{
}

\value{
When \code{return.parameters} is set to \code{TRUE} the registered point set is returned.

When \code{return.parameters} is set to \code{FALSE} a list with the following components is returned:
\item{opa}{the registered point set.}
\item{mean.from, mean.to}{the mean values of the \code{from} and \code{to} point sets.}
\item{rotation}{the rotation matrix of the registration.}
\item{scale}{the scale parameter of the registration.}
}

\examples{
\dontrun{
}
}
