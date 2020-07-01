\name{crossproduct}

\alias{crossproduct}

\title{Construct the crossproduct of pairs of vectors}

\description{
  This function constructs the crossproduct of pairs of vectors.
}

\usage{
crossproduct(a, b, scale = TRUE)
}

\arguments{
\item{a, b}{vectors of length 3 or matrices with 3 columns whose crossproducts are to be computed.}
\item{scale}{a logical variable which determines whether the computed crossproduct(s) have unit length.}
}

\details{
When matrices are specified in the arguments, crossproducts are computed for each row.
}

\value{
a vector of the computed crossproducts.}

% \examples{
% \dontrun{
% }}
