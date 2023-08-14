\name{make_weights}
\alias{make_weights}

\title{make_weights}

\description{
Create 2 weight matrices to capture asymmetric spatial effects for strip-type intercropping designs.
}

\usage{
make_weights(n)
}

\arguments{
  \item{n}{Number of observations.}
}

\value{
\item{W_BA}{A \eqn{n \times n} spatial weight matrix capturing the locations of type A that are adjacent to locations of type B.}
\item{W_AB}{A \eqn{n \times n} spatial weight matrix capturing the locations of type B that are adjacent to locations of type A.}
}

\references{
1. Hermes, S., van Heerwaarden, J., and Behrouzi, P. (2023). A Spatial Autoregressive Graphical Model with Applications in Intercropping. arXiv preprint, arXiv:2308.04325. \cr
}

\author{
  Sjoerd Hermes, Joost van Heerwaarden and Pariya Behrouzi\cr
  Maintainer: Sjoerd Hermes \email{sjoerd.hermes@wur.nl}
}

\examples{
make_weights(20)
}
