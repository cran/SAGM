\name{SAGM}
\alias{SAGM}

\title{SAGM}

\description{
This function applies the spatial autoregressive graphical model on a given dataset and array of spatial weight matrices. Different identifiability constraints can be imposed to estimate the \eqn{\Psi_k}. The method allows for both normal and normal-gamma priors, where the values for the hyperparameters can be specified by the user. Returns posterior samples for \eqn{\Theta_{E}} and the \eqn{\Psi_k}.
}
\usage{
SAGM(X, W, prior, constraint, triangular, zeta, kappa, c0, c1, nBurnin, nIter, verbose)}
\arguments{
  \item{X}{A \eqn{n \times p} matrix, where \eqn{n} is the number of observations and \eqn{p} is the number of variables.}
  \item{W}{A \eqn{n \times n \times 2} array, where the 2 generated spatial matrices are stacked. Note that the order in which the weight matrices are stacked corresponds to the order of the estimated spatial effect matrices.}
  \item{prior}{Prior choice on the spatial effects. Either normal (\code{"normal"}) or normal-gamma (\code{"ng"}).}
  \item{constraint}{Identifiability constraint on the spatial effects. Either symmetric (\code{"symmetric"}) or triangular (\code{"triangular"}).}
  \item{triangular}{Type of triangular restriction. Can be upper-triangular, or lower-triangular, or both, e.g. \code{triangular = c("upper","upper")}. Only has an effect whenever \code{constraint = "triangular"}}
  \item{zeta}{Value of hyperparameter \eqn{zeta}.}
  \item{kappa}{Value of hyperparameter \eqn{kappa}.}
  \item{c0}{Value of hyperparameter \eqn{c_0}.}
  \item{c1}{Value of hyperparameter \eqn{c_1}.}
  \item{nBurnin}{Number of burnin samples.}
  \item{nIter}{Number of post-burnin Gibbs samples.}
  \item{verbose}{Return progress of parameter estimation (\code{True}) or not (\code{False}).}

}

\value{
\item{Theta}{A \eqn{p \times p \times \code{nIter}} array consisting of the post-burnin samples for the within-plot dependencies.}
\item{Psi}{A \eqn{p \times p \times 2 \times \code{nIter}} array consisting of the post-burnin samples for the between-plot effects. The order of the third dimension of the array corresponds to that of the \code{W}.}
\item{lambda_sq}{A \eqn{p \times p \times \code{nIter}} array consisting of the post-burnin samples for \eqn{\Lambda^2}. This output is of secondary interest.}
\item{tau_sq}{A vector of length \code{nIter} consisting of the post-burnin samples for \eqn{\tau^2}. This output is of secondary interest.}
\item{accpt_rate}{Value of the acceptance rate of the Metropolis Hastings step.}
}

\references{
1. Hermes, S., van Heerwaarden, J., and Behrouzi, P. (2023). A Spatial Autoregressive Graphical Model with Applications in Intercropping. arXiv preprint, arXiv:2308.04325. \cr
}

\author{
  Sjoerd Hermes, Joost van Heerwaarden and Pariya Behrouzi\cr
  Maintainer: Sjoerd Hermes \email{sjoerd.hermes@wur.nl}
}

\examples{
\donttest{
data(intercrop)
n <- nrow(intercrop)
W <- make_weights(n)
W <- array(c(W$W_BA, W$W_AB), dim = c(n,n,2))
est <- SAGM(X = intercrop, W = W, prior = "normal", constraint = "symmetric",
triangular = c("upper","upper"), zeta = 0.1, kappa = 0.1, c0 = 0.01, c1 = 0.01,
nBurnin = 1000, nIter = 1000, verbose = TRUE)
}
}

