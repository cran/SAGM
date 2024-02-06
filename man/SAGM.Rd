\name{SAGM}
\alias{SAGM}

\title{SAGM}

\description{
This function applies the spatial autoregressive graphical model on a given dataset and array of spatial weight matrices. Different identifiability constraints can be imposed to estimate the \eqn{\Psi_k}. The method allows for both normal and normal-gamma priors, where the values for the hyperparameters can be specified by the user. Returns posterior samples for \eqn{\Theta_{E}} and the \eqn{\Psi_k}.
}
\usage{
SAGM(X, W, prior, constraint, triangular, idx_mat, zeta, kappa, b0, b1,
nBurnin, nIter, verbose)}
\arguments{
  \item{X}{A \eqn{n \times p} matrix, where \eqn{n} is the number of observations and \eqn{p} is the number of variables.}
  \item{W}{A \eqn{n \times n \times 2} array, where the 2 generated spatial matrices are stacked. Note that the order in which the weight matrices are stacked corresponds to the order of the estimated spatial effect matrices.}
  \item{prior}{Prior choice on the spatial effects. Either normal (\code{"normal"}) or normal-gamma (\code{"ng"}).}
  \item{constraint}{Identifiability constraint on the spatial effects. Either symmetric (\code{"symmetric"}), triangular (\code{"triangular"}) or informative (\code{"informative"}).}
  \item{triangular}{Type of triangular restriction. Can be upper-triangular, or lower-triangular, or both, e.g. \code{triangular = c("upper","upper")}. Only has an effect whenever \code{constraint = "triangular"}}
  \item{idx_mat}{A \eqn{\code{nknown} \times 5 \times} matrix, where \code{nknown} is the number of known spatial effects. This matrix contains the indices, means and standard deviations of the known spatial effects that is specified by the user. The matrix only needs to be specified whenever \code{constraint = "informative"} is entered.}
  \item{zeta}{Value of hyperparameter \eqn{\zeta}.}
  \item{kappa}{Value of hyperparameter \eqn{\kappa}.}
  \item{b0}{Value of hyperparameter \eqn{b_0}.}
  \item{b1}{Value of hyperparameter \eqn{b_1}.}
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

# Suppose we have 16 known effects. Here we assign informative normal
# priors to these effects
idx_mat <- matrix(NA, 16, 5)
idx_mat[1,] <- c(1,1,1,1, 0.1)
idx_mat[2,] <- c(1,2,1,1, 0.1)
idx_mat[3,] <- c(1,3,1,1, 0.1)
idx_mat[4,] <- c(1,1,2,1, 0.1)
idx_mat[5,] <- c(1,2,2,1, 0.1)
idx_mat[6,] <- c(1,3,2,1, 0.1)
idx_mat[7,] <- c(4,1,1,-1, 0.1)
idx_mat[8,] <- c(4,2,1,-1, 0.1)
idx_mat[9,] <- c(4,3,1,-1, 0.1)
idx_mat[10,] <- c(4,4,1,-1, 0.1)
idx_mat[11,] <- c(4,1,2,-1, 0.1)
idx_mat[12,] <- c(4,2,2,-1, 0.1)
idx_mat[13,] <- c(4,3,2,-1, 0.1)
idx_mat[14,] <- c(4,4,2,-1, 0.1)
idx_mat[15,] <- c(2,3,1,-1, 0.1)
idx_mat[16,] <- c(2,3,2,-1, 0.1)

W <- array(c(W$W_BA, W$W_AB), dim = c(n,n,2))
est <- SAGM(X = intercrop, W = W, prior = "normal", constraint = "informative",
triangular = c("upper","upper"), idx_mat = idx_mat, zeta = 0.1, kappa = 0.1,
b0 = 0.01, b1 = 0.01, nBurnin = 1000, nIter = 1000, verbose = TRUE)
}
}

