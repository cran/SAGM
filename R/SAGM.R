
SAGM <- function(X, W, prior = "ng", constraint = "triangular", triangular = c("upper","lower"), zeta = 0.1, kappa = 0.1, c0 = 0.01, c1 = 0.01,
                  nBurnin = 1000, nIter = 1000, verbose = TRUE)
{
  p <- ncol(X)
  accpt_rate <- 0
  S_obs <- nrow(X)
  K <- 0.5*dim(W)[3]
  # set initial values
  Theta <- diag(p)

  if(prior == "ng"){ # normal gamma prior (shrinkage)
    alpha <- array(rgamma(2*K*p^2, shape = kappa, rate = kappa),
                   dim = c(p, p, 2*K))
    omega <- rgamma(1, shape = c0, rate = c1)
    omega <- omega^2

    Results <- list(Theta = array(dim = c(p, p, nIter)),
                    Psi = array(dim = c(p, p, 2*K, nIter)),
                    alpha = array(dim = c(p, p, 2*K, nIter)),
                    lambda_sq = array(dim = c(p, p, nIter)),
                    tau_sq = numeric(nIter),
                    accpt_rate <- numeric(1))
  }

  else if(prior == "normal"){ # uninformative normal prior
    Results <- list(Theta = array(dim = c(p, p, nIter)),
                    Psi = array(dim = c(p, p, 2*K, nIter)),
                    lambda_sq = array(dim = c(p, p, nIter)),
                    tau_sq = numeric(nIter),
                    accpt_rate <- numeric(1))
  } else{ # uniform prior
    Results <- list(Theta = array(dim = c(p, p, nIter)),
                    Psi = array(dim = c(p, p, 2*K, nIter)),
                    lambda_sq = array(dim = c(p, p, nIter)),
                    tau_sq = numeric(nIter),
                    accpt_rate <- numeric(1))
  }

  # cannot set Psi to zero because it causes computation issues
  Psi <- array(0.00001, dim = c(p, p, 2*K))

  lambda_sq <- matrix(1,p,p)
  nu <- matrix(1,p,p)
  tau_sq <- 1
  xi <- 1

  if (verbose) {
    pb <- txtProgressBar(min = 2, max = nIter+nBurnin+1, style = 3,
                         initial = 0)
  }

  for (it in 2:(nBurnin+nIter+1)) {
    expsum <- matrix(0, S_obs, p)
    for(kk in 1:2*K){
      expsum <- expsum + W[, , kk] %*% X %*% Psi[, , kk]
    }
    S <- crossprod(X - expsum)
    # Update Theta
    for (i in 1:p) {
      Theta11 <- Theta[-i, -i]
      Theta12 <- Theta[i, -i]
      Theta21 <- Theta[-i, i]
      Theta22 <- Theta[i, i]
      S11 <- S[-i, -i]
      S12 <- S[i, -i]
      S21 <- S[-i, i]
      S22 <- S[i, i]
      lambda_sq_12 <- lambda_sq[-i, i]
      nu_12 <- nu[-i, i]
      gamma <- rgamma(1, shape = (S_obs/2 + 1), rate = S22/2)
      C <- solve((S22) * solve(Theta11) + diag(1/c(lambda_sq_12 *
                                                     tau_sq)))
      C <- (C + t(C))/2
      beta <- rmvnorm(1, -C %*% S21, C)
      Theta[-i, i] <- Theta[i, -i] <- beta
      Theta[i, i] <- gamma + t(t(beta)) %*% solve(Theta11) %*%
        t(beta)
      rate <- Theta[-i, i]^2/(2 * tau_sq) + 1/nu_12
      lambda_sq_12 <- 1/rgamma((p - 1), shape = 1, rate = rate)
      nu_12 <- 1/rgamma((p - 1), shape = 1, rate = 1 +
                          1/lambda_sq_12)
      lambda_sq[-i, i] <- lambda_sq[i, -i] <- lambda_sq_12
      nu[-i, i] <- nu[i, -i] <- nu_12
    }
    Theta.vec <- Theta[lower.tri(Theta)]
    lambda.sq.vec <- lambda_sq[lower.tri(lambda_sq)]
    rate <- 1/xi + sum(Theta.vec^2/(2 * lambda.sq.vec))
    tau_sq <- 1/rgamma(1, shape = (p * (p - 1)/2 + 1)/2, rate = rate)
    xi <- 1/rgamma(1, shape = 1, rate = 1 + 1/tau_sq)

    if(prior == "ng"){
      # Update omega
      omega <- rgamma(1, shape = c0 + kappa*2*K*p^2, rate = c1 + kappa/2 +
                        sum(alpha[, ,]))
    }

    if(constraint == "symmetric"){
      for(k in 1:2*K){
        for(i in 1:p){
          for(j in 1:p){
            if(i > j){

              if(prior == "ng"){
                # Update alpha
                alpha[i,j,k] <- rgig(n=1,lambda=kappa - 0.5, Psi[i, j, k]^2, kappa*omega)
              }


              # sample initial proposal value
              psi_p <- rnorm(1, Psi[i, j, k], zeta)

              # helper function to check stability conditions

              compdet <- function(psi_p){
                Psi_lik <- Psi[, , k]
                Psi_lik[i, j] <- psi_p
                Psi_lik[j, i] <- psi_p
                Psi_new <- Psi
                Psi_new[, , k] <- Psi_lik

                detsumnew <- matrix(0, S_obs*p, S_obs*p)

                for(kk in 1:2*K){
                  detsumnew <- detsumnew + kronecker.prod(t(Psi_new[, , kk]), W[, , kk])
                }

                return(log(det(diag(S_obs*p) - detsumnew)))
              }

              # faster to compute positive determinant than eigenvalues
              while(is.nan(suppressWarnings(compdet(psi_p)))){
                psi_p <- rnorm(1, Psi[i, j, k], zeta)
              }
              # change Psi in likelihood
              Psi_lik <- Psi[, , k]
              Psi_lik[i, j] <- psi_p
              Psi_lik[j, i] <- psi_p
              Psi_new <- Psi
              Psi_new[, , k] <- Psi_lik

              expsumnew <- matrix(0, S_obs, p)
              expsumold <- matrix(0, S_obs, p)
              detsumnew <- matrix(0, S_obs*p, S_obs*p)
              detsumold <- matrix(0, S_obs*p, S_obs*p)

              for(kk in 1:2*K){
                expsumnew <- expsumnew + W[, , kk] %*% X %*% Psi_new[, , kk]
                expsumold <- expsumold + W[, , kk] %*% X %*% Psi[, , kk]
                detsumnew <- detsumnew + kronecker.prod(t(Psi_new[, , kk]), W[, , kk])
                detsumold <- detsumold + kronecker.prod(t(Psi[, , kk]), W[, , kk])
              }

              # log likelihoods
              logliknew <- log(det(diag(S_obs*p) - detsumnew)) - 0.5*sum(diag(Theta %*% crossprod(X - expsumnew)))
              loglikold <- log(det(diag(S_obs*p) - detsumold)) - 0.5*sum(diag(Theta %*% crossprod(X - expsumold)))

              # log priors
              if(prior == "ng"){
                logpriornew <- dnorm(psi_p, 0, sqrt(2 * (1 / omega) * alpha[i, j, k]), log = TRUE)
                logpriorold <- dnorm(Psi[i, j, k], 0, sqrt(2 * (1 / omega) * alpha[i, j, k]), log = TRUE)
              } else if(prior == "normal"){
                logpriornew <- dnorm(psi_p, 0, 1, log = TRUE)
                logpriorold <- dnorm(Psi[i, j, k], 0, 1, log = TRUE)
              } else{
                logpriornew <- dunif(psi_p, -0.5, 0.5, log = TRUE)
                logpriorold <- dunif(Psi[i, j, k], -0.5, 0.5, log = TRUE)
              }

              if(logliknew + logpriornew > loglikold + logpriorold){
                accpt_rate <- accpt_rate+1
                Psi[i, j, k] <- Psi[j, i, k] <-  psi_p
              } }
            else if(i == j){
              psi_p <- rnorm(1, Psi[i, j, k], zeta)

              # helper function to check stability conditions

              compdet <- function(psi_p){
                Psi_lik <- Psi[, , k]
                Psi_lik[i, j] <- psi_p
                Psi_new <- Psi
                Psi_new[, , k] <- Psi_lik

                detsumnew <- matrix(0, S_obs*p, S_obs*p)

                for(kk in 1:2*K){
                  detsumnew <- detsumnew + kronecker.prod(t(Psi_new[, , kk]), W[, , kk])
                }

                return(log(det(diag(S_obs*p) - detsumnew)))
              }

              while(is.nan(suppressWarnings(compdet(psi_p)))){
                psi_p <- rnorm(1, Psi[i, j, k], zeta)
              }

              # change Psi in likelihood
              Psi_lik <- Psi[, , k]
              Psi_lik[i, j] <- psi_p
              Psi_new <- Psi
              Psi_new[, , k] <- Psi_lik

              expsumnew <- matrix(0, S_obs, p)
              expsumold <- matrix(0, S_obs, p)
              detsumnew <- matrix(0, S_obs*p, S_obs*p)
              detsumold <- matrix(0, S_obs*p, S_obs*p)

              for(kk in 1:2*K){
                expsumnew <- expsumnew + W[, , kk] %*% X %*% Psi_new[, , kk]
                expsumold <- expsumold + W[, , kk] %*% X %*% Psi[, , kk]
                detsumnew <- detsumnew + kronecker.prod(t(Psi_new[, , kk]), W[, , kk])
                detsumold <- detsumold + kronecker.prod(t(Psi[, , kk]), W[, , kk])
              }

              # log likelihoods
              logliknew <- log(det(diag(S_obs*p) - detsumnew)) - 0.5*sum(diag(Theta %*% crossprod(X - expsumnew)))
              loglikold <- log(det(diag(S_obs*p) - detsumold)) - 0.5*sum(diag(Theta %*% crossprod(X - expsumold)))

              # log priors
              logpriornew <- dnorm(psi_p, 0, 0.001, log = TRUE)
              logpriorold <- dnorm(Psi[i, j, k], 0, 0.001, log = TRUE)


              if(logliknew + logpriornew > loglikold + logpriorold){
                # do not count burnin for acceptance rate
                if(it > nBurnin+1){
                  accpt_rate <- accpt_rate+1
                }
                Psi[i, j, k] <-  psi_p
              }
            }
          }
        }
      }
    } else if(constraint == "triangular"){
      for(k in 1:2*K){
        for(i in 1:p){
          for(j in 1:p){
            if(triangular[k] == "lower"){
              if(i > j){

                if(prior == "ng"){
                  # Update alpha
                  alpha[i,j,k] <- rgig(n=1,lambda=kappa - 0.5, Psi[i, j, k]^2, kappa*omega)
                }


                # sample initial proposal value
                psi_p <- rnorm(1, Psi[i, j, k], zeta)

                # helper function to check stability conditions

                compdet <- function(psi_p){
                  Psi_lik <- Psi[, , k]
                  Psi_lik[i, j] <- psi_p
                  Psi_new <- Psi
                  Psi_new[, , k] <- Psi_lik

                  detsumnew <- matrix(0, S_obs*p, S_obs*p)

                  for(kk in 1:2*K){
                    detsumnew <- detsumnew + kronecker.prod(t(Psi_new[, , kk]), W[, , kk])
                  }

                  return(log(det(diag(S_obs*p) - detsumnew)))
                }

                # faster to compute positive determinant than eigenvalues
                while(is.nan(suppressWarnings(compdet(psi_p)))){
                  psi_p <- rnorm(1, Psi[i, j, k], zeta)
                }
                # change Psi in likelihood
                Psi_lik <- Psi[, , k]
                Psi_lik[i, j] <- psi_p
                Psi_new <- Psi
                Psi_new[, , k] <- Psi_lik

                expsumnew <- matrix(0, S_obs, p)
                expsumold <- matrix(0, S_obs, p)
                detsumnew <- matrix(0, S_obs*p, S_obs*p)
                detsumold <- matrix(0, S_obs*p, S_obs*p)

                for(kk in 1:2*K){
                  expsumnew <- expsumnew + W[, , kk] %*% X %*% Psi_new[, , kk]
                  expsumold <- expsumold + W[, , kk] %*% X %*% Psi[, , kk]
                  detsumnew <- detsumnew + kronecker.prod(t(Psi_new[, , kk]), W[, , kk])
                  detsumold <- detsumold + kronecker.prod(t(Psi[, , kk]), W[, , kk])
                }

                # log likelihoods
                logliknew <- log(det(diag(S_obs*p) - detsumnew)) - 0.5*sum(diag(Theta %*% crossprod(X - expsumnew)))
                loglikold <- log(det(diag(S_obs*p) - detsumold)) - 0.5*sum(diag(Theta %*% crossprod(X - expsumold)))

                # log priors
                if(prior == "ng"){
                  logpriornew <- dnorm(psi_p, 0, sqrt(2 * (1 / omega) * alpha[i, j, k]), log = TRUE)
                  logpriorold <- dnorm(Psi[i, j, k], 0, sqrt(2 * (1 / omega) * alpha[i, j, k]), log = TRUE)
                } else if(prior == "normal"){
                  logpriornew <- dnorm(psi_p, 0, 1, log = TRUE)
                  logpriorold <- dnorm(Psi[i, j, k], 0, 1, log = TRUE)
                } else{
                  logpriornew <- dunif(psi_p, -0.5, 0.5, log = TRUE)
                  logpriorold <- dunif(Psi[i, j, k], -0.5, 0.5, log = TRUE)
                }

                if(logliknew + logpriornew > loglikold + logpriorold){
                  accpt_rate <- accpt_rate+1
                  Psi[i, j, k] <-  psi_p
                } }
              else if(i == j){
                psi_p <- rnorm(1, Psi[i, j, k], zeta)

                # helper function to check stability conditions

                compdet <- function(psi_p){
                  Psi_lik <- Psi[, , k]
                  Psi_lik[i, j] <- psi_p
                  Psi_new <- Psi
                  Psi_new[, , k] <- Psi_lik

                  detsumnew <- matrix(0, S_obs*p, S_obs*p)

                  for(kk in 1:2*K){
                    detsumnew <- detsumnew + kronecker.prod(t(Psi_new[, , kk]), W[, , kk])
                  }

                  return(log(det(diag(S_obs*p) - detsumnew)))
                }

                while(is.nan(suppressWarnings(compdet(psi_p)))){
                  psi_p <- rnorm(1, Psi[i, j, k], zeta)
                }

                # change Psi in likelihood
                Psi_lik <- Psi[, , k]
                Psi_lik[i, j] <- psi_p
                Psi_new <- Psi
                Psi_new[, , k] <- Psi_lik

                expsumnew <- matrix(0, S_obs, p)
                expsumold <- matrix(0, S_obs, p)
                detsumnew <- matrix(0, S_obs*p, S_obs*p)
                detsumold <- matrix(0, S_obs*p, S_obs*p)

                for(kk in 1:2*K){
                  expsumnew <- expsumnew + W[, , kk] %*% X %*% Psi_new[, , kk]
                  expsumold <- expsumold + W[, , kk] %*% X %*% Psi[, , kk]
                  detsumnew <- detsumnew + kronecker.prod(t(Psi_new[, , kk]), W[, , kk])
                  detsumold <- detsumold + kronecker.prod(t(Psi[, , kk]), W[, , kk])
                }

                # log likelihoods
                logliknew <- log(det(diag(S_obs*p) - detsumnew)) - 0.5*sum(diag(Theta %*% crossprod(X - expsumnew)))
                loglikold <- log(det(diag(S_obs*p) - detsumold)) - 0.5*sum(diag(Theta %*% crossprod(X - expsumold)))

                # log priors
                logpriornew <- dnorm(psi_p, 0, 0.001, log = TRUE)
                logpriorold <- dnorm(Psi[i, j, k], 0, 0.001, log = TRUE)


                if(logliknew + logpriornew > loglikold + logpriorold){
                  # do not count burnin for acceptance rate
                  if(it > nBurnin+1){
                    accpt_rate <- accpt_rate+1
                  }
                  Psi[i, j, k] <-  psi_p
                }
              } else if(i < j){
                psi_p <- rnorm(1, Psi[i, j, k], zeta)

                # helper function to check stability conditions

                compdet <- function(psi_p){
                  Psi_lik <- Psi[, , k]
                  Psi_lik[i, j] <- psi_p
                  Psi_new <- Psi
                  Psi_new[, , k] <- Psi_lik

                  detsumnew <- matrix(0, S_obs*p, S_obs*p)

                  for(kk in 1:2*K){
                    detsumnew <- detsumnew + kronecker.prod(t(Psi_new[, , kk]), W[, , kk])
                  }

                  return(log(det(diag(S_obs*p) - detsumnew)))
                }

                while(is.nan(suppressWarnings(compdet(psi_p)))){
                  psi_p <- rnorm(1, Psi[i, j, k], zeta)
                }

                # change Psi in likelihood
                Psi_lik <- Psi[, , k]
                Psi_lik[i, j] <- psi_p
                Psi_new <- Psi
                Psi_new[, , k] <- Psi_lik

                expsumnew <- matrix(0, S_obs, p)
                expsumold <- matrix(0, S_obs, p)
                detsumnew <- matrix(0, S_obs*p, S_obs*p)
                detsumold <- matrix(0, S_obs*p, S_obs*p)

                for(kk in 1:2*K){
                  expsumnew <- expsumnew + W[, , kk] %*% X %*% Psi_new[, , kk]
                  expsumold <- expsumold + W[, , kk] %*% X %*% Psi[, , kk]
                  detsumnew <- detsumnew + kronecker.prod(t(Psi_new[, , kk]), W[, , kk])
                  detsumold <- detsumold + kronecker.prod(t(Psi[, , kk]), W[, , kk])
                }

                # log likelihoods
                logliknew <- log(det(diag(S_obs*p) - detsumnew)) - 0.5*sum(diag(Theta %*% crossprod(X - expsumnew)))
                loglikold <- log(det(diag(S_obs*p) - detsumold)) - 0.5*sum(diag(Theta %*% crossprod(X - expsumold)))

                # log priors
                logpriornew <- dnorm(psi_p, 0, 0.001, log = TRUE)
                logpriorold <- dnorm(Psi[i, j, k], 0, 0.001, log = TRUE)
              }
            } else if(triangular[k] == "upper"){
              if(i < j){

                if(prior == "ng"){
                  # Update alpha
                  alpha[i,j,k] <- rgig(n=1,lambda=kappa - 0.5, Psi[i, j, k]^2, kappa*omega)
                }


                # sample initial proposal value
                psi_p <- rnorm(1, Psi[i, j, k], zeta)

                # helper function to check stability conditions

                compdet <- function(psi_p){
                  Psi_lik <- Psi[, , k]
                  Psi_lik[i, j] <- psi_p
                  Psi_new <- Psi
                  Psi_new[, , k] <- Psi_lik

                  detsumnew <- matrix(0, S_obs*p, S_obs*p)

                  for(kk in 1:2*K){
                    detsumnew <- detsumnew + kronecker.prod(t(Psi_new[, , kk]), W[, , kk])
                  }

                  return(log(det(diag(S_obs*p) - detsumnew)))
                }

                # faster to compute positive determinant than eigenvalues
                while(is.nan(suppressWarnings(compdet(psi_p)))){
                  psi_p <- rnorm(1, Psi[i, j, k], zeta)
                }
                # change Psi in likelihood
                Psi_lik <- Psi[, , k]
                Psi_lik[i, j] <- psi_p
                Psi_new <- Psi
                Psi_new[, , k] <- Psi_lik

                expsumnew <- matrix(0, S_obs, p)
                expsumold <- matrix(0, S_obs, p)
                detsumnew <- matrix(0, S_obs*p, S_obs*p)
                detsumold <- matrix(0, S_obs*p, S_obs*p)

                for(kk in 1:2*K){
                  expsumnew <- expsumnew + W[, , kk] %*% X %*% Psi_new[, , kk]
                  expsumold <- expsumold + W[, , kk] %*% X %*% Psi[, , kk]
                  detsumnew <- detsumnew + kronecker.prod(t(Psi_new[, , kk]), W[, , kk])
                  detsumold <- detsumold + kronecker.prod(t(Psi[, , kk]), W[, , kk])
                }

                # log likelihoods
                logliknew <- log(det(diag(S_obs*p) - detsumnew)) - 0.5*sum(diag(Theta %*% crossprod(X - expsumnew)))
                loglikold <- log(det(diag(S_obs*p) - detsumold)) - 0.5*sum(diag(Theta %*% crossprod(X - expsumold)))

                # log priors
                if(prior == "ng"){
                  logpriornew <- dnorm(psi_p, 0, sqrt(2 * (1 / omega) * alpha[i, j, k]), log = TRUE)
                  logpriorold <- dnorm(Psi[i, j, k], 0, sqrt(2 * (1 / omega) * alpha[i, j, k]), log = TRUE)
                } else if(prior == "normal"){
                  logpriornew <- dnorm(psi_p, 0, 1, log = TRUE)
                  logpriorold <- dnorm(Psi[i, j, k], 0, 1, log = TRUE)
                } else{
                  logpriornew <- dunif(psi_p, -0.5, 0.5, log = TRUE)
                  logpriorold <- dunif(Psi[i, j, k], -0.5, 0.5, log = TRUE)
                }

                if(logliknew + logpriornew > loglikold + logpriorold){
                  accpt_rate <- accpt_rate+1
                  Psi[i, j, k] <-  psi_p
                } }
              else if(i == j){
                psi_p <- rnorm(1, Psi[i, j, k], zeta)

                # helper function to check stability conditions

                compdet <- function(psi_p){
                  Psi_lik <- Psi[, , k]
                  Psi_lik[i, j] <- psi_p
                  Psi_new <- Psi
                  Psi_new[, , k] <- Psi_lik

                  detsumnew <- matrix(0, S_obs*p, S_obs*p)

                  for(kk in 1:2*K){
                    detsumnew <- detsumnew + kronecker.prod(t(Psi_new[, , kk]), W[, , kk])
                  }

                  return(log(det(diag(S_obs*p) - detsumnew)))
                }

                while(is.nan(suppressWarnings(compdet(psi_p)))){
                  psi_p <- rnorm(1, Psi[i, j, k], zeta)
                }

                # change Psi in likelihood
                Psi_lik <- Psi[, , k]
                Psi_lik[i, j] <- psi_p
                Psi_new <- Psi
                Psi_new[, , k] <- Psi_lik

                expsumnew <- matrix(0, S_obs, p)
                expsumold <- matrix(0, S_obs, p)
                detsumnew <- matrix(0, S_obs*p, S_obs*p)
                detsumold <- matrix(0, S_obs*p, S_obs*p)

                for(kk in 1:2*K){
                  expsumnew <- expsumnew + W[, , kk] %*% X %*% Psi_new[, , kk]
                  expsumold <- expsumold + W[, , kk] %*% X %*% Psi[, , kk]
                  detsumnew <- detsumnew + kronecker.prod(t(Psi_new[, , kk]), W[, , kk])
                  detsumold <- detsumold + kronecker.prod(t(Psi[, , kk]), W[, , kk])
                }

                # log likelihoods
                logliknew <- log(det(diag(S_obs*p) - detsumnew)) - 0.5*sum(diag(Theta %*% crossprod(X - expsumnew)))
                loglikold <- log(det(diag(S_obs*p) - detsumold)) - 0.5*sum(diag(Theta %*% crossprod(X - expsumold)))

                # log priors
                logpriornew <- dnorm(psi_p, 0, 0.001, log = TRUE)
                logpriorold <- dnorm(Psi[i, j, k], 0, 0.001, log = TRUE)


                if(logliknew + logpriornew > loglikold + logpriorold){
                  # do not count burnin for acceptance rate
                  if(it > nBurnin+1){
                    accpt_rate <- accpt_rate+1
                  }
                  Psi[i, j, k] <-  psi_p
                }
              } else if(i > j){
                psi_p <- rnorm(1, Psi[i, j, k], zeta)

                # helper function to check stability conditions

                compdet <- function(psi_p){
                  Psi_lik <- Psi[, , k]
                  Psi_lik[i, j] <- psi_p
                  Psi_new <- Psi
                  Psi_new[, , k] <- Psi_lik

                  detsumnew <- matrix(0, S_obs*p, S_obs*p)

                  for(kk in 1:2*K){
                    detsumnew <- detsumnew + kronecker.prod(t(Psi_new[, , kk]), W[, , kk])
                  }

                  return(log(det(diag(S_obs*p) - detsumnew)))
                }

                while(is.nan(suppressWarnings(compdet(psi_p)))){
                  psi_p <- rnorm(1, Psi[i, j, k], zeta)
                }

                # change Psi in likelihood
                Psi_lik <- Psi[, , k]
                Psi_lik[i, j] <- psi_p
                Psi_new <- Psi
                Psi_new[, , k] <- Psi_lik

                expsumnew <- matrix(0, S_obs, p)
                expsumold <- matrix(0, S_obs, p)
                detsumnew <- matrix(0, S_obs*p, S_obs*p)
                detsumold <- matrix(0, S_obs*p, S_obs*p)

                for(kk in 1:2*K){
                  expsumnew <- expsumnew + W[, , kk] %*% X %*% Psi_new[, , kk]
                  expsumold <- expsumold + W[, , kk] %*% X %*% Psi[, , kk]
                  detsumnew <- detsumnew + kronecker.prod(t(Psi_new[, , kk]), W[, , kk])
                  detsumold <- detsumold + kronecker.prod(t(Psi[, , kk]), W[, , kk])
                }

                # log likelihoods
                logliknew <- log(det(diag(S_obs*p) - detsumnew)) - 0.5*sum(diag(Theta %*% crossprod(X - expsumnew)))
                loglikold <- log(det(diag(S_obs*p) - detsumold)) - 0.5*sum(diag(Theta %*% crossprod(X - expsumold)))

                # log priors
                logpriornew <- dnorm(psi_p, 0, 0.001, log = TRUE)
                logpriorold <- dnorm(Psi[i, j, k], 0, 0.001, log = TRUE)
              }
            }
          }
        }
      }
    }

    # save only iterations higher than burning + 1 because that is the point the burnin stops (and we start at iteration 2)
    if(prior == "ng"){
      if (it > nBurnin+1) {
        Results$Theta[, , it-(nBurnin+1)] <- Theta
        Results$lambda_sq[, , it-(nBurnin+1)] <- lambda_sq
        Results$tau_sq[it-(nBurnin+1)] <- tau_sq
        Results$omega[it-(nBurnin+1)] <- omega
        Results$alpha[, , , it-(nBurnin+1)] <- alpha
        Results$Psi[, , , it-(nBurnin+1)] <- Psi
      }
    } else{
      if (it > nBurnin+1) {
        Results$Theta[, , it-(nBurnin+1)] <- Theta
        Results$lambda_sq[, , it-(nBurnin+1)] <- lambda_sq
        Results$tau_sq[it-(nBurnin+1)] <- tau_sq
        Results$Psi[, , , it-(nBurnin+1)] <- Psi
      }
    }

    if (verbose) {
      setTxtProgressBar(pb, it)
    }

  }
  if (verbose) {
    close(pb)
  }
  # we have p(p+1) and 2kp^2 steps for a symmetric and
  # for a triangular constraint respectively
  if(constraint == "symmetric"){
    accpt_rate <- accpt_rate/((nIter)*(p*(p+1)))
  } else if(constraint == "triangular"){
    accpt_rate <- accpt_rate/((nIter)*2*K*p^2)
  }

  Results$accpt_rate <- accpt_rate
  return(Results)
}
