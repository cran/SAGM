make_weights <- function(n){
  W_BA <- matrix(0, n, n)
  W_AB <- matrix(0, n, n)
  for(i in 1:n){
    # W_BA
    if(i %% 2 != 0){
      # uneven number of n
      if(n %% 2 != 0){
        if(i == n){
          W_BA[i,i-1] <- 1
        } else if(i == 1){
          W_BA[i,2] <- 1
        } else{
          W_BA[i,c(i-1, i+1)] <- 0.5
        }
      }
      # even number of n
      else{
        if(i == 1){
          W_BA[i,2] <- 1
        } else{
          W_BA[i,c(i-1, i+1)] <- 0.5
        }
      }
    # W_BA
    } else{
      # uneven number of n
      if(n %% 2 != 0){
        W_AB[i,c(i-1, i+1)] <- 0.5
      }
      # even number of n
      else{
        if(i == n){
          W_AB[i,i-1] <- 1
        }
        else{
          W_AB[i,c(i-1, i+1)] <- 0.5
        }
      }
    }
  }
  return(list(W_BA = W_BA, W_AB = W_AB))
}

