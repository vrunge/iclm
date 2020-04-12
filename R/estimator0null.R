#' probaEstimation
#' @description Estimation of 4 probabilities
#' @param X1 1st data column
#' @param X2 2nd data column
#' @param sigma sd
#' @param n data length
#' @param nbSimu number of simulations
#' @return a vector of 4 estimated probabilities
probaEstimation <- function(X1, X2, sigma, n, nbSimu)
{
  counts <- rep(0,4)
  X1carre <- c(X1%*%X1)
  X2carre <- c(X2%*%X2)
  X12 <- c(X1%*%X2)

  for(i in 1:nbSimu)
  {
    epsilon <- rnorm(n, sd = sigma)
    A1 <- X2carre*c(X1%*%epsilon) - X12*c(X2%*%epsilon)
    A2 <- - X12*c(X1%*%epsilon) + X1carre*c(X2%*%epsilon)
    B1 <- c(X1%*%epsilon)
    B2 <- c(X2%*%epsilon)
    if(A1 >= 0 && A2 >= 0){counts[1] <- counts[1] + 1}
    if(B1 <= 0 && B2 <= 0){counts[2] <- counts[2] + 1}
    if(A1 <= 0 && B2 >= 0){counts[3] <- counts[3] + 1}
    if(B1 >= 0 && A2 <= 0){counts[4] <- counts[4] + 1}
  }
  return(counts/nbSimu)
}


#' densitybeta1
#' @description density of the estimator for beta1
#' @param X1 1st data column
#' @param X2 2nd data column
#' @param sigma sd
#' @param p the 4 probabilities obtained with the probaEstimation function
#' @param h step size
#' @param max max value
#' @return x y and d0 density at 0
densitybeta1 <- function(X1, X2, sigma, p, h = 0.01, max = 2)
{
  X1carre <- c(X1%*%X1)
  X2carre <- c(X2%*%X2)
  X12 <- c(X1%*%X2)
  Delta <- X1carre*X2carre - X12^2
  x <- seq(from = 0, to = max, by = h)
  y <- rep(0, length(x))
  P0 <- p[2] + p[3]
  P1 <- p[1]
  P2 <- p[4]
  y <- 2*P1*dnorm(x = x, mean = 0, sd = sigma*sqrt(X2carre/Delta)) + 2*P2*dnorm(x = x, mean = 0, sd = sigma/sqrt(X1carre))

  return(list(x = x,y = y, d0 = P0))
}


#' beta1estimator0
#' @description simulations for beta1
#' @param X1 1st data column
#' @param X2 2nd data column
#' @param sigma sd
#' @param n data length
#' @param nbSimu number of simulations
#' @return A vector of beta1 simulated
beta1estimator0 <- function(X1, X2, sigma, n, nbSimu)
{
  beta1 <- rep(0, nbSimu)
  X1carre <- c(X1%*%X1)
  X2carre <- c(X2%*%X2)
  X12 <- c(X1%*%X2)
  Delta <- X1carre*X2carre - X12^2

  for(i in 1:nbSimu)
  {
    epsilon <- rnorm(n, sd = sigma)
    A1 <- X2carre*c(X1%*%epsilon) - X12*c(X2%*%epsilon)
    A2 <- - X12*c(X1%*%epsilon) + X1carre*c(X2%*%epsilon)
    B1 <- c(X1%*%epsilon)
    B2 <- c(X2%*%epsilon)

    if(A1 >= 0 && A2 >= 0){beta1[i] <- A1/Delta}
    if(B1 <= 0 && B2 <= 0){beta1[i] <- 0}
    if(A1 <= 0 && B2 >= 0){beta1[i] <- 0}
    if(B1 >= 0 && A2 <= 0){beta1[i] <- B1/X1carre}
  }
  return(beta1)
}
