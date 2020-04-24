# nnlm : Non-negative linear model


Installation of the package from github
```r
devtools::install_github("vrunge/iclm")
```
## Data generation

We chose parameters: columns X1 and X2

```r
nbSimu <- 10^5
n <- 10
sigma <- 3
X1 <- runif(n, min=0, max = 1)
X2 <- runif(n, min=0, max = 5)
```

## Proba estimation


```r
#X1 <- X1 - mean(X1)
#X2 <- X2 - mean(X2)

p <- probaEstimation(X1, X2, sigma, n, nbSimu = nbSimu)
```


## Comparison simulation VS estimator 

```r
step <- 0.001
h <- densitybeta1(X1, X2, sigma, p, h = step, max = 10)
b1 <- beta1estimator0(X1, X2, sigma, n, nbSimu = nbSimu)
```


## Plots and tests

```r
his <- hist(b1$beta1[b1$beta1>0], breaks = 50, probability = TRUE)
par(new=TRUE)
plot(h$x,h$y/sum(step*h$y), type = 'l', xlim = c(0,max(his$breaks)), ylim = c(0, his$density[1]))


#verif density = 1
his$density%*%diff(his$breaks)
sum(step*h$y) + h$d0

#comparaison dirac en 0
sum(b1$beta1 == 0)/nbSimu
p[2] + p[3]
```

