## ----echo=FALSE, cache=FALSE--------------------------------------------------
set.seed(12345)
knitr::opts_chunk$set(
  cache=TRUE,
  comment = '', 
  fig.width = 6, 
  fig.height = 6,
  fig.align='center'
)

## ---- title="bank"------------------------------------------------------------
library(L2E)
y <- bank$y
X <- as.matrix(bank[,1:13])
X0 <- as.matrix(cbind(rep(1,length(y)), X))
tau_initial <- 1/mad(y)
beta_initial <- matrix(0, 14, 1)

## -----------------------------------------------------------------------------
sol <- l2e_regression(y, X0, tau_initial, beta_initial)

## -----------------------------------------------------------------------------
betaEstimate <- sol$beta
tauEstimate <- sol$tau

r <- y - X0 %*% betaEstimate
outliers <- which(abs(r) > 3/tauEstimate)
l2e_fit <- X0 %*% betaEstimate
plot(y, l2e_fit, ylab='Predicted values', pch=16, cex=0.8)
points(y[outliers], l2e_fit[outliers], pch=16, col='blue', cex=0.8)

## -----------------------------------------------------------------------------
set.seed(12345)
n <- 200
tau <- 1
x <- seq(-2.5, 2.5, length.out=n)
f <- x^3
y <- f + (1/tau)*rnorm(n)
plot(x, y, pch=16, col='gray', cex=0.8)
lines(x, f, col='black', lwd=3)

## -----------------------------------------------------------------------------
tau_initial <- 1
beta_initial <- y
sol <- l2e_regression_isotonic(y, beta_initial, tau_initial)
isotonic_LTE <- sol$beta

plot(x, y, pch=16, col='gray', cex=0.8)
lines(x, f, col='black', lwd=3)
isotonic_MLE <- gpava(1:n, y)$x
lines(x, isotonic_MLE, col='red', lwd=3)
lines(x, isotonic_LTE, col='blue', lwd=3)

## -----------------------------------------------------------------------------
ix <- 0:9
y[45 + ix] <- 14 + rnorm(10)

plot(x, y, pch=16, col='gray', cex=0.8)
lines(x, f, col='black', lwd=3)

## -----------------------------------------------------------------------------
plot(x, y, pch=16, col="gray", cex=0.8)
lines(x, f, col='black', lwd=3)
isotonic_MLE <- gpava(1:n, y)$x
lines(x, isotonic_MLE, col='red', lwd=3)
sol <- l2e_regression_isotonic(y, beta_initial, tau_initial)
isotonic_LTE <- sol$beta
lines(x, isotonic_LTE, col='blue', lwd=3)

## -----------------------------------------------------------------------------
set.seed(12345)
n <- 200
tau <- 1
x <- seq(-2, 2, length.out=n)
f <- x^4 + x
y <- f + (1/tau) * rnorm(n)
plot(x, y, pch=16, col='gray', cex=0.8)
lines(x, f, col='black', lwd=3)

## -----------------------------------------------------------------------------
tau_initial <- 1
beta_initial <- y
sol <- l2e_regression_convex(y, beta_initial, tau_initial)
convex_LTE <- sol$beta

plot(x, y, pch=16, col='gray', cex=0.8)
lines(x, f, col='black', lwd=3)
convex_MLE <- fitted(cobs::conreg(y, convex=TRUE))
lines(x, convex_MLE, col='red', lwd=3)
lines(x, convex_LTE, col='blue', lwd=3)

## -----------------------------------------------------------------------------
ix <- 0:9
y[45 + ix] <- 14 + rnorm(10)
 
plot(x, y, pch=16, col='gray', cex=0.8)
lines(x, f, col='black', lwd=3)

## -----------------------------------------------------------------------------
plot(x, y, pch=16, col='gray', cex=0.8)
lines(x, f, col='black', lwd=3)
convex_MLE <- fitted(cobs::conreg(y, convex=TRUE))
lines(x, convex_MLE, col='red', lwd=3)
sol <- l2e_regression_convex(y, beta_initial, tau_initial)
convex_LTE <- sol$beta
lines(x, convex_LTE, col='blue', lwd=3)

