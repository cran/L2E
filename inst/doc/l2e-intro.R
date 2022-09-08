## ----echo=FALSE, cache=FALSE--------------------------------------------------
set.seed(12345)
knitr::opts_chunk$set(
  cache=TRUE,
  comment = '', 
  fig.width = 5, 
  fig.height = 5,
  fig.align='center'
)

## -----------------------------------------------------------------------------
set.seed(12345)
n <- 200
tau <- 1
x <- seq(-2.5, 2.5, length.out=n)
f <- x^3
y <- f + (1/tau)*rnorm(n)
plot(x, y, pch=16, col='gray')
lines(x, f, lwd=3)

## -----------------------------------------------------------------------------
library(L2E)
library(isotone)
tau <- 1/mad(y)
b <- y
# LS method
iso <- gpava(1:n, y)$x
# MM method
sol_mm <- L2E_isotonic(y, b, tau, method = "MM")
# PG method
sol_pg <- L2E_isotonic(y, b, tau, method = 'PG')


# Plots
plot(x, y, pch=16, col='gray')
lines(x, f, lwd=3)
lines(x, iso, col='blue', lwd=3) ## LS
lines(x, sol_mm$beta, col='red', lwd=3) ## MM
lines(x, sol_pg$beta, col='green', lwd=3) ## PG
legend("bottomright", legend = c("LS", "MM", "PG"), col = c('blue','red', 'green'), lwd=3)

## -----------------------------------------------------------------------------
num <- 20
ix <- 1:num
y[45 + ix] <- 14 + rnorm(num)
plot(x, y, pch=16, col='gray')
lines(x, f, lwd=3)

## -----------------------------------------------------------------------------

tau <- 1/mad(y)
b <- y
# LS method
iso <- gpava(1:n, y)$x
# MM method
sol_mm <- L2E_isotonic(y, b, tau, method = "MM")
# PG method
sol_pg <- L2E_isotonic(y, b, tau, method = 'PG')


# Plots
plot(x, y, pch=16,  col='gray')
lines(x, f, lwd=3)
lines(x, iso, col='blue', lwd=3) ## LS
lines(x, sol_mm$beta, col='red', lwd=3) ## MM
lines(x, sol_pg$beta, col='green', lwd=3) ## PG
legend("bottomright", legend = c("LS", "MM", "PG"), col = c('blue','red', 'green'), lwd=3)

## -----------------------------------------------------------------------------
set.seed(12345)
n <- 300
tau <- 1
x <- seq(-2, 2, length.out=n)
f <- x^4 + x
y <- f + (1/tau) * rnorm(n)
plot(x, y, pch=16, col='gray', cex=0.8)
lines(x, f, col='black', lwd=3)

## -----------------------------------------------------------------------------
library(cobs)
tau <- 1/mad(y)
b <- y
## LS method
cvx <- fitted(conreg(y, convex=TRUE))
## MM method
sol_mm <- L2E_convex(y, b, tau, method = "MM")
## PG method
sol_pg <- L2E_convex(y, b, tau, method = 'PG')

plot(x, y, pch=16, col='gray')
lines(x, f, lwd=3)
lines(x, cvx, col='blue', lwd=3) ## LS
lines(x, sol_mm$beta, col='red', lwd=3) ## MM
lines(x, sol_pg$beta, col='green', lwd=3) ## PG
legend("bottomright", legend = c("LS", "MM", "PG"), col = c('blue','red', 'green'), lwd=3)

## -----------------------------------------------------------------------------
num <- 50
ix <- 1:num
y[45 + ix] <- 14 + rnorm(num)
 
plot(x, y, pch=16, col='gray', cex=0.8)
lines(x, f, col='black', lwd=3)

## -----------------------------------------------------------------------------
tau <- 1/mad(y)
b <- y
## LS method
cvx <- fitted(conreg(y, convex=TRUE))
## MM method
sol_mm <- L2E_convex(y, b, tau, method = "MM")
## PG method
sol_pg <- L2E_convex(y, b, tau, method = 'PG')

plot(x, y, pch=16, col='gray')
lines(x, f, lwd=3)
lines(x, cvx, col='blue', lwd=3) ## LS
lines(x, sol_mm$beta, col='red', lwd=3) ## MM
lines(x, sol_pg$beta, col='green', lwd=3) ## PG
legend("bottomright", legend = c("LS", "MM", "PG"), col = c('blue','red', 'green'), lwd=3)


## ----echo=FALSE, cache=FALSE--------------------------------------------------
set.seed(12345)
knitr::opts_chunk$set(
  cache=TRUE,
  comment = '', 
  fig.width = 5, 
  fig.height = 4,
  fig.align='center'
)

## ---- title="star"------------------------------------------------------------
library(robustbase)
data(starsCYG)
plot(starsCYG)

## -----------------------------------------------------------------------------
y <- starsCYG[, "log.light"]
x <- starsCYG[, "log.Te"]
X0 <- cbind(rep(1, length(y)), x)

# LS method
mle <- lm(log.light ~ log.Te, data = starsCYG)
r_lm <- y - X0 %*% mle$coefficients

# L2E+MM method
tau <- 1/mad(y)
b <- c(0, 0)
# Fit the regression model
sol_mm <- L2E_multivariate(y, X0, b, tau, method="MM")
l2e_fit_mm <- X0 %*% sol_mm$beta
# compute limit weights
r_mm <- y - l2e_fit_mm


data <- data.frame(x, y, l2e_fit_mm)


d_lines <- data.frame(int = c(sol_mm$beta[1], mle$coefficients[1]),
                      sl = c(sol_mm$beta[2], mle$coefficients[2]),
                      col = c("red", "blue"),
                      lty = c("solid", "dashed"),
                      method = c("L2E", "LS"))

ltys <- as.character(d_lines$lty)
names(ltys) <- as.character(d_lines$lty)

cols <- as.character(d_lines$col)
cols <- cols[order(as.character(d_lines$lty))]
method <-  as.character(d_lines$method)


library(ggplot2)
library(latex2exp)
p <- ggplot() +
  geom_point(data = data, aes(x, y), size=2.5) + ylim(2, 6.5)+
  geom_abline(data = d_lines[d_lines$col == "red", ], 
              aes(intercept = int, slope = sl, lty = lty), color = "red", size=1) +
  geom_abline(data = d_lines[d_lines$col == "blue", ], 
              aes(intercept = int, slope = sl, lty = lty), color = "blue", size=1) +
  scale_linetype_manual(name = "Method", values = ltys, breaks=c("dashed", "solid"),
                        labels = c("LS  ", 
                                   expression(L[2]~E)),
                        guide = guide_legend(override.aes = list(colour = cols), legend=method))+
  theme_bw()

print(p)

## -----------------------------------------------------------------------------
w <- as.vector(exp(-0.5* (sol_mm$tau*r_mm)**2 ))
data <- data.frame(x, y, l2e_fit_mm, w)
ggplot(data, aes(x=log10(w))) + geom_histogram()+
  labs(
    y="Count", x=expression(log[10]~'(w)'))+theme_bw()

## -----------------------------------------------------------------------------
outlier_mm <- rep("yes", length(y))
for (k in 1:length(y)) {
  if(w[k]>1e-5) # the threshold value can range from 1e-3 to 1e-14 according to the histogram
    outlier_mm[k] <- "no"
}
outlier_mm  <- factor(outlier_mm , levels=c("yes", "no"))
data <- data.frame(x, y, l2e_fit_mm, outlier_mm)


p+
  geom_point(data = data, aes(x, y, color=outlier_mm), size=2.5) +
  scale_color_manual(values = c(2,1), name="Outlier")+
  labs(
    y="Light Intensity", x="Temperature")+ theme_bw() 


