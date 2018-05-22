library("tseries") 
library("timsac") # Akaike's final prediction error for AR
library("forecast")
library("MASS")

my.spectrum <- function(phi.of.b, theta.of.b, variance=1)
{
  p <- length(phi.of.b)
  q <- length(theta.of.b)
  omega <- seq(from=0, to=pi, by=.001)
  phi.of.e.minus.i.omega <- 1
  phi.of.e.i.omega <- 1
  if(p>1)
  {   for(i in 2:p)
  {
    phi.of.e.minus.i.omega <-  phi.of.e.minus.i.omega + phi.of.b[i]*exp(complex(imaginary = -(i-1))*omega)
    phi.of.e.i.omega <-  phi.of.e.i.omega + phi.of.b[i]*exp(complex(imaginary = (i-1))*omega)
  }
  }
  theta.of.e.minus.i.omega <- 1
  theta.of.e.i.omega <- 1
  if(q>1)
  {
    for(i in 2:q)
    {
      theta.of.e.minus.i.omega <-  theta.of.e.minus.i.omega + theta.of.b[i]*exp(complex(imaginary = -(i-1))*omega)
      theta.of.e.i.omega <-  theta.of.e.i.omega + theta.of.b[i]*exp(complex(imaginary = (i-1))*omega)
    }
  }
  my.spectrum <- (variance/(2*pi))*Re(theta.of.e.minus.i.omega*theta.of.e.i.omega/(phi.of.e.minus.i.omega*phi.of.e.i.omega))
  plot(omega, 10*log10(my.spectrum), ylab="spectrum (in decibels)", type="l")   
}



data1 <- read.csv("proj1.txt", header = FALSE)
data2 <- read.csv("proj2.txt", header = FALSE)
data3 <- read.csv("proj3.txt", header = FALSE)
data4 <- read.csv("deposits.txt", header = FALSE)

auto.arima(data1, trace = TRUE, stepwise = FALSE)
auto.arima(data2, trace = TRUE, stepwise = FALSE, max.p = 5, max.q = 5)
auto.arima(data3, trace = TRUE, stepwise = FALSE)
auto.arima(data4, trace = TRUE, stepwise = FALSE)
auto.arima(data4t, trace = TRUE, stepwise = FALSE)

auto.arima(data2, trace = TRUE, stepwise = FALSE)


data1 <- as.ts(data1)
data1t <- as.ts(data1[1:564])
data1v <- as.ts(data1[565:584])
data2 <- as.ts(data2)
data2t <- as.ts(data2[1:780])
data2v <- as.ts(data2[781:800])
data3 <- as.ts(data3)
data3t <- as.ts(data3[1:480])
data3v <- as.ts(data3[481:500])
data4 <- as.ts(data4)
data4t <- as.ts(data4[1:604])
data4v <- as.ts(data4[605:624])

################################################################################
# proj1
par(mfrow = c(3,1), mar = c(4,4,3,1))
plot(data1t, type="b", pch = 20, main = "Scatterplot", ylab = "Data")
acf(data1t, type = c("correlation"), main = "ACF Plot")
acf(data1t, type = c("partial"), main = "PACF Plot")
par(mfrow = c(1,1), mar = c(4,4.5,3,1))
spectrum(data1t, method = "ar")
spectrum(data1t, method = "pgram", ylab = expression(paste("I(", omega, ")")), main = "Periodogram")

adf.test(data1t)

ar(data1t, aic=FALSE, order.max = 4, method=c("yule-walker"))
ar(data1t, aic=FALSE, order.max = 4, method=c("burg"))
ar(data1t, aic=FALSE, order.max = 4, method=c("ols"))
ar(as.numeric(data1t), aic=FALSE, order.max = 4, method=c("mle"))

ar(data1t, method=c("yule-walker"))
ar(data1t, method=c("burg"))
ar(data1t, method=c("ols"))
ar(data1t, method=c("mle"))

fit1.4 <- arima(data1t, order = c(4, 0, 0))
fit1.5 <- arima(data1t, order = c(5, 0, 0))
fit1.6 <- arima(data1t, order = c(6, 0, 0))
fit1.7 <- arima(data1t, order = c(7, 0, 0))

fit1.4
fit1.5
fit1.6
fit1.7

tsdiag(fit1.4)
tsdiag(fit1.5)
tsdiag(fit1.6)
tsdiag(fit1.7)

sum((as.numeric(data1v) - predict(fit1.4, n.ahead = 20, se.fit = FALSE))^2)
sum((as.numeric(data1v) - predict(fit1.5, n.ahead = 20, se.fit = FALSE))^2)
sum((as.numeric(data1v) - predict(fit1.6, n.ahead = 20, se.fit = FALSE))^2)
sum((as.numeric(data1v) - predict(fit1.7, n.ahead = 20, se.fit = FALSE))^2)

plot(as.numeric(data1v), type = "b")
points(as.numeric(predict(fit1.4, n.ahead = 20, se.fit = FALSE)), type = "b", col = "red")
points(as.numeric(predict(fit1.5, n.ahead = 20, se.fit = FALSE)), type = "b", col = "yellow")
points(as.numeric(predict(fit1.6, n.ahead = 20, se.fit = FALSE)), type = "b", col = "green")
points(as.numeric(predict(fit1.7, n.ahead = 20, se.fit = FALSE)), type = "b", col = "blue")

fit1.4.1 <- arima(data1t, order = c(4, 0, 1))
fit1.4.2 <- arima(data1t, order = c(4, 0, 2))
fit1.4.3 <- arima(data1t, order = c(4, 0, 3))
fit1.4.4 <- arima(data1t, order = c(4, 0, 4))
fit1.4.5 <- arima(data1t, order = c(4, 0, 5))
fit1.4.6 <- arima(data1t, order = c(4, 0, 6))

tsdiag(fit1.4.1)
tsdiag(fit1.4.2)
tsdiag(fit1.4.3)
tsdiag(fit1.4.4)
tsdiag(fit1.4.5)
tsdiag(fit1.4.6)

fit1.4.3
fit1.4.4
fit1.4.5
fit1.4.6

sum((as.numeric(data1v) - predict(fit1.4.3, n.ahead = 20, se.fit = FALSE))^2)
sum((as.numeric(data1v) - predict(fit1.4.4, n.ahead = 20, se.fit = FALSE))^2)

plot(as.numeric(data1v), type = "b")
points(as.numeric(predict(fit1.4.3, n.ahead = 20, se.fit = FALSE)), type = "b", col = "green")
points(as.numeric(predict(fit1.4.4, n.ahead = 20, se.fit = FALSE)), type = "b", col = "blue")

fit1.5.1 <- arima(data1t, order = c(5, 0, 1))
fit1.5.2 <- arima(data1t, order = c(5, 0, 2))
fit1.5.3 <- arima(data1t, order = c(5, 0, 3))
fit1.5.4 <- arima(data1t, order = c(5, 0, 4))
fit1.5.5 <- arima(data1t, order = c(5, 0, 5))
fit1.5.6 <- arima(data1t, order = c(5, 0, 6))
tsdiag(fit1.5.1)
tsdiag(fit1.5.2)
tsdiag(fit1.5.3)
tsdiag(fit1.5.4)
tsdiag(fit1.5.5)
tsdiag(fit1.5.6)
fit1.5.1 
fit1.5.2
fit1.5.3 
fit1.5.4 
fit1.5.5 
fit1.5.6

par(mfrow=c(2,1))
spectrum(data1,ylab = expression(paste("I(", omega, ")")), main = "Periodogram")
my.spectrum(phi.of.b=c(1, -final1$coef[1:5]), theta.of.b=c(1), variance=final1$sigma2)
par(mfrow=c(1,1))

BIC(arima(data1, order = c(5, 0, 0)))

arima(data1, order = c(5, 0, 0), include.mean = FALSE)


final1 <- arima(data1, order = c(5, 0, 0), include.mean = FALSE)
my.preds <- predict(final1, n.ahead = 13)
my.lower.limits  <- my.preds$pred - 1.96*my.preds$se
my.upper.limits <- my.preds$pred + 1.96*my.preds$se

par(mfrow = c(2,1), mar = c(3.9,4,0.5,1))
plot(1:584, data1[1:584], xlim = c(1, 598), type = "l", xlab = "Time", ylab = "Data")
points(585:597, as.numeric(my.preds$pred), type = "b", pch = 4)
points(585:597, as.numeric(my.lower.limits), type = "l", lty = 2)
points(585:597, as.numeric(my.upper.limits), type = "l", lty = 2)
plot(564:584, data1[564:584], xlim = c(564, 598), ylim = c(-300,490), type = "l", xlab = "Time", ylab = "Data")
points(585:597, as.numeric(my.preds$pred), type = "b", pch = 4)
points(584:597, c(data1[584],as.numeric(my.lower.limits)), type = "l",lty = 2)
points(584:597, c(data1[584],as.numeric(my.upper.limits)), type = "l",lty = 2)
par(mfrow = c(1,1))

qqnorm(final1$residuals)
tsdiag(final1)
final1

Mod(polyroot(c(1,-2.2609, 1.0560, 1.3543, -1.6397, 0.5187)))
Arg(polyroot(c(1,-2.2609, 1.0560, 1.3543, -1.6397, 0.5187)))
plot(polyroot(c(1,-2.2609, 1.0560, 1.3543, -1.6397, 0.5187)),
     xlim=c(-1,1), ylim=c(-1,1), 
     xlab="Real axis", ylab="Imaginary axis",
     main="Roots of Autoregressive Polynomial", asp=1)
symbols(0, 0, circles=1, add=T, inches=F, col=5) 
spectrum(data1,ylab = expression(paste("I(", omega, ")")), main = "Periodogram")

plot(data1[50:150])
par(mfrow=c(2,1))
spectrum(data1,ylab = expression(paste("I(", omega, ")")), main = "Periodogram")
my.spectrum(phi.of.b=c(1, -fit1.4.3$coef[1:4]), theta.of.b=c(1, fit1.4.3$coef[5:7]), variance=fit1.4.3$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
spectrum(data1,ylab = expression(paste("I(", omega, ")")), main = "Periodogram")
my.spectrum(phi.of.b=c(1, -fit1.4.4$coef[1:4]), theta.of.b=c(1, fit1.4.4$coef[5:8]), variance=fit1.4.4$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
spectrum(data1,ylab = expression(paste("I(", omega, ")")), main = "Periodogram")
my.spectrum(phi.of.b=c(1, -fit1.5$coef[1:5]), theta.of.b=c(1), variance=fit1.5$sigma2)
par(mfrow=c(1,1))


ARMAacf(ar=c(final1$coef[1:5]), lag.max=10, pacf=FALSE)
ARMAacf(ar=c(final1$coef[1:5]), lag.max=10, pacf=TRUE)

# 
# my.preds <- predict(fit1.5, n.ahead = 59)
# my.lower.limits  <- my.preds$pred - 1.96*my.preds$se
# my.upper.limits <- my.preds$pred + 1.96*my.preds$se
# 
# plot(as.numeric(data1v), type = "b")
# points(as.numeric(my.preds$pred), type = "b", col = "red")
# points(as.numeric(my.lower.limits), type = "l")
# points(as.numeric(my.upper.limits), type = "l")
# 
# fit1.55 <- arima(data1[1:571], order = c(5, 0, 0))
# my.preds <- predict(fit1.55, n.ahead = 13)
# my.lower.limits  <- my.preds$pred - 1.96*my.preds$se
# my.upper.limits <- my.preds$pred + 1.96*my.preds$se
# plot(as.numeric(data1[572:584]), type = "b", ylim = c(-350, 300))
# points(as.numeric(my.preds$pred), type = "b", col = "red")
# points(as.numeric(my.lower.limits), type = "l")
# points(as.numeric(my.upper.limits), type = "l")
# 
# fit1.77 <- arima(data1[1:571], order = c(7, 0, 0))
# my.preds <- predict(fit1.77, n.ahead = 13)
# my.lower.limits  <- my.preds$pred - 1.96*my.preds$se
# my.upper.limits <- my.preds$pred + 1.96*my.preds$se
# plot(as.numeric(data1[572:584]), type = "b", ylim = c(-350, 300))
# points(as.numeric(my.preds$pred), type = "b", col = "red")
# points(as.numeric(my.lower.limits), type = "l")
# points(as.numeric(my.upper.limits), type = "l")
# 
# 
# 
# fit3 <- auto.arima(data1t)
# sum((as.numeric(data1v) - predict(fit3, n.ahead = 59))^2)
# plot(as.numeric(data1v), type = "b")
# points(asdf$x, col = "red")
# points(as.numeric(predict(fit3, n.ahead = 59, se.fit = FALSE)), type = "b", col = "red")
# 
# par(mfrow=c(2,1))
# acf(fit1.7$residuals, type = c("correlation"), ylab="Sample ACF")
# acf(fit1.7$residuals, type = c("partial"), ylab="Sample PACF")
# par(mfrow=c(1,1)) 
# 
# qqnorm(fit1.7$residuals)
# qqnorm(fit1.4$residuals)

################################################################################
# proj2

par(mfrow = c(1,1), mar = c(4,4,3,1))
plot(data2t, type="b", pch = 20, main = "Scatterplot", ylab = "Data")
acf(data2t, type = c("correlation"), main = "ACF Plot")
acf(data2t, type = c("partial"), main = "PACF Plot")
par(mfrow = c(1,1), mar = c(4,4.5,3,1))

adf.test(data2t)

auto.arima(data2t)

spectrum(data2t, method = "ar")
spectrum(data2t, method = "pgram", ylab = expression(paste("I(", omega, ")")), main = "Periodogram")

tsdiag(fit2.1)

fit2.1 <- arima(data2t, order = c(5, 0, 5))
fit2.2 <- arima(data2t, order = c(5, 0, 3))
fit2.3 <- arima(data2t, order = c(5, 0, 2))
fit2.4 <- arima(data2t, order = c(5, 0, 1))
fit2.1
fit2.2
fit2.3
fit2.4
tsdiag(fit2.1)

fit2.1 <- arima(data2t, order = c(4, 0, 4)) #keep
fit2.2 <- arima(data2t, order = c(4, 0, 3))
fit2.3 <- arima(data2t, order = c(4, 0, 2))
fit2.4 <- arima(data2t, order = c(4, 0, 1))
fit2.1
fit2.4
tsdiag(fit2.1)

fit2.1 <- arima(data2t, order = c(3, 0, 4))
fit2.2 <- arima(data2t, order = c(3, 0, 3)) #keep
fit2.3 <- arima(data2t, order = c(3, 0, 2))
fit2.4 <- arima(data2t, order = c(3, 0, 1))
fit2.2
fit2.4
tsdiag(fit2.2)

fit2.1 <- arima(data2t, order = c(2, 0, 4))
fit2.2 <- arima(data2t, order = c(2, 0, 3))
fit2.3 <- arima(data2t, order = c(2, 0, 2)) #keep
fit2.4 <- arima(data2t, order = c(2, 0, 1))
fit2.3
fit2.4
tsdiag(fit2.3)

fit2.1 <- arima(data2t, order = c(4, 0, 4))
fit2.2 <- arima(data2t, order = c(3, 0, 3)) 
fit2.3 <- arima(data2t, order = c(2, 0, 2)) 
fit2.4 <- arima(data2t, order = c(5, 0, 5)) 
fit2.5 <- arima(data2t, order = c(6, 0, 6)) 

tsdiag(fit2.5)

sum((as.numeric(data2v) - predict(fit2.1, n.ahead = 20, se.fit = FALSE))^2)
sum((as.numeric(data2v) - predict(fit2.2, n.ahead = 20, se.fit = FALSE))^2)
sum((as.numeric(data2v) - predict(fit2.3, n.ahead = 20, se.fit = FALSE))^2)
sum((as.numeric(data2v) - predict(fit2.4, n.ahead = 20, se.fit = FALSE))^2)
sum((as.numeric(data2v) - predict(fit2.5, n.ahead = 20, se.fit = FALSE))^2)

plot(700:800, as.numeric(data2[700:800]), type = "b")
plot(770:800, as.numeric(data2[770:800]), type = "b")
points(781:800,as.numeric(predict(fit2.1, n.ahead = 20, se.fit = FALSE)), type = "b", col = "red")
points(781:800,as.numeric(predict(fit2.2, n.ahead = 20, se.fit = FALSE)), type = "b", col = "green")
points(781:800,as.numeric(predict(fit2.3, n.ahead = 20, se.fit = FALSE)), type = "b", col = "yellow")
points(781:800,as.numeric(predict(fit2.4, n.ahead = 20, se.fit = FALSE)), type = "b", col = "blue")
points(781:800,as.numeric(predict(fit2.5, n.ahead = 20, se.fit = FALSE)), type = "b", col = "pink")
points(781:800,as.numeric(predict(fit2.6, n.ahead = 20, se.fit = FALSE)), type = "b", col = "magenta")

par(mfrow=c(2,1))
spectrum(data2t,ylab = expression(paste("I(", omega, ")")), main = "Periodogram")
my.spectrum(phi.of.b=c(1, -fit2.1$coef[1:4]), theta.of.b=c(1, fit2.1$coef[5:8]), variance=fit2.1$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
spectrum(data2t,ylab = expression(paste("I(", omega, ")")), main = "Periodogram")
my.spectrum(phi.of.b=c(1, -fit2.2$coef[1:3]), theta.of.b=c(1, fit2.2$coef[4:7]), variance=fit2.2$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
spectrum(data2t,ylab = expression(paste("I(", omega, ")")), main = "Periodogram")
my.spectrum(phi.of.b=c(1, -fit2.3$coef[1:2]), theta.of.b=c(1, fit2.3$coef[3:4]), variance=fit2.3$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
spectrum(data2t,ylab = expression(paste("I(", omega, ")")), main = "Periodogram")
my.spectrum(phi.of.b=c(1, -fit2.4$coef[1:5]), theta.of.b=c(1, fit2.4$coef[6:10]), variance=fit2.4$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
spectrum(data2t,ylab = expression(paste("I(", omega, ")")), main = "Periodogram")
my.spectrum(phi.of.b=c(1, -fit2.5$coef[1:6]), theta.of.b=c(1, fit2.5$coef[7:12]), variance=fit2.5$sigma2)
par(mfrow=c(1,1))

plot(as.numeric(data2v), type = "b", pch = 20)
points(as.numeric(predict(fit2.1, n.ahead = 20, se.fit = FALSE)), type = "b", pch = 4)
points(as.numeric(predict(fit2.3, n.ahead = 20, se.fit = FALSE)), type = "b", pch = 1)

par(mfrow = c(2,2), mar = c(4,4,3,1))
acf(data2t, type = c("correlation"), main = "ACF Plot")
points(0:27, ARMAacf(ar = as.numeric(c(fit2.1$coef[1:4])), ma = as.numeric(c(fit2.1$coef[5:8])), lag.max=27, pacf=FALSE), pch = 4)
acf(data2t, type = c("partial"), main = "PACF Plot")
points(1:27, ARMAacf(ar = as.numeric(c(fit2.1$coef[1:4])), ma = as.numeric(c(fit2.1$coef[5:8])), lag.max=27, pacf=TRUE), pch = 4)
acf(data2t, type = c("correlation"), main = "ACF Plot")
points(0:27, ARMAacf(ar = as.numeric(c(fit2.3$coef[1:2])), ma = as.numeric(c(fit2.3$coef[3:4])), lag.max=27, pacf=FALSE), pch = 4)
acf(data2t, type = c("partial"), main = "PACF Plot")
points(1:27, ARMAacf(ar = as.numeric(c(fit2.3$coef[1:2])), ma = as.numeric(c(fit2.3$coef[3:4])), lag.max=27, pacf=TRUE), pch = 4)
par(mfrow = c(1,1), mar = c(4,4.5,3,1))

fit2.1 <- arima(data2, order = c(4, 0, 4))
fit2.3 <- arima(data2, order = c(2, 0, 2)) 

trueacf <- acf(data2, type = c("correlation"), plot = FALSE, lag.max = 20)
truepacf <- acf(data2, type = c("partial"), plot = FALSE, lag.max = 20)

acf1<- ARMAacf(ar = as.numeric(c(fit2.1$coef[1:4])), ma = as.numeric(c(fit2.1$coef[5:8])), lag.max=20, pacf=FALSE)
pacf1<-ARMAacf(ar = as.numeric(c(fit2.1$coef[1:4])), ma = as.numeric(c(fit2.1$coef[5:8])), lag.max=20, pacf=TRUE)

acf2<-ARMAacf(ar = as.numeric(c(fit2.3$coef[1:2])), ma = as.numeric(c(fit2.3$coef[3:4])), lag.max=20, pacf=FALSE)
pacf2<-ARMAacf(ar = as.numeric(c(fit2.3$coef[1:2])), ma = as.numeric(c(fit2.3$coef[3:4])), lag.max=20, pacf=TRUE)

acf3<-ARMAacf(ar = as.numeric(c(fit2.2$coef[1:3])), ma = as.numeric(c(fit2.2$coef[4:6])), lag.max=20, pacf=FALSE)
pacf3<-ARMAacf(ar = as.numeric(c(fit2.2$coef[1:3])), ma = as.numeric(c(fit2.2$coef[4:6])), lag.max=20, pacf=TRUE)



sum((as.numeric(acf1) - trueacf$acf)^2)*1000
sum((as.numeric(pacf1) - truepacf$acf)^2)*1000

sum((as.numeric(acf2) - trueacf$acf)^2)*1000
sum((as.numeric(pacf2) - truepacf$acf)^2)*1000

sum((as.numeric(acf3) - trueacf$acf)^2)*1000
sum((as.numeric(pacf3) - truepacf$acf)^2)*1000


final2 <- arima(data2, order = c(4, 0, 4))

Arg(polyroot(c(1, 0.2581, 0.4421, 0.0300, 0.823)))


##############ar1      ar2      ar3     ar4      ma1      ma2      ma3     ma4  intercept
##########-0.2581  -0.4421  -0.0300  -0.823  -0.2112  -0.1004  -0.3618  0.8804    99.7526
####s.e.   0.0227   0.0259   0.0259   0.022      NaN   0.0142      NaN     NaN     0.1624

initval <- runif(9, 0, 0.8)
arima(data2, order = c(4, 0, 4), init = initval)

final2 <- arima(data2, order = c(4, 0, 4))

final2
tsdiag(final2)
qqnorm(final2$residuals)

my.preds <- predict(final2, n.ahead = 13)
my.lower.limits  <- my.preds$pred - 1.96*my.preds$se
my.upper.limits <- my.preds$pred + 1.96*my.preds$se

par(mfrow = c(2,1), mar = c(3.9,4,0.5,1))
plot(1:800, data2[1:800], xlim = c(1, 813), type = "l", xlab = "Time", ylab = "Data")
points(801:813, as.numeric(my.preds$pred), type = "b", pch = 4)
points(801:813, as.numeric(my.lower.limits), type = "l", lty = 2)
points(801:813, as.numeric(my.upper.limits), type = "l", lty = 2)
plot(780:800, data2[780:800], xlim = c(780, 813), ylim = c(60,140), type = "l", xlab = "Time", ylab = "Data")
points(801:813, as.numeric(my.preds$pred), type = "b", pch = 4)
points(800:813, c(data2[800],as.numeric(my.lower.limits)), type = "l",lty = 2)
points(800:813, c(data2[800],as.numeric(my.upper.limits)), type = "l",lty = 2)
points(800:801, c(data2[800], my.preds$pred[1]), type = "l")
par(mfrow = c(1,1))

x<-arima.sim(model=list(ma=c(-0.2112, -0.1004, -0.3618, 0.8804), ar=c(-0.2581, -0.4421, -0.0300, -0.823)) , n=800 )
par(mfrow = c(2,1))
spectrum(x)
spectrum(data2)
par(mfrow = c(1,1))

x<-arima.sim(model=list(ma=c(-1.5202, 0.8835), ar=c(1.0464, -0.8144)) , n=800 )
par(mfrow = c(2,1))
spectrum(x)
spectrum(data2)
par(mfrow = c(1,1))

x<-arima.sim(model=list(ma=c(-0.5408, -0.5999, 0.8597 ), ar=c(0.0608,  0.2186, -0.8070)) , n=800 )
par(mfrow = c(2,1))
spectrum(x)
spectrum(data2)
par(mfrow = c(1,1))



################################################################################
# proj3

par(mfrow = c(3,1), mar = c(4,4,3,1))
plot(data3t, type="b", pch = 20, main = "Scatterplot", ylab = "Data")
acf(data3t, type = c("correlation"), main = "ACF Plot")
acf(data3t, type = c("partial"), main = "PACF Plot")
par(mfrow = c(1,1), mar = c(4,4.5,3,1))

auto.arima(data3)

adf.test(data3t)
BoxCox.lambda(data3 + 18034.8)

data3l <- ((data3+18034.8)^BoxCox.lambda(data3 + 18034.8)-1)/BoxCox.lambda(data3 + 18034.8)
data3lt <- data3l[1:480]
data3lv <- data3l[481:500]

par(mfrow = c(3,1), mar = c(4,4,3,1))
plot(data3lt, type="b", pch = 20, main = "Scatterplot", ylab = "Data")
acf(data3lt, type = c("correlation"), main = "ACF Plot")
acf(data3lt, type = c("partial"), main = "PACF Plot")
par(mfrow = c(1,1), mar = c(4,4.5,3,1))

data3ltd <- diff(data3lt, differences = 1)
par(mfrow = c(3,1), mar = c(4,4,3,1))
plot(data3ltd, type="b", pch = 20, main = "Scatterplot", ylab = "Data")
acf(data3ltd, type = c("correlation"), main = "ACF Plot")
acf(data3ltd, type = c("partial"), main = "PACF Plot")
par(mfrow = c(1,1), mar = c(4,4.5,3,1))

spectrum(data3ltd, method = "pgram", ylab = expression(paste("I(", omega, ")")), main = "Periodogram")

data3td1 <- diff(data3t, differences = 1)
par(mfrow = c(3,1), mar = c(4,4,3,1))
plot(data3td1, type="b", pch = 20, main = "Scatterplot", ylab = "Data")
acf(data3td1, type = c("correlation"), main = "ACF Plot")
acf(data3td1, type = c("partial"), main = "PACF Plot")
par(mfrow = c(1,1), mar = c(4,4.5,3,1))

spectrum(data3td1, method = "ar")
spectrum(data3td1, method = "pgram", ylab = expression(paste("I(", omega, ")")), main = "Periodogram")

adf.test(data3td1)

# test ars
fit3.1 <- arima(data3t, order = c(7, 1, 0), xreg = 1:length(data3t)) #keep

# test ars with ma 2
fit3.2 <- arima(data3t, order = c(2, 1, 2), xreg = 1:length(data3t))
fit3.3 <- arima(data3t, order = c(4, 1, 2), xreg = 1:length(data3t))
fit3.4 <- arima(data3t, order = c(5, 1, 2), xreg = 1:length(data3t))

# test ars with ma 4
fit3.9 <- arima(data3t, order = c(1, 1, 4), xreg = 1:length(data3t))

par(mfrow=c(2,1))
spectrum(data3td1,ylab = expression(paste("I(", omega, ")")), main = "Periodogram")
my.spectrum(phi.of.b=c(1, -fit3.1$coef[1:7]), theta.of.b=c(1), variance=fit3.1$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
spectrum(data3td1,ylab = expression(paste("I(", omega, ")")), main = "Periodogram")
my.spectrum(phi.of.b=c(1, -fit3.2$coef[1:2]), theta.of.b=c(1, fit3.2$coef[3:4]), variance=fit3.2$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
spectrum(data3td1,ylab = expression(paste("I(", omega, ")")), main = "Periodogram")
my.spectrum(phi.of.b=c(1, -fit3.3$coef[1:4]), theta.of.b=c(1, fit3.3$coef[5:6]), variance=fit3.3$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
spectrum(data3td1,ylab = expression(paste("I(", omega, ")")), main = "Periodogram")
my.spectrum(phi.of.b=c(1, -fit3.4$coef[1:5]), theta.of.b=c(1, fit3.4$coef[6:7]), variance=fit3.4$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
spectrum(data3td1,ylab = expression(paste("I(", omega, ")")), main = "Periodogram")
my.spectrum(phi.of.b=c(1, -fit3.9$coef[1]), theta.of.b=c(1, fit3.9$coef[2:5]), variance=fit3.9$sigma2)
par(mfrow=c(1,1))

# test ars wiht ma 5


# final candidates

data3t <- data3[1:450]
data3v <- data3[451:500]

fit3.1 <- arima(data3t, order = c(7, 1, 0), xreg = 1:length(data3t))
fit3.2 <- arima(data3t, order = c(2, 1, 2), xreg = 1:length(data3t))
fit3.3 <- arima(data3t, order = c(4, 1, 2), xreg = 1:length(data3t))
fit3.4 <- arima(data3t, order = c(5, 1, 2), xreg = 1:length(data3t))
fit3.5 <- arima(data3t, order = c(1, 1, 4), xreg = 1:length(data3t))

sum((as.numeric(data3v) - predict(fit3.1, n.ahead = 50, newxreg = ((length(data3t)+1):(length(data3t)+50)), se.fit = FALSE))^2)
sum((as.numeric(data3v) - predict(fit3.2, n.ahead = 50, newxreg = ((length(data3t)+1):(length(data3t)+50)), se.fit = FALSE))^2)
sum((as.numeric(data3v) - predict(fit3.3, n.ahead = 50, newxreg = ((length(data3t)+1):(length(data3t)+50)), se.fit = FALSE))^2)
sum((as.numeric(data3v) - predict(fit3.4, n.ahead = 50, newxreg = ((length(data3t)+1):(length(data3t)+50)), se.fit = FALSE))^2)
sum((as.numeric(data3v) - predict(fit3.5, n.ahead = 50, newxreg = ((length(data3t)+1):(length(data3t)+50)), se.fit = FALSE))^2)

my.preds <- predict(fit3.2, n.ahead = 20, newxreg = ((length(data3t)+1):(length(data3t)+20)), se.fit = TRUE)
plot(450:500, data3[450:500], xlim=c(450,500), ylim = c(-19000, -16000), type="b")
lines(481:500, my.preds$pred, type="b", col=2)
lines(481:500, my.preds$pred + 2*my.preds$se, type="l", col=2)
lines(481:500, my.preds$pred - 2*my.preds$se, type="l", col=2)

data3t <- data3[1:490]
data3v <- data3[491:500]

fit3.1 <- arima(data3t, order = c(7, 1, 0), xreg = 1:length(data3t))
fit3.2 <- arima(data3t, order = c(2, 1, 2), xreg = 1:length(data3t))
fit3.3 <- arima(data3t, order = c(4, 1, 2), xreg = 1:length(data3t))
fit3.4 <- arima(data3t, order = c(5, 1, 2), xreg = 1:length(data3t))
fit3.5 <- arima(data3t, order = c(1, 1, 4), xreg = 1:length(data3t))

sum((as.numeric(data3v) - predict(fit3.1, n.ahead = 10, newxreg = ((length(data3t)+1):(length(data3t)+10)), se.fit = FALSE))^2)
sum((as.numeric(data3v) - predict(fit3.2, n.ahead = 10, newxreg = ((length(data3t)+1):(length(data3t)+10)), se.fit = FALSE))^2)
sum((as.numeric(data3v) - predict(fit3.3, n.ahead = 10, newxreg = ((length(data3t)+1):(length(data3t)+10)), se.fit = FALSE))^2)
sum((as.numeric(data3v) - predict(fit3.4, n.ahead = 10, newxreg = ((length(data3t)+1):(length(data3t)+10)), se.fit = FALSE))^2)
sum((as.numeric(data3v) - predict(fit3.5, n.ahead = 10, newxreg = ((length(data3t)+1):(length(data3t)+10)), se.fit = FALSE))^2)

my.preds <- predict(fit3.2, n.ahead = 10, newxreg = ((length(data3t)+1):(length(data3t)+10)), se.fit = TRUE)
plot(450:500, data3[450:500], xlim=c(450,500), ylim = c(-19000, -16000), type="b")
lines(491:500, my.preds$pred, type="b", col=2)
lines(491:500, my.preds$pred + 2*my.preds$se, type="l", col=2)
lines(491:500, my.preds$pred - 2*my.preds$se, type="l", col=2)

fit3.1 <- arima(data3, order = c(7, 1, 0), xreg = 1:length(data3))
fit3.2 <- arima(data3, order = c(2, 1, 2), xreg = 1:length(data3))
fit3.3 <- arima(data3, order = c(4, 1, 2), xreg = 1:length(data3))
fit3.4 <- arima(data3, order = c(5, 1, 2), xreg = 1:length(data3))
fit3.5 <- arima(data3, order = c(1, 1, 4), xreg = 1:length(data3))
tsdiag(fit3.1)
tsdiag(fit3.2)
tsdiag(fit3.3)
tsdiag(fit3.4)
tsdiag(fit3.5)

data3d <- diff(data3, differences = 1)

par(mfrow=c(2,1))
spectrum(data3d,ylab = expression(paste("I(", omega, ")")), main = "Periodogram")
my.spectrum(phi.of.b=c(1, -fit3.1$coef[1:7]), theta.of.b=c(1), variance=fit3.1$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
spectrum(data3d,ylab = expression(paste("I(", omega, ")")), main = "Periodogram")
my.spectrum(phi.of.b=c(1, -fit3.2$coef[1:2]), theta.of.b=c(1, fit3.2$coef[3:4]), variance=fit3.2$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
spectrum(data3d,ylab = expression(paste("I(", omega, ")")), main = "Periodogram")
my.spectrum(phi.of.b=c(1, -fit3.3$coef[1:4]), theta.of.b=c(1, fit3.3$coef[5:6]), variance=fit3.3$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
spectrum(data3d,ylab = expression(paste("I(", omega, ")")), main = "Periodogram")
my.spectrum(phi.of.b=c(1, -fit3.4$coef[1:5]), theta.of.b=c(1, fit3.4$coef[6:7]), variance=fit3.4$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
spectrum(data3d,ylab = expression(paste("I(", omega, ")")), main = "Periodogram")
my.spectrum(phi.of.b=c(1, -fit3.5$coef[1]), theta.of.b=c(1, fit3.5$coef[2:5]), variance=fit3.5$sigma2)
par(mfrow=c(1,1))

ARMAacf(ar = as.numeric(c(fit3.1$coef[1:7])), ma = numeric())

plot(0:27, ARMAacf(ar = as.numeric(c(fit3.1$coef[1:7])), ma = numeric(), lag.max=27, pacf=FALSE), type="h")
abline(h=0)

plot(1:27, ARMAacf(ar = as.numeric(c(fit3.1$coef[1:7])), ma = numeric(), lag.max=27, pacf=TRUE), type="h")
abline(h=0)

par(mfrow = c(2,1), mar = c(4,4,3,1))
acf(data3d, type = c("correlation"), main = "ACF Plot")
points(0:27, ARMAacf(ar = as.numeric(c(fit3.1$coef[1:7])), ma = numeric(), lag.max=27, pacf=FALSE), pch = 4)
acf(data3d, type = c("partial"), main = "PACF Plot")
points(1:27, ARMAacf(ar = as.numeric(c(fit3.1$coef[1:7])), ma = numeric(), lag.max=27, pacf=TRUE), pch = 4)
par(mfrow = c(1,1), mar = c(4,4.5,3,1))

par(mfrow = c(2,1), mar = c(4,4,3,1))
acf(data3d, type = c("correlation"), main = "ACF Plot")
points(0:27, ARMAacf(ar = as.numeric(c(fit3.2$coef[1:2])), ma = as.numeric(c(fit3.2$coef[3:4])), lag.max=27, pacf=FALSE), pch = 4)
acf(data3d, type = c("partial"), main = "PACF Plot")
points(1:27, ARMAacf(ar = as.numeric(c(fit3.2$coef[1:2])), ma = as.numeric(c(fit3.2$coef[3:4])), lag.max=27, pacf=TRUE), pch = 4)
par(mfrow = c(1,1), mar = c(4,4.5,3,1))

par(mfrow = c(2,1), mar = c(4,4,3,1))
acf(data3d, type = c("correlation"), main = "ACF Plot")
points(0:27, ARMAacf(ar = as.numeric(c(fit3.3$coef[1:4])), ma = as.numeric(c(fit3.3$coef[5:6])), lag.max=27, pacf=FALSE), pch = 4)
acf(data3d, type = c("partial"), main = "PACF Plot")
points(1:27, ARMAacf(ar = as.numeric(c(fit3.3$coef[1:4])), ma = as.numeric(c(fit3.3$coef[5:6])), lag.max=27, pacf=TRUE), pch = 4)
par(mfrow = c(1,1), mar = c(4,4.5,3,1))

par(mfrow = c(2,1), mar = c(4,4,3,1))
acf(data3d, type = c("correlation"), main = "ACF Plot")
points(0:27, ARMAacf(ar = as.numeric(c(fit3.4$coef[1:5])), ma = as.numeric(c(fit3.4$coef[6:7])), lag.max=27, pacf=FALSE), pch = 4)
acf(data3d, type = c("partial"), main = "PACF Plot")
points(1:27, ARMAacf(ar = as.numeric(c(fit3.4$coef[1:5])), ma = as.numeric(c(fit3.4$coef[6:7])), lag.max=27, pacf=TRUE), pch = 4)
par(mfrow = c(1,1), mar = c(4,4.5,3,1))

par(mfrow = c(2,1), mar = c(4,4,3,1))
acf(data3d, type = c("correlation"), main = "ACF Plot")
points(0:27, ARMAacf(ar = as.numeric(c(fit3.5$coef[1])), ma = as.numeric(c(fit3.5$coef[2:5])), lag.max=27, pacf=FALSE), pch = 4)
acf(data3d, type = c("partial"), main = "PACF Plot")
points(1:27, ARMAacf(ar = as.numeric(c(fit3.5$coef[1])), ma = as.numeric(c(fit3.5$coef[2:5])), lag.max=27, pacf=TRUE), pch = 4)
par(mfrow = c(1,1), mar = c(4,4.5,3,1))

fit3.1 <- arima(data3, order = c(7, 1, 0), xreg = 1:length(data3))
fit3.2 <- arima(data3, order = c(2, 1, 2), xreg = 1:length(data3))
fit3.3 <- arima(data3, order = c(4, 1, 2), xreg = 1:length(data3))
fit3.4 <- arima(data3, order = c(5, 1, 2), xreg = 1:length(data3))
fit3.5 <- arima(data3, order = c(1, 1, 4), xreg = 1:length(data3))

fit3.1
summary(fit3.1)


data3d <- diff(data3, differences = 1)

trueacf <- acf(data3d, type = c("correlation"), plot = FALSE, lag.max = 20)
truepacf <- acf(data3d, type = c("partial"), plot = FALSE, lag.max = 7)

acf1<- ARMAacf(ar = as.numeric(c(fit3.1$coef[1:7])), ma = as.numeric(), lag.max=20, pacf=FALSE)
pacf1<-ARMAacf(ar = as.numeric(c(fit3.1$coef[1:7])), ma = as.numeric(), lag.max=7, pacf=TRUE)
acf2<-ARMAacf(ar = as.numeric(c(fit3.2$coef[1:2])), ma = as.numeric(c(fit3.2$coef[3:4])), lag.max=20, pacf=FALSE)
pacf2<-ARMAacf(ar = as.numeric(c(fit3.2$coef[1:2])), ma = as.numeric(c(fit3.2$coef[3:4])), lag.max=7, pacf=TRUE)
acf3<-ARMAacf(ar = as.numeric(c(fit3.3$coef[1:4])), ma = as.numeric(c(fit3.3$coef[5:6])), lag.max=20, pacf=FALSE)
pacf3<-ARMAacf(ar = as.numeric(c(fit3.3$coef[1:4])), ma = as.numeric(c(fit3.3$coef[5:6])), lag.max=7, pacf=TRUE)
acf4<-ARMAacf(ar = as.numeric(c(fit3.4$coef[1:5])), ma = as.numeric(c(fit3.4$coef[6:7])), lag.max=20, pacf=FALSE)
pacf4<-ARMAacf(ar = as.numeric(c(fit3.4$coef[1:5])), ma = as.numeric(c(fit3.4$coef[6:7])), lag.max=7, pacf=TRUE)
acf5<-ARMAacf(ar = as.numeric(c(fit3.5$coef[1])), ma = as.numeric(c(fit3.5$coef[2:5])), lag.max=20, pacf=FALSE)
pacf5<-ARMAacf(ar = as.numeric(c(fit3.5$coef[1])), ma = as.numeric(c(fit3.5$coef[2:5])), lag.max=7, pacf=TRUE)

sum((as.numeric(acf1) - trueacf$acf)^2)*1000
sum((as.numeric(pacf1) - truepacf$acf)^2)*1000

sum((as.numeric(acf2) - trueacf$acf)^2)*1000
sum((as.numeric(pacf2) - truepacf$acf)^2)*1000

sum((as.numeric(acf3) - trueacf$acf)^2)*1000
sum((as.numeric(pacf3) - truepacf$acf)^2)*1000

sum((as.numeric(acf4) - trueacf$acf)^2)*1000
sum((as.numeric(pacf4) - truepacf$acf)^2)*1000

sum((as.numeric(acf5) - trueacf$acf)^2)*1000
sum((as.numeric(pacf5) - truepacf$acf)^2)*1000

fit3.1
fit3.2
fit3.3
fit3.5

tsdiag(fit3.1)
tsdiag(fit3.2)
tsdiag(fit3.3)
tsdiag(fit3.5)

fit3.5 <- arima(data3, order = c(1, 1, 4), xreg = 1:length(data3))

BIC(arima(data3, order = c(1, 1, 4), xreg = 1:length(data3)))
BIC(arima(data3, order = c(1, 1, 4)))

arima(data3, order = c(1, 1, 4), xreg = 1:length(data3))
final3 <- arima(data3, order = c(1, 1, 4))
tsdiag(final3)
qqnorm(final3$residuals)


my.preds <- predict(final3, n.ahead = 13)
my.lower.limits  <- my.preds$pred - 1.96*my.preds$se
my.upper.limits <- my.preds$pred + 1.96*my.preds$se

par(mfrow = c(2,1), mar = c(3.9,4,0.5,1))
plot(1:500, data3[1:500], xlim = c(1, 513), type = "l", xlab = "Time", ylab = "Data")
points(501:513, as.numeric(my.preds$pred), type = "b", pch = 4)
points(501:513, as.numeric(my.lower.limits), type = "l", lty = 2)
points(501:513, as.numeric(my.upper.limits), type = "l", lty = 2)
plot(480:500, data3[480:500], xlim = c(480, 513), ylim = c(-19000, -15000), type = "l", xlab = "Time", ylab = "Data")
points(501:513, as.numeric(my.preds$pred), type = "b", pch = 4)
points(500:513, c(data3[500],as.numeric(my.lower.limits)), type = "l",lty = 2)
points(500:513, c(data3[500],as.numeric(my.upper.limits)), type = "l",lty = 2)
points(500:501, c(data3[500], my.preds$pred[1]), type = "l")
par(mfrow = c(1,1))

my.upper.limits/1000
my.preds$pred/1000
my.lower.limits/1000

################################################################################
# proj4 deposits.txt

par(mfrow = c(3,1), mar = c(4,4,3,1))
plot(data4, type="b", pch = 20, main = "Scatterplot", ylab = "Data")
acf(data4, type = c("correlation"), main = "ACF Plot")
acf(data4, type = c("partial"), main = "PACF Plot")
par(mfrow = c(1,1), mar = c(4,4.5,3,1))

BoxCox.lambda(data4)

data4l <- log(data4)

par(mfrow = c(3,1), mar = c(4,4,3,1))
plot(data4l, type="b", pch = 20, main = "Scatterplot", ylab = "Data")
acf(data4l, type = c("correlation"), main = "ACF Plot")
acf(data4l, type = c("partial"), main = "PACF Plot")
par(mfrow = c(1,1), mar = c(4,4.5,3,1))

auto.arima(data4l, )

data4ld <- diff(data4l, differences = 1)

par(mfrow = c(3,1), mar = c(4,4,3,1))
plot(data4ld, type="b", pch = 20, main = "Scatterplot", ylab = "Data")
acf(data4ld, type = c("correlation"), main = "ACF Plot")
acf(data4ld, type = c("partial"), main = "PACF Plot")
par(mfrow = c(1,1), mar = c(4,4.5,3,1))

adf.test(data4ld)
spectrum(data4ld, method = "pgram", ylab = expression(paste("I(", omega, ")")), main = "Periodogram")

data4lt <- data4l[1:604]
data4lv <- data4l[605:624]

fit4.1 <- arima(data4lt, order = c(0, 1, 2), xreg = 1:length(data4lt))
fit4.2 <- arima(data4lt, order = c(1, 1, 1), xreg = 1:length(data4lt))

fit4.1
fit4.2

sum((as.numeric(data4lv) - predict(fit4.1, n.ahead = 20, newxreg = ((length(data4lt)+1):(length(data4lt)+20)), se.fit = FALSE))^2)
sum((as.numeric(data4lv) - predict(fit4.2, n.ahead = 20, newxreg = ((length(data4lt)+1):(length(data4lt)+20)), se.fit = FALSE))^2)

arima(data4l, order = c(0, 1, 2), xreg = 1:length(data4l))
arima(data4l, order = c(1, 1, 1), xreg = 1:length(data4l))

par(mfrow=c(2,1))
spectrum(data4ltd,ylab = expression(paste("I(", omega, ")")), main = "Periodogram")
my.spectrum(phi.of.b=c(1), theta.of.b=c(1, fit4.1$coef[1:2]), variance=fit4.1$sigma2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
spectrum(data4ltd,ylab = expression(paste("I(", omega, ")")), main = "Periodogram")
my.spectrum(phi.of.b=c(1, -fit4.2$coef[1]), theta.of.b=c(1, fit4.2$coef[2]), variance=fit4.2$sigma2)
par(mfrow=c(1,1))

arima(data4lt, order = c(1, 1, 2), xreg = 1:length(data4lt))

fit4.1 <- arima(data4lt, order = c(0, 1, 2), xreg = 1:length(data4lt))
fit4.2 <- arima(data4lt, order = c(1, 1, 1), xreg = 1:length(data4lt))
fit4.3 <- arima(data4lt, order = c(5, 1, 5), xreg = 1:length(data4lt))
fit4.4 <- arima(data4lt, order = c(4, 1, 4), xreg = 1:length(data4lt))
fit4.5 <- arima(data4lt, order = c(3, 1, 3), xreg = 1:length(data4lt))

sum((as.numeric(data4lv) - predict(fit4.5, n.ahead = 20, newxreg = ((length(data4lt)+1):(length(data4lt)+20)), se.fit = FALSE))^2)
sum((as.numeric(data4lv) - predict(fit4.4, n.ahead = 20, newxreg = ((length(data4lt)+1):(length(data4lt)+20)), se.fit = FALSE))^2)
sum((as.numeric(data4lv) - predict(fit4.3, n.ahead = 20, newxreg = ((length(data4lt)+1):(length(data4lt)+20)), se.fit = FALSE))^2)

tsdiag(fit4.4)

data4ltd <- diff(data4lt, differences = 1)

par(mfrow = c(2,2), mar = c(4,4,3,1))
acf(data4ltd, type = c("correlation"), main = "ACF Plot")
points(0:27, ARMAacf(ar = as.numeric(), ma = as.numeric(c(fit4.1$coef[1:2])), lag.max=27, pacf=FALSE), pch = 4)
acf(data4ltd, type = c("partial"), main = "PACF Plot")
points(1:27, ARMAacf(ar = as.numeric(), ma = as.numeric(c(fit4.1$coef[1:2])), lag.max=27, pacf=TRUE), pch = 4)
acf(data4ltd, type = c("correlation"), main = "ACF Plot")
points(0:27, ARMAacf(ar = as.numeric(c(fit4.2$coef[1])), ma = as.numeric(c(fit4.2$coef[2])), lag.max=27, pacf=FALSE), pch = 4)
acf(data4ltd, type = c("partial"), main = "PACF Plot")
points(1:27, ARMAacf(ar = as.numeric(c(fit4.2$coef[1])), ma = as.numeric(c(fit4.2$coef[2])), lag.max=27, pacf=TRUE), pch = 4)
par(mfrow = c(1,1), mar = c(4,4.5,3,1))

data4lt <- data4l[1:614]
data4lv <- data4l[615:624]
fit4.1 <- arima(data4lt, order = c(0, 1, 2), xreg = 1:length(data4lt))
fit4.2 <- arima(data4lt, order = c(1, 1, 1), xreg = 1:length(data4lt))
sum((as.numeric(data4lv) - predict(fit4.1, n.ahead = 10, newxreg = ((length(data4lt)+1):(length(data4lt)+10)), se.fit = FALSE))^2)
sum((as.numeric(data4lv) - predict(fit4.2, n.ahead = 10, newxreg = ((length(data4lt)+1):(length(data4lt)+10)), se.fit = FALSE))^2)

data4lt <- data4l[1:574]
data4lv <- data4l[575:624]
fit4.1 <- arima(data4lt, order = c(0, 1, 2), xreg = 1:length(data4lt))
fit4.2 <- arima(data4lt, order = c(1, 1, 1), xreg = 1:length(data4lt))
sum((as.numeric(data4lv) - predict(fit4.1, n.ahead = 50, newxreg = ((length(data4lt)+1):(length(data4lt)+50)), se.fit = FALSE))^2)
sum((as.numeric(data4lv) - predict(fit4.2, n.ahead = 50, newxreg = ((length(data4lt)+1):(length(data4lt)+50)), se.fit = FALSE))^2)

data4lt <- data4l[1:619]
data4lv <- data4l[620:624]
fit4.1 <- arima(data4lt, order = c(0, 1, 2))
fit4.2 <- arima(data4lt, order = c(1, 1, 1))
sum((as.numeric(data4lv) - predict(fit4.1, n.ahead = 5, se.fit = FALSE))^2)
sum((as.numeric(data4lv) - predict(fit4.2, n.ahead = 5, se.fit = FALSE))^2)

final4 <- arima(data4l, order = c(1, 1, 1), xreg = 1:length(data4l))
final4

final4 <- arima(data4l, order = c(1, 1, 1))
tsdiag(final4)

my.preds <- predict(final4, n.ahead = 13)
my.lower.limits  <- exp(my.preds$pred - 1.96*my.preds$se)
my.upper.limits <- exp(my.preds$pred + 1.96*my.preds$se)
my.preds <- exp(my.preds$pred)

par(mfrow = c(2,1), mar = c(3.9,4,0.5,1))
plot(1:624, data4[1:624], xlim = c(1, 637), type = "l", xlab = "Year", ylab = "Data")
points(625:637, as.numeric(my.preds), type = "b", pch = 4)
points(625:637, as.numeric(my.lower.limits), type = "l", lty = 2)
points(625:637, as.numeric(my.upper.limits), type = "l", lty = 2)
plot(604:624, data4[604:624], xlim = c(604, 637), ylim = c(0,35), type = "l", xlab = "Year", ylab = "Data")
points(625:637, as.numeric(my.preds), type = "b", pch = 4)
points(624:637, c(data4[624],as.numeric(my.lower.limits)), type = "l",lty = 2)
points(624:637, c(data4[624],as.numeric(my.upper.limits)), type = "l",lty = 2)
par(mfrow = c(1,1))

my.upper.limits
my.lower.limits
my.preds
qqnorm(final4$residuals)

