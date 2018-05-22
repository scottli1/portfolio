library("tseries") 
library("timsac") # Akaike's final prediction error for AR
library("forecast")
library("MASS")
library('fGarch')
library('xts')
library('roll')
library(astsa)
library(zoo)
library(zyp)

autofit <- auto.arima(data, trace = TRUE, stepwise = FALSE)
tsdiag(autofit)

data <- read.csv("HomePrice.csv", header = TRUE)
ddate <- as.Date(data$date, "%m/%d/%Y")
data <- ts(data$index, start = c(1987,1), frequency = 12)

plot(data, type = "l")
plot(decompose(data, type = "a"))
plot(decompose(data, type = "m"))
decomp <- decompose(data, type = "m")
plot(decomp$trend)
plot(stldecomp)
plot(stldecomp$time.series[,2])

diffdata <- diff(data, differences = 1)
plot(diffdata)

par(mfrow = c(3,1), mar = c(4,4,3,1))
plot(data, type="b", pch = 20, main = "Scatterplot", ylab = "Data")
acf(data, type = c("correlation"), main = "ACF Plot")
acf(data, type = c("partial"), main = "PACF Plot")
par(mfrow = c(1,1), mar = c(4,4.5,3,1))

par(mfrow = c(3,1), mar = c(4,4,3,1))
plot(diffdata, type="b", pch = 20, main = "Scatterplot", ylab = "Data")
acf(diffdata, type = c("correlation"), main = "ACF Plot")
acf(diffdata, type = c("partial"), main = "PACF Plot")
par(mfrow = c(1,1), mar = c(4,4.5,3,1))

sdiffdata <- diff(diffdata, lag = 12, differences = 1)
par(mfrow = c(3,1), mar = c(4,4,3,1))
plot(sdiffdata, type="b", pch = 20, main = "Scatterplot", ylab = "Data")
acf(sdiffdata, type = c("correlation"), main = "ACF Plot")
acf(sdiffdata, type = c("partial"), main = "PACF Plot")
par(mfrow = c(1,1), mar = c(4,4.5,3,1))

adf.test(sdiffdata)

dsdiffdata <- diff(sdiffdata, differences = 1)
par(mfrow = c(3,1), mar = c(4,4,3,1))
plot(dsdiffdata, type="b", pch = 20, main = "Scatterplot", ylab = "Data")
acf(dsdiffdata, type = c("correlation"), main = "ACF Plot")
acf(dsdiffdata, type = c("partial"), main = "PACF Plot")
par(mfrow = c(1,1), mar = c(4,4.5,3,1))

fit1 <- auto.arima(data)
fit2 <- auto.arima(sdiffdata)
plot(simulate(fit1))
plot(simulate(fit2))
arima.sim(model=list(ma=c(-0.5408, -0.5999, 0.8597 ), ar=c(0.8610 , 0.1035,  -0.2207)) , n=800 )

fit <- garchFit(~garch(1, 1), data = dsdiffdata, trace = FALSE)
plot(fit)
plot(volatility(fit))
plot(residuals(fit))

par(mfrow = c(5,1), mar = c(4,4,3,1))
plot(data, type="b", pch = 20, main = "Raw Data", ylab = "Data")
plot(diffdata, type="b", pch = 20, main = "Differenced Data", ylab = "Data")
plot(sdiffdata, type="b", pch = 20, main = "Seasonally Differenced Differenced Data", ylab = "Data")
plot(dsdiffdata, type="b", pch = 20, main = "Three Times Differenced", ylab = "Data")
plot(volatility(fit))
par(mfrow = c(1,1), mar = c(4,4.5,3,1))

#roll_lm()
#plot(apply.monthly(window(data, start=1) , function(x) sd(x)), ylim = c(0,4), type = "b")

sarima(data, 4,1,0,2,2,0,12)
par(mfrow = c(2,1))
acf(dsdiffdata, type = c("correlation"), main = "ACF Plot")
acf(dsdiffdata, data, type = c("partial"), main = "PACF Plot")
par(mfrow = c(1,1))

plot(diffdata)
plot(ddate, rollapply(data, width = 6, FUN = var, fill = NA), xlab = "Year", ylab = "Variance", type = "l", main = "6 Month Rolling Variance of the Home Index")
abline(v = ddate[c(175,192,217,236,332)], lty = 2)

plot(ddate, rollapply(decomp$trend, width = 6, FUN = var, fill = NA), xlab = "Year", ylab = "Variance", type = "l", main = "6 Month Rolling Variance of the Home Index")



plot(rollapply(data, width = 2, FUN = var, fill = NA))
plot(rollapply(data, width = 13, FUN = var, fill = NA))
plot(rollapply(data, width = 20, FUN = sd, fill = NA))

library(pspline)
predict(sm.spline(1:length(data), data), 1:length(data), 1)
plot(predict(sm.spline(1:length(data), data), 1:length(data), 1), type = "l")

trendslope <- diff(decomp$trend, differences = 1)
par(mfrow = c(2,1), mar = c(4,4,3,1))
plot(ddate, data, xlab = "Year", ylab = "Index", main = "Home Index", type = "l")
text(ddate[c(175,192,217,236,332)], data[c(175,192,217,236,332)],
     labels = ddate[c(175,192,217,236,332)], cex= 1 , pos = 1, offset =0.25)
plot(ddate[-1], trendslope, xlab = "Year", ylab = "Slope", main = "Estimated Slope of Trend", type = "l")
text(ddate[c(175,192,217,236,332)], trendslope[c(175,192,217,236,332)],
     labels = ddate[c(175,192,217,236,332)], cex= 1 , pos = 1, offset =0.25)
abline(0,0, lty = 2)

plot(trendslope)
abline(0,0, lty = 2)
identify(trendslope)
abline(0,0, lty = 2)

library(prophet)
datadf <- as.data.frame(data)
colnames(datadf) <- "y"
datadf$ds <- ddate
pmodel <- prophet(datadf[1:230,], n.changepoints = 1, changepoint.prior.scale = 0.01)
pmodel <- prophet(datadf[200:270,], n.changepoints = 2, changepoint.prior.scale = 0.01)
pmodel <- prophet(datadf[270:368,], n.changepoints = 3, changepoint.prior.scale = 0.01)


plot(data)
#abline(v = (as.POSIXct(gsub(" GMT", "", pmodel$changepoints))))

prebubble <- ts(data[1:181], c(1987,1), c(2002,1), frequency = 12)

prebmodel <- auto.arima(prebubble, stepwise = FALSE, trace = TRUE)
prebmodel
tsdiag(prebmodel)
qqnorm((prebmodel$residuals))


plot(forecast(prebubble, h = 187, level = 80), ylim = c(60, 300), ylab = "Index", xlab = "Date", main = "Pre-Bubble SARIMA Forecast", lwd = 3)
points(data, type = "l", lty = 3, lwd = 3)


fit1 <-hw(prebubble, seasonal="multiplicative")
fit2 <- hw(prebubble, seasonal="additive")
hwprebmodel <- hw(prebubble, seasonal="multiplicative", h = 187, level = 80)
plot(hwprebmodel, xlab = "Year", ylab = "Index", lwd = 3, ylim = c(60, 300), main = "Holt-Winters Forecast")
points(data, type = "l", lty = 3, lwd = 3)

hwprebmodel <- hw(prebubble, seasonal="multiplicative", h = 187, level = 80)
plot(hwprebmodel, xlab = "Year", ylab = "Index", lwd = 3, ylim = c(60, 300))
points(data, type = "l", lty = 3, lwd = 3)

plot(rwf(prebubble, h = 187, drift=TRUE, level = 80), xlab = "Year", ylab = "Index", lwd = 3, ylim = c(60, 300))
points(data, type = "l", lty = 3, lwd = 3)

fullmodel <- auto.arima(data, stepwise = FALSE, trace = TRUE)
fullmodel
tsdiag(fullmodel)

plot(forecast(fullmodel, h = 60, level = 80), ylim = c(60, 300), ylab = "Index", xlab = "Date", main = "SARIMA Forecast for Present and Future", lwd = 3)

auto.arima(data, allowmean = FALSE)

sdata <- diff(data, lag = 12, differences = 1)
auto.arima(sdata)

prebmodel <- auto.arima(prebubble, stepwise = FALSE, trace = TRUE, lambda = 0)
prebmodel
tsdiag(prebmodel)
qqnorm((prebmodel$residuals))

plot(forecast(prebmodel, h = 187, level = 80, lambda = 0), ylim = c(60, 300), ylab = "Index", xlab = "Date", main = "Pre-Bubble SARIMA Forecast", lwd = 3)
points(data, type = "l", lty = 3, lwd = 3)

plot(forecast(prebmodel, h = 187, level = 80, lambda = 0), ylim = c(60, 300), ylab = "Index", xlab = "Date", main = "Pre-Bubble SARIMA Forecast", lwd = 3)
points(data, type = "l", lty = 3, lwd = 3)

plot(data[1:50])
plot(data[131:181])

ddate <- ddate - 2

par(mar = c(4,4,2,2))
plot(ddate, data, xlab = "Year", ylab = "Index", main = "Case-Shiller Index", type = "l", lwd = 2)
#text(ddate[c(182,225,234)], data[c(182,225,234)],
#     labels = ddate[c(182,225,234)], cex= 1 , pos = 1, offset =0.25)
abline(v = ddate[c(182,225,234)], lty = 2)


postbubble <- ts(data[301:368], c(2012,1), c(2017,8), frequency = 12)

hwprebmodel <- hw(data, seasonal="multiplicative", h = 60, level = 80)
plot(hwprebmodel, xlab = "Year", ylab = "Index", lwd = 3, ylim = c(60, 300))
plot(forecast(fullmodel, h = 60, level = 80), ylim = c(60, 300), ylab = "Index", xlab = "Date", main = "SARIMA Forecast for Present and Future", lwd = 3)
qqnorm(fullmodel$residuals)


