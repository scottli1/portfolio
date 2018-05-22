# data <- left_join(sp500, claims, by = "DATE")
# data <- left_join(data, corecpi, by = "DATE")
# data <- left_join(data, credit, by = "DATE")
# data <- left_join(data, disposincome, by = "DATE")
# data <- left_join(data, durable, by = "DATE")
# data <- left_join(data, fedrate, by = "DATE")
# data <- left_join(data, housestarts, by = "DATE")
# data <- left_join(data, IPI, by = "DATE")
# data <- left_join(data, newhome, by = "DATE")
# data <- left_join(data, neworders, by = "DATE")
# data <- left_join(data, payroll, by = "DATE")
# data <- left_join(data, perincome, by = "DATE")
# data <- left_join(data, PMI, by = "DATE")
# data <- left_join(data, retail, by = "DATE")
# data <- left_join(data, unemployment, by = "DATE")
# names(data) <- c("Date", "Open", "High", "Low", "Close", "Volume", "Adjclose",
#                  "claims", "corecpi", "credit", "disposincome", "durable", 
#                  "fedrate", "housestarts", "IPI", "newhome", "neworders", "payroll",
#                  "perincome", "PMI", "retail", "unemployment")
# 
# # 327, 212 to 190, 90 to 70
# normal <- rep(1, length(data[,1]))
# normal[327] <- 0
# normal[190:212] <- 0
# normal[70:90] <- 0
# data <- cbind(data, normal)
# 
# sp500dataset <- data[ , c(-2,-3,-4,-5,-6)]
# names(sp500dataset)[c(1,2)] <- c("date", "sp500")
#_______________________________________________________________________________
# fit <- with(data, lm(sp500 ~ claims + corecpi + credit + disposincome + durable 
#                      + fedrate + housestarts + IPI + newhome + neworders + payroll
#                      + perincome + PMI + retail + unemployment + normal))
# tempdata <- data[-276 : -780 ,]
# tempdata <- tempdata[, -2:-7]
# tempdata <- left_join(tempdata, sp500[, -2:-6], by = "Date")
# tempdata <- tempdata[, c(1, 18, 2:17)]
# names(tempdata)[2] <- "sp500"
# names(usd)[1] <- "Date"
# tempdata <- left_join(tempdata, usd, by = "Date")
# names(tempdata)[19] <- "usd"
# tempdata <- tempdata[, c(1:17, 19, 18)]
# 
# cor(tempdata[3:19])
# 
# tempdata1 <- tempdata[, c(-4, -6, -7, -9, -16, -10, -14)]
# tempdata1 <- left_join(tempdata1, sporigin[, -2:-6], by = "Date")
# names(tempdata1)[13] <- "spclose"
# names(tempdata1)[2] <- "sppercent"
# tempdata1 <- tempdata1[, c(1, 13, 2:12)]
data <- data[, c(-2,-3,-4,-5,-6)]
names(data)[2] <- "spclose"
names(usd)[1] <- "Date"
data <- left_join(data, usd, by = "Date")
names(data)[19] <- "usd"
data <- data[, c(1:17, 19, 18)]


fit2 <- with(data, lm(spclose ~ claims + credit + fedrate + newhome + neworders + payroll + PMI + unemployment + usd))

par(mfrow = c(2, 4), oma = c(0,0,0,0))
with(data, plot(spclose ~ claims))
with(data = , plot(spclose ~ credit))
with(data = , plot(spclose ~ fedrate))
with(data, plot(spclose ~ newhome))
with(data, plot(spclose ~ neworders))
with(data, plot(spclose ~ payroll))
with(data, plot(spclose ~ PMI))
with(data, plot(spclose ~ unemployment))
with(data, plot(spclose ~ usd)) #added

par(mfrow = c(2, 4), oma = c(0,0,0,0))
with(data, plot(sppercent ~ claims))
with(data, plot(sppercent ~ credit))
with(data, plot(sppercent ~ fedrate))
with(data, plot(sppercent ~ newhome))
with(data, plot(sppercent ~ neworders))
with(data, plot(sppercent ~ payroll))
with(data, plot(sppercent ~ PMI))
with(data, plot(sppercent ~ unemployment))
with(data, plot(sppercent ~ usd)) #added

plot(fitted.values(fit2), rstandard(fit2))
qqnorm(rstandard(fit2))
abline(0,1)
plot(fitted.values(fit1), rstandard(fit1))
qqnorm(rstandard(fit1))
abline(0,1)

#################################################################################

library(MASS)
boxcox(fit2)
rstandard(fit2)
rstudent(fit2)
library(MPV)
PRESS(fit2)


plot(fitted.values(fit2), rstudent(fit2), xlab = "fitted values", ylab = "Studentized residuals", main = "Residual plot")
fit3 <- lm(I(log(spclose))~claims + credit + fedrate + newhome + neworders + payroll + PMI + unemployment +usd, data = data)
summary(fit3)

plot(fitted.values(fit3), rstandard(fit3))
qqnorm(rstandard(fit3))
abline(0,1)




X <- with(data, as.matrix(cbind(claims, credit, fedrate, newhome, neworders, payroll, PMI, unemployment, usd)))
H <- X%*%solve(t(X)%*%X)%*%t(X)
H.order <- cbind(c(1:275), diag(H))
order(H.order)
2*(8+1)/275  #cut off of hii: 0.06545455
d <- cooks.distance(fit2)
d[order(d, decreasing = TRUE)]
pf(0.5,9,266) #cut off of Cook's D: 0.1260753
dfbetas(fit2)
dffits(fit2)
2/sqrt(275) #cut off of DFBETAS: 0.1206045

#largest influence points
reg.rstand <- abs(rstandard(fit2))
influence.pts <- order(reg.rstand)

#cook's distance
d <- cooks.distance(fit2)
cook.d <- order(d)

#dfbetas
dfbeta.fit2 <- dfbetas(fit2)
#cutoff = 2/sqrt(n) = 2/sqrt(275) = 0.1206045

#dffits
dffits.fit2 <- dffits(fit2)
#cutoff = 2sqrt(p/n) = 2sqrt(10/275) = .381385

#covariance ratio
covratio.fit2 <- covratio(fit2)
#cutoff = 1 +- 3p/n = 1 +- 3(10)/275 = (0.8909091, 1.109091)