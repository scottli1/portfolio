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
# data <- data[, c(-2,-3,-4,-5,-6)]
# names(data)[2] <- "spclose"
# names(usd)[1] <- "Date"
# data <- left_join(data, usd, by = "Date")
# names(data)[19] <- "usd"
# data <- data[, c(1:17, 19, 18)]
# data <- data[-276 : -780 ,]


# Computing variance inflation factors
# First, standardize matrix with predictors in its columns:
#cor(data[3:18])
#Xvif <- solve(t(as.matrix(data[3:18])) %*% as.matrix(data[3:18]))
Xvif <- as.matrix(data[3:18])

ssi <- function(x){
    return(sqrt(sum((x - mean(x))^2)))
}
allssi <- apply(Xvif, 2, ssi)
W <- scale(Xvif, scale = allssi) # subtracts column means and divides by the square root of the sums of squares
t(W)%*%W
cor(Xvif)
VIF <- solve(cor(Xvif))
# Compare to variance inflation factors via R^2
fitall <- with(data, lm(spclose ~ claims + corecpi + credit + disposincome + durable + fedrate + housestarts + newhome + IPI + newhome + neworders + payroll + perincome + PMI + retail + unemployment + usd))
r2 <- summary(fitall)$r.squared
1/(1-r2)

fit2 <- with(data, lm(spclose ~ claims + credit + fedrate + newhome + neworders + payroll + PMI + unemployment + usd))
X2vif <- as.matrix(data[3:18][c(-4, -6, -7, -9, -16, -10, -14)])
VIF2 <- solve(cor(X2vif))

par(mfrow = c(2, 4), oma = c(0,0,0,0))
with(data, plot(spclose ~ claims))
with(data, plot(spclose ~ corecpi))
with(data, plot(spclose ~ credit))
with(data, plot(spclose ~ disposincome))
with(data, plot(spclose ~ durable))
with(data, plot(spclose ~ fedrate))
with(data, plot(spclose ~ housestarts))
with(data, plot(spclose ~ IPI))
with(data, plot(spclose ~ newhome))
with(data, plot(spclose ~ neworders))
with(data, plot(spclose ~ payroll))
with(data, plot(spclose ~ perincome))
with(data, plot(spclose ~ PMI))
with(data, plot(spclose ~ retail))
with(data, plot(spclose ~ unemployment))
with(data, plot(spclose ~ usd)) #added

with(data, plot(log(spclose) ~ claims))
with(data, plot(log(spclose) ~ corecpi))
with(data, plot(log(spclose) ~ credit))
with(data, plot(log(spclose) ~ disposincome))
with(data, plot(log(spclose) ~ durable))
with(data, plot(log(spclose) ~ fedrate))
with(data, plot(log(spclose) ~ housestarts))
with(data, plot(log(spclose) ~ IPI))
with(data, plot(log(spclose) ~ newhome))
with(data, plot(log(spclose) ~ neworders))
with(data, plot(log(spclose) ~ payroll))
with(data, plot(log(spclose) ~ perincome))
with(data, plot(log(spclose) ~ PMI))
with(data, plot(log(spclose) ~ retail))
with(data, plot(log(spclose) ~ unemployment))
with(data, plot(log(spclose) ~ usd)) #added

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


X <- with(data, as.matrix(cbind(neworders, PMI, unemployment, usd, retail)))
H <- X%*%solve(t(X)%*%X)%*%t(X)
H.order <- cbind(c(1:275), diag(H))
sort(H.order)
order(H.order)
2*(8+1)/275  #cut off of hii: 0.06545455
H.order[,1][H.order[,2] > 2*(5+1)/275]
data[,1][c(117, 203, 204, 205, 206)]

fit2 <- lm(spclose ~ neworders + PMI + usd
          + retail + unemployment, data = data)

d <- cooks.distance(fit2)
d.order <- cbind(c(1:275),d)
d.order[,1][d.order[,2] > pf(0.5,9,266)]
pf(0.5,9,266) #cut off of Cook's D: 0.1260753


#Judging from Cook's distance, only observation 270 seems to be influential 

dfbetas(fit2)
dffits(fit2)
2/sqrt(275) #cut off of DFBETAS: 0.1206045
dfbeta.fit2.order[,1][dfbeta.fit2.order[,c(2:9)]]   #???
which(dfbetas(fit2)[,2]> 2/sqrt(275))
data[,1][c(83,  95,  98, 103, 128, 129, 131, 132)]
which(dfbetas(fit2)[,3]> 2/sqrt(275))
data[,1][c(201, 202, 205, 271, 273, 274)] 
which(dfbetas(fit2)[,4]> 2/sqrt(275))
data[,1][103] 
which(dfbetas(fit2)[,5]> 2/sqrt(275))
data[,1][c(67, 236, 239, 261, 267, 268, 269, 270, 271, 273, 274, 275)]
which(dfbetas(fit2)[,6]> 2/sqrt(275))

#largest influence points
reg.rstand <- abs(rstandard(fit2))
influence.pts <- order(reg.rstand)

#cook's distance
d <- cooks.distance(fit2)
cook.d <- order(d)
which(d > 1)
#dfbetas
dfbeta.fit2 <- dfbetas(fit2)
#cutoff = 2/sqrt(n) = 2/sqrt(275) = 0.1206045

#dffits
dffits.fit2 <- dffits(fit2)
#cutoff = 2sqrt(p/n) = 2sqrt(10/275) = .381385

#covariance ratio
covratio.fit2 <- covratio(fit2)

#cutoff = 1 +- 3p/n = 1 +- 3(10)/275 = (0.8909091, 1.109091)