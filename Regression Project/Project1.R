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

fit1 <- with(tempdata1, lm(sppercent ~ claims + credit + fedrate + newhome + neworders + payroll + PMI + unemployment + usd))
fit2 <- with(tempdata1, lm(spclose ~ claims + credit + fedrate + newhome + neworders + payroll + PMI + unemployment + usd))

par(mfrow = c(2, 4), oma = c(0,0,0,0))
with(tempdata1, plot(spclose ~ claims))
with(tempdata1, plot(spclose ~ credit))
with(tempdata1, plot(spclose ~ fedrate))
with(tempdata1, plot(spclose ~ newhome))
with(tempdata1, plot(spclose ~ neworders))
with(tempdata1, plot(spclose ~ payroll))
with(tempdata1, plot(spclose ~ PMI))
with(tempdata1, plot(spclose ~ unemployment))
with(tempdata1, plot(spclose ~ usd))


par(mfrow = c(2, 4), oma = c(0,0,0,0))
with(tempdata1, plot(sppercent ~ claims))
with(tempdata1, plot(sppercent ~ credit))
with(tempdata1, plot(sppercent ~ fedrate))
with(tempdata1, plot(sppercent ~ newhome))
with(tempdata1, plot(sppercent ~ neworders))
with(tempdata1, plot(sppercent ~ payroll))
with(tempdata1, plot(sppercent ~ PMI))
with(tempdata1, plot(sppercent ~ unemployment))

plot(fitted.values(fit2), rstandard(fit2))
qqnorm(rstandard(fit2))
abline(0,1)
plot(fitted.values(fit1), rstandard(fit1))
qqnorm(rstandard(fit1))
abline(0,1)

fit3 <- with(tempdata1, lm(spclose ~ claims + fedrate + newhome + payroll + PMI + unemployment))
plot(fitted.values(fit3), rstandard(fit3))
qqnorm(rstandard(fit3))
abline(0,1)


