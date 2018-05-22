set.seed(50)

set1 <- sample(1:275, size = 50)
set2 <- sample(1:275, size = 50)
set3 <- sample(1:275, size = 50)

fit1.0 <- lm(log(spclose) ~ claims + fedrate + newhome + neworders + PMI + usd + normal +
               claims*normal + newhome*normal + neworders*normal +
               PMI*normal + usd*normal, data = data)

fit2.0 <- lm(I(log(spclose)) ~ claims + fedrate + newhome + neworders + PMI + usd, data = data)

fit3.0 <- lm(spclose ~ newhome + neworders + PMI + usd
                  + retail + unemployment + normal + normal*newhome, data = data)

fit4.0 <- lm(I(sqrt(spclose)) ~ newhome + neworders + PMI + usd + retail + unemployment
                        + normal + normal*newhome, data = data)

fit5.0 <- lm(spclose ~ neworders + PMI + usd
             + retail + unemployment, data = data)



fit1.1 <- lm(log(spclose) ~ claims + fedrate + newhome + neworders + PMI + usd + normal +
                 claims*normal + newhome*normal + neworders*normal +
                 PMI*normal + usd*normal, data = data[-set1,])

fit2.1 <- lm(I(log(spclose)) ~ claims + fedrate + newhome + neworders + PMI + usd, data = data[-set1,])

fit3.1 <- lm(spclose ~ newhome + neworders + PMI + usd
             + retail + unemployment + normal + normal*newhome, data = data[-set1,])

fit4.1 <- lm(I(sqrt(spclose)) ~ newhome + neworders + PMI + usd + retail + unemployment
             + normal + normal*newhome, data = data[-set1,])


fit5.1 <- lm(spclose ~ neworders + PMI + usd
             + retail + unemployment, data = data[-set1,])


p1 <- predict(fit1.1, data[set1, ]) - log(data[set1, ][, 2])
pp1 <- exp(predict(fit1.1, data[set1, ])) - (data[set1, ][, 2])
table(abs(p1) > 0.1068 * 2)
f1 <- predict(fit1.0, val) - log(val[, 2])
ff1 <- exp(predict(fit1.0, val)) - (val[, 2])

p2 <- predict(fit2.1, data[set1, ]) - log(data[set1, ][, 2])
pp2 <- exp(predict(fit2.1, data[set1, ])) - (data[set1, ][, 2])
table(abs(p2) > 0.1155 * 2)
f2 <- predict(fit2.0, val) - log(val[, 2])
ff2 <- exp(predict(fit2.0, val)) - (val[, 2])

p3 <- predict(fit3.1, data[set1, ]) - data[set1, ][, 2]
table(abs(p3) > 100.7 * 2)
f3 <- predict(fit3.0, val) - (val[, 2])
table(abs(f3) > 100.7 * 2)

p4 <- (predict(fit4.1, data[set1, ])) - sqrt(data[set1, ][, 2])
pp4 <- (predict(fit4.1, data[set1, ]))^2 - (data[set1, ][, 2])
table(abs(p4) > 1.44 * 2)
f4 <- predict(fit4.0, val) - sqrt(val[, 2])
ff4 <- (predict(fit4.0, val))^2 - (val[, 2])

p5 <- predict(fit5.1, data[set1, ]) - data[set1, ][, 2]
table(abs(p5) > 114.4 * 2)
f5 <- predict(fit5.0, val) - (val[, 2])
table(abs(f5) > 114.4 * 2)

cbind(sum(pp1^2)/50, sum(pp2^2)/50, sum(p3^2)/50, sum(pp4^2)/50, sum(p5^2)/50)
cbind(sum(ff1^2)/21, sum(ff2^2)/21, sum(f3^2)/21, sum(ff4^2)/21, sum(f5^2)/21)

cp <- 3671363/13062-275+2*(12-1)
press.value1 <- PRESS(fit5.0)
summary(fit5.0)
anova(fit5.0)
SSt1 <- sum(anova(fit5.0)$Sum)
rsquare.prediction1 <- 1- press.value1/SSt1

plot(fitted.values(fit5.0), rstandard(fit5.0), main = "fit5.0")
qqnorm(rstandard(fit5.0))
abline(0,1)

fit3.3 <- lm(spclose ~ newhome + neworders + PMI + usd
             + retail + unemployment, data = data)

predict(fit, val) - (val[, 2])
