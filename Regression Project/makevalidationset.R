# create validation val frame

names(sp500)[1] <- "DATE"
names(PMI)[1] <- "DATE"
val <- left_join(sp500[c(1,5)], claims, by = "DATE")
val <- left_join(val, corecpi, by = "DATE")
val <- left_join(val, credit, by = "DATE")
val <- left_join(val, disposincome, by = "DATE")
val <- left_join(val, durable, by = "DATE")
val <- left_join(val, fedrate, by = "DATE")
val <- left_join(val, housestarts, by = "DATE")
val <- left_join(val, IPI, by = "DATE")
val <- left_join(val, newhome, by = "DATE")
val <- left_join(val, neworders, by = "DATE")
val <- left_join(val, payroll, by = "DATE")
val <- left_join(val, perincome, by = "DATE")
val <- left_join(val, PMI, by = "DATE")
val <- left_join(val, retail, by = "DATE")
val <- left_join(val, unemployment, by = "DATE")
val <- left_join(val, usd, by = "DATE")
names(val) <- c("Date", "spclose",
                 "claims", "corecpi", "credit", "disposincome", "durable",
                 "fedrate", "housestarts", "IPI", "newhome", "neworders", "payroll",
                 "perincome", "PMI", "retail", "unemployment", "usd")

val <- val[-24,]
val <- val[c(-1,-2),]

val <- val[order(val[[1]]), ]
rownames(val) <- 1:length(val[,1])

with(val, plot(spclose))
with(val, plot(log(spclose)))


data[,20] <- rep(0, length.out = length(data[,19]))
data[,20] <- 1*(data[,19] == 0)
data[,19] <- data[,20]
data <- data[, -20]

val[, 19] <- rep(0, length.out = length(val[,1]))
names(val)[19] <- "normal"

data[,19] <- as.factor(data[,19])
val[, 19] <- as.factor(val[, 19])


