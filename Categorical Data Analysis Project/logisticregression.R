library("DescTools") # statistical stuff

# read the sote data 
data <- fread('SOTE-data.csv')

# take variables of interest
dataset <- data[,c(2,5,6,7,8,9,10,12,26,27,29,30)]
dataset$`Grading Basis` <- as.factor(dataset$`Grading Basis`)
dataset$College <- as.factor(dataset$College)
dataset$`Student Level` <- as.factor(dataset$`Student Level`)
dataset$Component <- as.factor(dataset$Component)
dataset$`Intruction Mode` <- as.factor(dataset$`Intruction Mode`)
dataset$diff <- grade_data$difference

summary(dataset)

# remove no response and N/A responses from Question 13
dataset <- dataset[dataset$`Question 13` != 0 & dataset$`Question 13` != 6,]

# convert Question 13 into binary 
dataset$`Question 13` <- ifelse(dataset$`Question 13` < 5, 0 , 1)

plot(dataset$`Total Enrollment`, dataset$`Question 13`)

dataset$prop <- dataset$`Total Enrollment`/dataset$`Enrollment Cap`
plot(dataset$prop, dataset$`Question 13`)

plot(dataset$`Question 13`, dataset$`Student Level`)

# fit a glm model 
model <- glm(`Question 13` ~ `Student Level` + College 
             , data = dataset, family = binomial)

glm.diag.plots(model)

# standardized deviance residuals
sdres <- rstandard(model)

par(mar = c(4, 4, 1, 1))
plot(sdres, pch = 20, col = "palevioletred", cex = 0.5, ylim = c(-2, 2), axes = FALSE, ann = FALSE)
box(col = "gray60")
axis(1, col = "gray70")
axis(2, col = "gray70")
title(xlab = "Index", ylab = "Standardized deviance residual")
#text(650, sdres[650], "650", adj = 1.2, cex = 1)
identify(sdres)

outliers <- dataset[c(15691, 17556, 31805, 31806, 38494, 39545, 44486, 49691, 53224, 65477, 65540, 79885, 89443, 96639, 102468, 136726, 165116, 165161), c(2,3,9)]
sum(dataset$College == "Underg")

prop.table(table(dataset[dataset$College=="Underg",]$`Question 13`))

prop.table(table(dataset[dataset$College=="Engine",]$`Question 13`))

prop.table(table(dataset[dataset$College=="Engine",]$`Question 13`))

library(hnp)
model.hnp <- hnp(model, sim = 99, conf = 0.95, plot = FALSE)


par(mar = c(4, 4, 1, 1), mfrow = c(1,2))
plot(sdres, pch = 20, col = "palevioletred", cex = 0.5, ylim = c(-2, 2), axes = FALSE, ann = FALSE)
box(col = "gray60")
axis(1, col = "gray70")
axis(2, col = "gray70")
title(xlab = "Index", ylab = "Standardized deviance residual")

plot(model.hnp, lty = c(2, 1, 2), pch = 20, cex = 0.3, col=c(2,2,2,"gray50"),
     ylim = c(0.5, 2), axes = FALSE, ann = FALSE)
box(col = "gray60")
axis(1, col = "gray70")
axis(2, col = "gray70")
title(xlab = "Half-normal quantile", ylab = "Residual quantiles")

HosmerLemeshowTest(fitted(model), obs = dataset$`Question 13`)
HosmerLemeshowTest(fitted(model), obs = dataset$`Question 13`, ngr = 10)

logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

logit2prob(-0.333951 + -0.664595) 
logit2prob(-0.333951 + -0.301341) 
logit2prob(-0.333951 + 0.102748) 
logit2prob(-0.333951 + -0.351901) 

logit2prob(-0.333951 + -0.664595 + 0.265950) 
logit2prob(-0.333951 + -0.301341+ 0.265950) 
logit2prob(-0.333951 + 0.102748+ 0.265950) 
logit2prob(-0.333951 + -0.351901+ 0.265950) 

logit2prob(-0.333951 + -0.664595+ 0.352476) 
logit2prob(-0.333951 + -0.301341+ 0.352476) 
logit2prob(-0.333951 + 0.102748+ 0.352476) 
logit2prob(-0.333951 + -0.351901+ 0.352476) 

logit2prob(-0.333951 + -0.664595 + 0.494796) 
logit2prob(-0.333951 + -0.301341 + 0.494796) 
logit2prob(-0.333951 + 0.102748 + 0.494796) 
logit2prob(-0.333951 + -0.351901 + 0.494796) 

logit2prob(-0.333951 + -0.664595 + 0.481934) 
logit2prob(-0.333951 + -0.301341 + 0.481934) 
logit2prob(-0.333951 + 0.102748 + 0.481934) 
logit2prob(-0.333951 + -0.351901 + 0.481934) 

logit2prob(-0.333951) 
logit2prob(-0.333951+0.265950) 
logit2prob(-0.333951+0.352476) 
logit2prob(-0.333951+0.494796) 
logit2prob(-0.333951+0.481934) 

logit2prob(model$coefficients)

library(ryouready)
nom.uncertainty(table(data$`Student Level`, data$College)[c(-5,-6),-8])
