require(ggplot2)
require(vcd)

# read data as csv file
mydata <- read.csv("HW1.csv", header = TRUE)

# convert question 13 response to factor
mydata$Question.13 <- as.factor(mydata$Question.13)

# view summary of data, there are no missing values
summary(mydata)
ftable(with(mydata, table(College, Question.13)))
ftable(with(mydata, table(Student.Level, Question.13)))

ftable(round(prop.table(with(mydata, table(College, Question.13)),2), digits = 3)) # column frequency
ftable(round(prop.table(with(mydata, table(Student.Level, Question.13)),2), digits = 3)) # column frequency

ftable(with(mydata, table(College, Student.Level, Question.13))) # frequency table
ftable(xtabs(~College + Student.Level + Question.13, data = mydata)) # same frequency table

ftable(round(prop.table(xtabs(~College + Student.Level + Question.13, data = mydata)), digits = 3))

# some barplots of the variables
ggplot(mydata, aes(Question.13)) + geom_bar()
ggplot(mydata, aes(College)) + geom_bar()
ggplot(mydata, aes(Student.Level)) + geom_bar()

ggplot(mydata, aes(x = Question.13, fill = College)) + geom_bar()
ggplot(mydata, aes(x = Question.13, fill = Student.Level)) + geom_bar()

ggplot(mydata, aes(x = Question.13, fill = College)) + geom_bar(position = "dodge")
ggplot(mydata, aes(x = Question.13, fill = Student.Level)) + geom_bar(position = "dodge")

ggplot(mydata, aes(x = Student.Level , fill = College)) + geom_bar(position = "dodge")
ggplot(mydata, aes(x = College, fill = Student.Level)) + geom_bar(position = "dodge")


# create frequency table as data frame
fdata <- as.data.frame(xtabs(~College + Student.Level + Question.13, data = mydata)) 

# some random plots to try
assocplot(with(mydata, table(Student.Level, Question.13)))
mosaicplot(with(mydata, table(Student.Level, Question.13)))
spineplot(with(mydata, table(Student.Level, Question.13)))
mosaic(Freq ~ Student.Level + Question.13, data = fdata, shade = TRUE, legend = FALSE)
mosaic(Freq ~ Question.13 + College, data = fdata, shade = TRUE, legend = FALSE)
