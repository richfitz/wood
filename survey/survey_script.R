d <- read.csv(file="Woodyness_survey.csv", as.is=TRUE)

# remove timestamp column
d <- d[,-c(1,5)]

# change the colnames
colnames(d) <- c("Estimate", "Familiarity", "Training", "Country")

# change descriptions to numerical

lvl.familiarity <- c("Very Familiar", "Familiar", "Somewhat Familiar",
                     "What's a Plant?")
lvl.training <-
  c("Postgraduate degree in botany or a related field",
    "Partially complete postgraduate degree in botany or a related field",
    "Undergraduate degree in botany or a related field",
    "Some botany courses at either an undergraduate or postgraduate level",
    "No formal training in botany")             
             
d$Familiarity <- factor(d$Familiarity, lvl.familiarity, ordered=TRUE)
d$Training <- factor(d$Training, lvl.training, ordered=TRUE)

d$Estimate[d$Estimate < 1] <- d$Estimate[d$Estimate < 1] * 100

# get mean
m <- mean(d$Estimate)

# get median
med <- median(d$Estimate)

sd <- sd(d$Estimate)

## Convert estimates to normal using logit transformation
d$Estimate.logit <- boot::logit(d$Estimate / 100)
res <- lm(Estimate.logit ~ Training + Familiarity, data=d)


summary(res)


