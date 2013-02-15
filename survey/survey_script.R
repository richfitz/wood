

if (Sys.getenv("USER")=='wcornwell') setwd("/Users/wcornwell/Documents/data/howMuchWoodiness/wood/survey/")
d <- read.csv(file="Plant_survey_final.csv", as.is=TRUE)


# remove timestamp column
d <- d[,-c(1,5)]

# change the colnames
colnames(d) <- c("Estimate", "Familiarity", "Training", "Country")

levels(factor(d$Estimate))
as.numeric(d$Estimate)

d$Estimate
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

#d$Estimate[d$Estimate < 1] <- d$Estimate[d$Estimate < 1] * 100

d$Estimate

# get mean
m <- mean(d$Estimate)

# get median
med <- median(d$Estimate)

sd <- sd(d$Estimate)

## Convert estimates to normal using logit transformation
d$Estimate.logit <- boot::logit(d$Estimate / 100)
res <- lm(Estimate.logit ~ Training + Familiarity, data=d)


summary(res)
anova(res)
plot(d$Estimate~d$Familiarity)


