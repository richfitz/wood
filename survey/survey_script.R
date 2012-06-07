d <- read.csv(file="Woodyness_survey.csv", as.is=TRUE)

# remove timestamp column
d <- d[,-c(1,5)]

# change the colnames
colnames(d) <- c("Estimate", "Familiarity", "Training", "Country")

# change descriptions to numerical

changeFamiliarity <- function(x){
	
	if (x == "Very Familiar"){
		x <- 1
	} 
	if (x == "Familiar"){
		x <- 2
	} 
	if (x == "Somewhat Familiar"){
		x <- 3
	} 
	if (x == "What's a Plant?"){
		x <- 4
	} 
	return(x)
}

changeTraining <- function(x){
	
	if (x == "Postgraduate degree in botany or a related field"){
		x <- 1
	}
	if (x == "Partially complete postgraduate degree in botany or a related field"){
		x <- 2
	}
	if (x == "Undergraduate degree in botany or a related field"){
		x <- 3
	}
	if (x == "Some botany courses at either an undergraduate or postgraduate level"){
		x <- 4
	}
	if (x == "No formal training in botany"){
		x <- 5
	}
	return(x)
}

d[,"Familiarity"] <- sapply(d[,"Familiarity"], function(x) changeFamiliarity(x))

d[,"Training"] <- sapply(d[,"Training"], function(x) changeTraining(x))

# get mean
m <- mean(d$Estimate)

# get median
med <- median(d$Estimate)

sd <- sd(d$Estimate)


# Convert estimates to normal using logit transformation

d[,"Estimate"] <- sapply(d[,"Estimate"], function(x) return(x/100))
d[,"Estimate"] <- logit(d[,"Estimate"])

# fit glm
formula <- Estimate ~ Training + Familiarity

res <- glm(formula, data=d)

summary(res)


