

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

pdf("SurveyByExperience.pdf",width=11,height=8.5)
par(mfrow=c(1,2),mar=c(6.5,5,2,2))
plot(d$Estimate~d$Familiarity,col="lightgrey",cex.axis=0.75,ylab="Estimate of proportion woodiness",xlab="",axes=FALSE)
axis(side=2,cex=1.5,cex.lab=2)
text(seq(1,4,by=1),y=-5, labels = levels(d$Familiarity),srt = -45, xpd = TRUE,adj = c(0, NA),cex=0.9)
abline(h=47.5,col="red",lty=1)
polygon(x=c(0,0,5,5),y=c(46,49,49,46),col=rgb(red=1,green=0,blue=0,alpha=0.2),border=NA)


plot(d$Estimate~d$Training,col="lightgrey",cex.axis=0.75,ylab="",xlab="",axes=FALSE)
axis(side=2)
par(srt=2)
#axis(side=1,at=1:5,labels=levels(d$Training))
xl<-c("Postgrad","Part Postgrad","Undergrad","Part Undergrad","None")
text(1:5,y=-3, labels = xl, srt = 0, ,srt = -45, xpd = TRUE,adj = c(0, NA),cex=0.9)
abline(h=47.5,col="red",lty=1)
polygon(x=c(0,0,6,6),y=c(46,49,49,46),col=rgb(red=1,green=0,blue=0,alpha=0.2),border=NA)

dev.off()



unique(d$In.what.country.did.you.receive.your.biology.botany.training.)