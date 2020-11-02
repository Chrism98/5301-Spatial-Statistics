# Loading libraries
library(spdep)
library(rgdal)
library(RColorBrewer)
library(car)
library(MASS)
library(robustbase)
library(olsrr)
library(classInt)

# Getting election data from files
elec<-readOGR(dsn="C:/Users/Chris/Documents/OSU/GEOG/5103/elec2016")

# Creating subsets for roman catholics and evangelic protestants and ignoring all potential zeros in data set
elecRC<-subset(elec,PctRC>0)
summary(elecRC)
elecEV<-subset(elec,PctEvang>0)
summary(elecEV)
elec.df<-data.frame(elec)

# setting each data set as a variable to simplify and graphing as Q-Q plots
RC<-elecRC$PctRC
summary(RC)
qqnorm(RC)
qqline(RC)
EV<-elecEV$PctEvang
summary(EV)
qqnorm(EV)
qqline(EV, col="red")

# Plots and histograms to visualize data, transforming data to normalize
plot(elec)
plot(elecRC)
plot(elecEV)
hist(RC)

hist(elecEV$PctEvang)

# Plotting relationships between confounding variables against religion (Roman-Catholics)
hist(elecRC$PctWht)
cor(RC,elecRC$PctWht) #r=0.1069
lm(elecRC$PctWht~RC)
plot(RC,elecRC$PctWht,
     main="Percent Roman Catholic vs. Percent White", xlab="% Roman Catholic", 
     ylab="% White", pch=20)
abline(82.48829,0.09333, col="red") #all abline values below taken from the linear model ran per each unique variable (example line 42)

hist(elecRC$PctHisp)
cor(RC,elecRC$PctHisp) #r=0.2690
lm(elecRC$PctHisp~RC)
     main="Percent Roman Catholic vs. Percent Hispanic", xlab="% Roman Catholic", 
     ylab="% Hispanic", pch=20)
abline(5.0198,0.2152, col="red")

hist(elecRC$EduBach)
LRCEB<-log(elecRC$EduBach) #removes tail from curve
cor(RC,LRCEB) #r=0.3046
lm(LRCEB~RC)
plot(RC,LRCEB,
     main="Percent Roman Catholic vs. log(Percent with at least a Bachelor's Degree)", xlab="% Roman Catholic", 
     ylab="log(% with at least a B.A. or B.S.)", pch=20)
abline(2.850074,0.006726, col="red") 


hist(elecRC$MedIncome)
LRCMI<-log(elecRC$MedIncome) #removes tail from curve
cor(RC,LRCMI) #r=0.3632
lm(LRCMI~RC)
plot(RC,LRCMI,
     main="Percent Roman Catholic vs. log(Median Income)", xlab="% Roman Catholic", 
     ylab="log(Median Income)", pch=20)
abline(10.69737,0.004897,col="red")

# Plotting a map of Roman Catholics per county
rc.jenks<-classIntervals(elec$PctRC, n=5, style="jenks")
rc.jenks
paletteRC <- brewer.pal(5, "Reds")
col.rc <- findColours(rc.jenks, paletteRC)
windows()
plot(elec,
     col=col.rc,
     pch=20)

legend("topleft", fill = attr(col.rc, "palette"),
       legend = names(attr(col.rc, "table")), bty = "n")
mtext("Percent Roman Catholic", cex=1.2, side=3, line=1)

cor(elecRC$PctHisp,elecRC$PctGOP) #checking for correlation between hispanic voters and GOP voters (r=-0.2071)

# Plotting relationships between confounding variables against religion (Evangelical Protestants)
hist(elecEV$PctPov) 
LEVPP<-log(elecEV$PctPov) #removes tail from curve
cor(EV,LEVPP) #r=0.4823
lm(LEVPP~EV)
plot(EV,LEVPP,
     main="Percent Evangelical Protestants vs. log(Percent in Poverty)", xlab="% Evangelical Protestant", 
     ylab="log(% in Poverty)", pch=20)
abline(2.322123, 0.008427, col="blue")

hist(elecEV$Unemp)
LEVU<-log(elecEV$Unemp) #removes tail from curve
cor(EV,LEVU) #r=0.3133
lm(LEVU~EV)
plot(EV,LEVU,
     main="Percent Evangelical Protestants vs. log(Percent Unemployed)", xlab="% Evangelical Protestant", 
     ylab="log(% Unemployed)", pch=20)
abline(1.394563,0.004647, col="blue")

hist(elecEV$EduLTHS)
LEVDL<-log(elecEV$EduLTHS) #removes tail from curve
cor(EV,LEVDL) #r=0.5434
lm(LEVDL~EV)
plot(EV,LEVDL,
     main="Percent Evangelical Protestants vs. log(Percent with Less than a High School Education)", xlab="% Evangelical Protestant", 
     ylab="log(% < High School Education)", pch=20)
abline(2.0447,0.01143, col="blue")

hist(elecEV$MedIncome)
LEVMI<-log(elecEV$MedIncome) #removes tail from curve
cor(EV,LEVMI) #r=-0.4910
lm(LEVMI~EV)
plot(EV,LEVMI,
     main="Percent Evangelical Protestants vs. log(Median Income)", xlab="% Evangelical Protestant", 
     ylab="log(Median Income)", pch=20)
abline(11.015083,-0.005354, col="blue")

hist(elecEV$PctWht)
cor(EV,elecEV$PctWht) #r=-0.2400
lm(elecEV$PctWht~EV)
plot(EV,elecEV$PctWht,
     main="Percent Evangelical Protestants vs. Percent White", xlab="% Evangelical Protestant", 
     ylab="% White", pch=20)
abline(91.5392,-0.1751, col="blue")

# Plotting a map of Evangelical Protestants per county
ev.jenks<-classIntervals(elec$PctEvang, n=5, style="jenks")
ev.jenks
paletteEV <- brewer.pal(5, "Blues")
col.ev <- findColours(ev.jenks, paletteEV)
windows()
plot(elec,
     col=col.ev,
     pch=20)

legend("topleft", fill = attr(col.ev, "palette"),
       legend = names(attr(col.ev, "table")), bty = "n")
mtext("Percent Evangelical Protestant", cex=1.2, side=3, line=1)

cor(elecEV$PctBlk,elecEV$PctGOP) #checking for correlation between black voters and GOP voters (r=-0.4271)

# Building a model to relate the factors that exaplin voting habits of Roman Catholics
rc1<-lm(PctRC~PctWht+PctHisp+LRCEB+LRCMI, data=elecRC)
summary(rc1) # all values significant!
anova(rc1) # all values significant!
plot(rc1)

windows()
# Check for heteroskedasticity
plot(rc1$fitted.values, rc1$residuals,
     main="Roman Catholic Residuals vs. Fitted", xlab="Fitted values", 
     ylab="Residuals", pch=20); abline(h=0, col="red") # +heteroskedasticity

# Checking for high VIF values
rc.mult<-lm(PctRC~PctWht+PctHisp+LRCEB+LRCMI, data=elecRC)
summary(rc.mult)
vif(rc.mult) # all values at or below 2!

# Building a model to relate the factors that exaplin voting habits of Evangelical Protestants
ev1<-lm(PctEvang~LEVPP+LEVU+LEVDL+LEVMI+PctWht, data=elecEV)
summary(ev1) # LEVU not significant
anova(ev1)
# Checking for heteroskedasticity
plot(ev1$fitted.values, ev1$residuals,
     main="Evangelical Protestant Residuals vs. Fitted", xlab="Fitted values", 
     ylab="Residuals", pch=20); abline(h=0, col="blue") 

ev2<-lm(PctEvang~LEVDL+LEVMI+LEVPP+PctWht, data=elecEV)
summary(ev2) #all values significant, better model!
anova(ev2)
anova(ev1,ev2)
# Checking for heteroskedasticity
plot(ev2$fitted.values, ev2$residuals,
     main="Evangelical Protestant Residuals vs. Fitted", xlab="Fitted values", 
     ylab="Residuals", pch=20); abline(h=0, col="blue") # +heteroskedasticity
# Checking for high VIF values
ev.mult<-lm(PctEvang~LEVDL+LEVMI+LEVPP+PctWht, data=elecEV)
summary(ev.mult)
vif(ev.mult) #LEVPP value is concerning, but workable. LEVMI value a little less concerning but still some concern.

# Performing robust regression on RC data
rc.rob<-lmrob(PctRC~PctWht+PctHisp+LRCEB+LRCMI, data=elecRC)
confint(rc.rob)
confint(rc1)
plot(rc.rob$fitted.values, rc.rob$residuals, main="Roman Catholic Residuals vs. Predictor", 
     xlab="Predictor", ylab="Residuals", pch=20); abline(h=0, col="red") # similar results to previous regression

# Performing robust regression on EV data
ev.rob<-lmrob(PctEvang~LEVDL+LEVMI+LEVPP+PctWht, data=elecEV)
confint(ev.rob)
confint(ev2)
plot(ev.rob$fitted.values, ev.rob$residuals, main="Evangelic Protestant Residuals vs. Predictor",
     xlab="Predictor", ylab="Residuals", pch=20); abline(h=0, col="blue") #similar results to previous regression

# Performing OSL regression for outliers (Roman Catholics)
lmRC<-lm(PctRC~PctWht+LRCEB+LRCMI, data=elecRC)
windows()
ols_plot_resid_stand(lmRC)
windows()
ols_plot_resid_stud(lmRC)
windows()
ols_plot_dfbetas(lmRC)
windows()
ols_plot_cooksd_chart(lmRC)

# Performing OSL regression for outliers (Evangelical Protestants)
lmEV<-lm(PctEvang~LEVDL+LEVMI+LEVPP+PctWht, data=elecEV)
windows()
ols_plot_resid_stand(lmEV)
windows()
ols_plot_resid_stud(lmEV)
windows()
ols_plot_dfbetas(lmEV)
