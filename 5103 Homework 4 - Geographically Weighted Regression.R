#loading libraries and setting working directory
library(spdep)
library(rgdal)
library(RColorBrewer)
library(car)
library(MASS)
library(robustbase)
library(olsrr)
library(pscl)
library(classInt)
library(graphicsQC)
library(lmtest)
library(sp)
library(spatstat)
library(spatialreg)
library(GWmodel)
library(ggplot2)

setwd("~/OSU/GEOG/5103/Homework 4 Files")
cov<-readOGR(dsn="C:/Users/Chris/Documents/OSU/GEOG/5103/COVID19")
cov<-cov[-c(1412),] #removing NA value
cov.df<-data.frame(cov)
summary(cov)

#creating a subset without counties with zero reported cases
cc<-subset(cov,confrmd>0)
cc.df<-data.frame(cc)
ccc<-cc$confrmd
windows()
hist(ccc, breaks="fd",col="Gray", prob=T, xlim = c(0, 400),
     main = "Distribution of Confirmed Cases",
     xlab = "Confirmed Cases", ylab = "Density")
lines(density(ccc, adjust=1.7), col="red", lwd=2, xlim=(c(1,400))) #poisson distribution; needs generalized linear model


#finding correlation coefficients for data sets
cor(cc.df[,c(10,4:6)],method="spearman")
cor(cc.df[,c(10,11:45)],method="spearman") 

#building a model based on correlations (poisson)
covid1<-glm(ccc~density+pctwht+pctblk+pctasin+pcttrns+migratn+RUCC,
            family="poisson", offset=log(pop_tot), data=cc)
summary(covid1)
vif(covid1)

#negative binomial model
covid1.nb<-glm.nb(covid1)
summary(covid1.nb)
sigma(covid1.nb)

#building a second model
covid2<-glm(ccc~ pct65+pctCOPD+pctchls+pcthypr+log(density+0.001)+
                 pcttrns+log(incm_cp+0.001)+pctsmok+log(pctvccn+0.001)+log(migratn+0.001)+RUCC,
                 family="poisson", offset=log(pop_tot), data=cc)
summary(covid2)
vif(covid2)

#negative binomial model
covid.nb<-glm.nb(covid2)
summary(covid.nb)
sigma(covid.nb)

covid.nb$deviance/covid.nb$df.residual
mean(residuals(covid.nb))
plot(covid.nb)

#regression diagnostics
me1 <- mean(residuals(covid.nb))
me1    #-0.355
sd1 <- sd(residuals(covid.nb))
sd1    #0.986
summary(residuals(covid.nb)) #not very symmetric

hist(residuals(covid.nb), col=8, probability=T,
     ylab='Density', main='Histogram of Residuals (covid.nb)',
     xlab='Residuals(covid.nb)')
box()
curve(dnorm(x, mean=me1, sd=sd1), from=-3, to=7, add=T,
      col='red', lwd=2)
leg.txt <- c("Residuals(covid.nb)", "Min. =-123984.4",
             "Max.=319607.9",
             "Mean =-0.364", "Median =-0.514", "Std. dev. =1.007")
#slight positive skew

#looking at qqplot for outliers
qqPlot(residuals(covid.nb), distribution="norm",
       xlab='', main='Quantile Comparison (Plot of covid.nb residuals)',
       envelope=.95, las=0, pch=NA, lwd=2, col="red",
       line="quartiles")
par(new=TRUE)
qqPlot(residuals(covid.nb), distribution="norm", envelope=FALSE,
       pch=1, cex=1, col="black") #definite postitive skew
par(new=FALSE)

#plotting to look for heteroskedasticity
plot(fitted(covid.nb), residuals(covid.nb), xlab="Fitted y", ylab= "Residuals",
     main="Plot of Residuals against Fitted y")
abline(h=0)

#looking for spatial autocorrelation
cov_NAD<-spTransform(cc, CRS("+init=epsg:3358"))
proj4string(cov_NAD)

##queen and rook
cov_nbq<-poly2nb(cov_NAD)
cov_nbr<-poly2nb(cov_NAD, queen=FALSE)

##nearest neighbors
coords<-coordinates(cov_NAD)
IDs<-row.names(as(cov_NAD, "data.frame"))
cov_kn1<-knn2nb(knearneigh(coords, k=1), row.names=IDs)
cov_kn2<-knn2nb(knearneigh(coords, k=2), row.names=IDs)
cov_kn4<-knn2nb(knearneigh(coords, k=4), row.names=IDs)

cov_kn1_w<- nb2listw(cov_kn1, zero.policy = TRUE)
cov_kn2_w<- nb2listw(cov_kn2, zero.policy = TRUE)
cov_kn4_w<- nb2listw(cov_kn4, zero.policy = TRUE)

##creating weights matricies
cov_nbq_w<- nb2listw(cov_nbq, zero.policy = TRUE)
cov_nbq_wb<-nb2listw(cov_nbq, style="B", zero.policy = TRUE)

cov_nbr_w<- nb2listw(cov_nbr, zero.policy = TRUE)
cov_kn1_w<- nb2listw(cov_kn1, zero.policy = TRUE)
cov_kn2_w<- nb2listw(cov_kn2, zero.policy = TRUE)
cov_kn4_w<- nb2listw(cov_kn4, zero.policy = TRUE)

moran.test(cov_NAD$confrmd, listw=cov_nbq_w, zero.policy = TRUE) #0.765
moran.test(cov_NAD$confrmd, listw=cov_nbr_w, zero.policy = TRUE) #0.766
moran.test(cov_NAD$confrmd, listw=cov_kn1_w, zero.policy = TRUE) #0.781
moran.test(cov_NAD$confrmd, listw=cov_kn2_w, zero.policy = TRUE) #0.759
moran.test(cov_NAD$confrmd, listw=cov_kn4_w, zero.policy = TRUE) #0.688

moran.test(covid.nb$resid, cov_kn1_w, zero.policy=TRUE) #0.239; certainly better but not great
#moran scatterplot and LISA map generated in GeoDa

#searching for average number of neighbors
covtest<- poly2nb(cc)
lapply(covtest, function(x) attr(covtest, "region.id")[x])

covspat1.gal<-nb2listw(poly2nb(cc),zero.policy=TRUE)
summary(covspat1.gal,zero.policy=TRUE)

#attaching data set
attach(cc.df)
class(cc.df)
names(cc.df)
row.names(cc.df)<-cc.df$COMM

#running moran test on residuals from above regression
moran.test(covid.nb$residuals,covspat1.gal, alternative="two.sided",zero.policy=TRUE) #I=0.189

#looking for potential spatial model
lm.LMtests(covid.nb, cov_kn2_w, test=c("LMerr", "LMlag", "RLMerr",
                             "RLMlag", "SARMA"))
#LMerr selected

##spatial regression
#spatial error
cov.err.eig<-errorsarlm(ccc~ pct65+pctCOPD+pctchls+pcthypr+log(density+0.001)+
                    pcttrns+log(incm_cp+0.001)+pctsmok+log(pctvccn+0.001)+log(migratn+0.001)+RUCC,
                    data=cc,cov_kn2_w, method="eigen", quiet=FALSE)
summary(cov.err.eig) #some spatial error

#looking at residuals
me2<-mean(residuals(cov.err.eig))
me2 #effectlively zero
sd2<-sd(residuals(cov.err.eig))
sd2 #1090.784
summary(residuals(cov.err.eig))
plot(residuals(cov.err.eig));abline(h=0, col="red") #reasonably symmetric

#building a histogram to visualize residual data
hist(residuals(cov.err.eig), col=8, probability = T,xlim = c(-5e+03, 5e+03),
     ylab="Density", main="Histogram of Residuals (Spatial Error Model)", xlab="Residuals (Spatial Error Model")
box()
curve(dnorm(x, mean=me2, sd=sd2), from=-5e+03, to=5e+03, add=T,
      col="red", lwd=2) #fairly normal

#looking at outliers
(sd(residuals(cov.err.eig)))*3 #3272.352
which(residuals(cov.err.eig)<=3272.352) #more outliers on the low side, to be expected
which(residuals(cov.err.eig)>=3272.352)

#checking residuals for spatial autocorrelation using Q1 matrix
names(cov.err.eig)
moran.test(cov.err.eig$residuals, covspat1.gal, alternative="two.sided",zero.policy=TRUE) #I=0.289

#spatial lag #####not necessary, just curious!
cov.lag.eig<-lagsarlm(ccc~ pct65+pctCOPD+pctchls+pcthypr+log(density+0.001)+
                         pcttrns+log(incm_cp+0.001)+pctsmok+log(pctvccn+0.001)+log(migratn+0.001)+RUCC,
                        data=cc,cov_kn1_w, method="eigen", quiet=FALSE)
summary(cov.lag.eig) #some issues with lag

#looking at residuals
me3<-mean(residuals(cov.lag.eig))
me3 #effectlively zero
sd3<-sd(residuals(cov.lag.eig))
sd3 #1053.28
summary(residuals(cov.lag.eig))
plot(residuals(cov.lag.eig));abline(h=0, col="red") #reasonably symmetric

#building a histogram to visualize residual data
hist(residuals(cov.lag.eig), col=8, probability = T,xlim = c(-5e+03, 5e+03),
     ylab="Density", main="Histogram of Residuals (Spatial Lag Model)", xlab="Residuals (Spatial Lag Model")
box()
curve(dnorm(x, mean=me3, sd=sd3), from=-5e+03, to=5e+03, add=T,
      col="red", lwd=2) #fairly normal
#almost identical results to the spatial error model

#looking at outliers
(sd(residuals(cov.lag.eig)))*3 #3159.84
which(residuals(cov.lag.eig)<=3159.84) #more outliers on the low side, to be expected
which(residuals(cov.lag.eig)>=3159.84)

#checking residuals for spatial autocorrelation using Q1 matrix
names(cov.lag.eig)
moran.test(cov.lag.eig$residuals, covspat1.gal, alternative="two.sided",zero.policy=TRUE) #I=0.210
#lag residuals less spatially autocorrelated than error residuals

#comparing autocorrelation in residuals
moran.test(covid.nb$residuals, covspat1.gal, alternative="two.sided",zero.policy = TRUE) #I=0.189
moran.test(cov.err.eig$residuals, covspat1.gal, alternative="two.sided",zero.policy = TRUE) #I=0.289
moran.test(cov.lag.eig$residuals, covspat1.gal, alternative="two.sided",zero.policy = TRUE) #I=0.210

#comparing AIC values
AIC(covid.nb) #25318.58
AIC(covid1.nb)#25237.38
AIC(cov.err.eig) #46284.01
AIC(cov.lag.eig) #46076.73

#plotting residuals across the different models
windows()
par(mfrow=c(1,2))
par(mar=c(1.2,1.2,1.2,1.2))
colors<-brewer.pal(5, "Blues")
color.cat.err<-classIntervals(cov.err.eig$residuals, n=5, style="quantile", dataPrecision=2)
colerr <- findColours(color.cat.err, colors)
plot(cov, col=colerr)
title('Map of Error Model Residuals')
legend('topleft', legend=c(names(attr(colerr, 'table'))), fill=c(attr(colerr, 'palette')), 
       title='Regression Residuals')
box()

color.cat.lag<-classIntervals(cov.lag.eig$residuals, n=5, style="quantile", dataPrecision=2)
collag<-findColours(color.cat.lag, colors)
plot(cov, col=collag)
title('Map of Lag Model Residuals')
legend('topleft', legend=c(names(attr(collag, 'table'))), fill=c(attr(collag, 'palette')), 
       title='Regression Residuals')
box()

writeOGR(cov, getwd(), "COVID19", "ESRI Shapefile")
writeOGR(cc, getwd(), "COVID19cc", "ESRI Shapefile")

#Examining source of spatial effects in lag model
effects.SAR<-impacts(cov.lag.eig, listw=covspat1.gal)
effects.SAR
