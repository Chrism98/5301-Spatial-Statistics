#loading libraries and setting working directory
library(spdep)
library(rgdal)
library(RColorBrewer)
library(car)
library(MASS)
library(robustbase)
library(olsrr)
library(classInt)
library(graphicsQC)
library(lmtest)
library(sp)
library(spatstat)
library(spatialreg)

clev<-readOGR(dsn="C:/Users/Chris/Documents/OSU/GEOG/5103/EJ_Cleveland")
clev.df<-data.frame(clev)
summary(clev)

cmin<-clev$MINORPCT
clow<-clev$LOWINCPCT
cvul<-clev$VULEOPCT

#finding correlations and plotting with minority percentage as the dependent variable
##VULEOPCT will be used later for model building
hist(clev$D_LDPNT_2)
cor(cmin,clev$D_LDPNT_2) #r=0.8162
lm(clev$D_LDPNT_2~cmin)
plot(cmin,clev$D_LDPNT_2,
     main="Minority Population vs Lead Paint Indicator", xlab="% Minority",
     ylab="Lead Paint Indicator", pch=20)
abline(-148.9,434.5, col="red") #all abline values obtained from linear model (example line 28)

hist(clev$PM25)
cor(cmin,clev$PM25) #r=0.4387
lm(clev$PM25~cmin)
plot(cmin,clev$PM25,
     main="Minority Population vs Particulate Matter Levels", xlab="% Minority",
     ylab="Particulate Matter Levels", pch=20)
abline(9.4325,0.2094, col="red")

hist(clev$OZONE)
cor(cmin,clev$OZONE) #r=0.2214
hist(log(clev$OZONE))
cor(cmin,log(clev$OZONE)) #r=0.2241

hist(clev$PTSDF)
hist(log(clev$PTSDF))
cor(cmin,clev$PTSDF) #r=0.6358
cor(cmin,(log(clev$PTSDF))) #r=0.6516
lm(log(clev$PTSDF)~cmin)
plot(cmin,log(clev$PTSDF),
     main="Minority Population vs Proximity to TSDF Facilities", xlab="% Minority",
     ylab="log(Proximity to TSDF Facilities)", pch=20)
abline(0.1674,2.1737, col="red")    

hist(clev$PRMP)    
hist(log(clev$PRMP))
cor(cmin,clev$PRMP) #r=0.2338
cor(cmin,log(clev$PRMP)) #r=0.3574
lm(log(clev$PRMP)~cmin)
plot(cmin,log(clev$PRMP),
     main="Minority Population vs Count of Risk Management Plans", xlab="% Minority",
     ylab="log(Count of Risk Management Plans)", pch=20)
abline(-0.8853,1.0805, col="red")    

hist(clev$PNPL)
hist(log(clev$PNPL))
cor(cmin,clev$PNPL) #r=-0.1303
cor(cmin,log(clev$PNPL)) #r=-0.1054

hist(clev$PTRAF)
cor(cmin,clev$PTRAF) #r=0.066

hist(clev$RESP)
hist(log(clev$RESP))
cor(cmin,clev$RESP) #r=0.5915
cor(cmin,log(clev$RESP)) #r=0.6004
lm(log(clev$RESP)~cmin)
plot(cmin,log(clev$RESP),
     main="Minority Population vs log(Air Toxics Respiratory Hazard Index)", xlab="% Minority",
     ylab="log(Air Toxics Respiratory Hazard Index)", pch=20)
abline(-1.1256,0.1796, col="red") 

#finding correlations and plotting with low income percentage as the dependent variable
hist(clev$D_LDPNT_2)
cor(clow,clev$D_LDPNT_2) #r=0.7698
lm(clev$D_LDPNT_2~clow)
plot(clow,clev$D_LDPNT_2,
     main="Low Income Population vs Lead Paint Indicator", xlab="% Minority",
     ylab="Lead Paint Indicator", pch=20)
abline(-178.0,569.9, col="red")

hist(clev$PM25)
cor(clow,clev$PM25) #r=0.4556
lm(clev$PM25~clow)
plot(clow,clev$PM25,
     main="Low Income Population vs Particulate Matter Levels", xlab="% Minority",
     ylab="Particulate Matter Levels", pch=20)
abline(9.4072,0.3024, col="red")

hist(clev$OZONE)
cor(clow,clev$OZONE) #r=0.0836
hist(log(clev$OZONE))
cor(clow,log(clev$OZONE)) #r=0.08639

hist(clev$PTSDF)
hist(log(clev$PTSDF))
cor(clow,clev$PTSDF) #r=0.7193
cor(clow,(log(clev$PTSDF))) #r=0.6748
lm(clev$PTSDF~clow)
plot(clow,clev$PTSDF,
     main="Low Income Population vs Proximity to TSDF Facilities", xlab="% Minority",
     ylab="Proximity to TSDF Facilities", pch=20)
abline(0.4895,11.5384, col="red")    

hist(clev$PRMP)    
hist(log(clev$PRMP))
cor(clow,clev$PRMP) #r=0.4608
cor(clow,log(clev$PRMP)) #r=0.5365
lm(log(clev$PRMP)~clow)
plot(clow,log(clev$PRMP),
     main="Low Income Population vs Count of Risk Management Plans", xlab="% Minority",
     ylab="log(Count of Risk Management Plans)", pch=20)
abline(-1.297,2.256, col="red")    

hist(clev$PNPL)
hist(log(clev$PNPL))
cor(clow,clev$PNPL) #r=-0.2649
cor(clow,log(clev$PNPL)) #r=-0.2482


hist(clev$PTRAF)
cor(clow,clev$PTRAF) #r=0.1545

hist(clev$RESP)
hist(log(clev$RESP))
cor(clow,clev$RESP) #r=0.6878
cor(clow,log(clev$RESP)) #r=0.6959
lm(log(clev$RESP)~clow)
plot(clow,log(clev$RESP),
     main="Low Income Population vs log(Air Toxics Respiratory Hazard Index)", xlab="% Minority",
     ylab="log(Air Toxics Respiratory Hazard Index)", pch=20)
abline(-1.1595,0.2895, col="red")

# finding correlations between VULEOPCT and independent variables significant with both minority % and low income %
cor(cvul, clev$RESP) #r=0.6803
cor(cvul, log(clev$RESP)) # r=0.6895

cor(cvul, clev$PRMP) #r=0.3540
cor(cvul,log(clev$PRMP)) #r=0.4655

cor(cvul, clev$PTSDF) #r=0.7222
cor(cvul,log(clev$PTSDF)) #r=0.7121

cor(cvul, clev$PM25) #r=0.4800

cor(cvul, clev$D_LDPNT_2) #r=0.8580

hist(cvul, col="Gray", prob=T,
     main = "Distribution of VULEOPCT",
     xlab = "VULEOPCT", ylab = "Density")
lines(density(cvul, bw=.05), col="red", lwd=2,)

#model building based on correlations, using VULEOPCT as the dependent variable
clev1<-lm(cvul~D_LDPNT_2+PM25+PTSDF+log(PRMP)+log(RESP), data=clev)
summary(clev1)

clev2<-lm(cvul~D_LDPNT_2+log(PM25)+PTSDF+log(PRMP)+log(RESP), data=clev)
summary(clev2)

clev3<-lm(cvul~D_LDPNT_2+log(PM25)+log(PTSDF)+log(PRMP)+log(RESP), data=clev)
summary(clev3) #best model so far

#checking for high VIF values
vif(clev3) #some concern about log(PTSDF), but workable

#running a robust regression
clev.rob<-lmrob(cvul~D_LDPNT_2+log(PM25)+log(PTSDF)+log(PRMP)+log(RESP), data=clev)
confint(clev.rob)
confint(clev3)
plot(clev.rob$fitted.values, clev.rob$residuals, main="Residuals vs. Fitted Values",
     xlab="Fitted Values", ylab="Residuals", pch=21); abline(h=0, col="red")

#performing OLS regression for outliers
windows()
ols_plot_resid_stand(clev3)
windows()
ols_plot_resid_stud(clev3)
windows()
ols_plot_dfbetas(clev3)
windows()
ols_plot_cooksd_chart(clev3)

#assigning weights to run Moran's I tests on data
clev_NAD<-spTransform(clev, CRS("+init=epsg:3358"))
proj4string(clev_NAD)

##queen and rook
clev_nbq<-poly2nb(clev_NAD)
clev_nbr<-poly2nb(clev_NAD, queen=FALSE)

##nearest neighbors
coords<-coordinates(clev_NAD)
IDs<-row.names(as(clev_NAD, "data.frame"))
clev_kn1<-knn2nb(knearneigh(coords, k=1), row.names=IDs)
clev_kn2<-knn2nb(knearneigh(coords, k=2), row.names=IDs)
clev_kn4<-knn2nb(knearneigh(coords, k=4), row.names=IDs)

clev_kn1_w<- nb2listw(clev_kn1, zero.policy = TRUE)
clev_kn2_w<- nb2listw(clev_kn2, zero.policy = TRUE)
clev_kn4_w<- nb2listw(clev_kn4, zero.policy = TRUE)

##creating weights matricies
clev_nbq_w<- nb2listw(clev_nbq, zero.policy = TRUE)
clev_nbq_wb<-nb2listw(clev_nbq, style="B", zero.policy = TRUE)

clev_nbr_w<- nb2listw(clev_nbr, zero.policy = TRUE)
clev_kn1_w<- nb2listw(clev_kn1, zero.policy = TRUE)
clev_kn2_w<- nb2listw(clev_kn2, zero.policy = TRUE)
clev_kn4_w<- nb2listw(clev_kn4, zero.policy = TRUE)

moran.test(clev_NAD$VULEOPCT, listw=clev_nbq_w, zero.policy = TRUE) #0.8396
moran.test(clev_NAD$VULEOPCT, listw=clev_nbr_w, zero.policy = TRUE) #0.8450
moran.test(clev_NAD$VULEOPCT, listw=clev_kn1_w, zero.policy = TRUE) #0.8680
moran.test(clev_NAD$VULEOPCT, listw=clev_kn2_w, zero.policy = TRUE) #0.8549
moran.test(clev_NAD$VULEOPCT, listw=clev_kn4_w, zero.policy = TRUE) #0.8499

#performing spatial regression due to presence of outliers
##searching for areas with zero neighbors
clevtest<- poly2nb(clev)
lapply(clevtest, function(x) attr(clevtest, "region.id")[x]) #none found

clevspat1.gal<-nb2listw(poly2nb(clev))
summary(clevspat1.gal)

#creating and attaching new data set
attach(clev.df)
class(clev.df)
names(clev.df)
row.names(clev.df)<-clev.df$COMM

lvul<-logit(cvul)

#running moran and geary tests on residuals from above regression
moran.test(clev3$residuals, clevspat1.gal, alternative="two.sided")
geary(clev3$residuals,clevspat1.gal,n=length(lvul), n1=length(lvul)-1,
      S0=length(lvul))

#running lagrange multiplier tests
lm.LMtests(clev3,clevspat1.gal, test=c("LMerr", "LMlag", "RLMerr", "RLMlag", "SARMA"))

#LMlag model was the highest significant value, it will be used going forward. 
##extracting and returning eigenvalues of the weights matrix
w.eig<-eigenw(clevspat1.gal,quiet=NULL)
1/min(w.eig)
1/max(w.eig)

clev.lag.eig<-lagsarlm(lvul~D_LDPNT_2+log(PM25)+log(PTSDF)+log(PRMP)+log(RESP), data=clev,
                       clevspat1.gal, method="eigen",quiet=FALSE)
summary(clev.lag.eig)

#looking at residuals
me1<-mean(residuals(clev.lag.eig))
me1 #effectively zero
sd1<-sd(residuals(clev.lag.eig))
sd1 #0.4737
summary(residuals(clev.lag.eig))
plot(residuals(clev.lag.eig)); abline(h=0, col="red") #reasonably symetric

#building a histogram to visualize residual data
hist(residuals(clev.lag.eig), breaks=seq(-3, 3, .2),
     col=8, probability=T,
     ylab='Density',
     main='Histogram of Residuals(Spatial Lag Model)',
     xlab='Residuals(Spatial Lag Model)')
box()
curve(dnorm(x, mean=me1, sd=sd1), from=-3, to=3, add=T,
      col='red', lwd=2)

#looking at outliers
(sd(residuals(clev.lag.eig)))*3 #1.421
which(residuals(clev.lag.eig)<=-1.421)
which(residuals(clev.lag.eig)>=1.421)

#estimating spatial lag model using 2SLS estimator
clev2SLS<-stsls(lvul~D_LDPNT_2+log(PM25)+log(PTSDF)+log(PRMP)+log(RESP), 
                data=clev, clevspat1.gal)
summary(clev2SLS)

clev2SLSR<-stsls(lvul~D_LDPNT_2+log(PM25)+log(PTSDF)+log(PRMP)+log(RESP), 
                 robust=T, data=clev, clevspat1.gal)
summary(clev2SLSR)

#spatial eror model
clev.err.eig<-errorsarlm(lvul~D_LDPNT_2+log(PM25)+log(PTSDF)+log(PRMP)+log(RESP), data=clev,
                         clevspat1.gal, method="eigen",quiet=FALSE)
summary(clev.err.eig)

#looking at residuals of spatial error model
me2<-mean(residuals(clev.err.eig))
me2 #effectively zero
sd2<-sd(residuals(clev.err.eig))
sd2 #0.4923
summary(residuals(clev.err.eig))
plot(residuals(clev.err.eig)); abline(h=0, col="red") #reasonably symetric

#building a histogram to visualize residual data
hist(residuals(clev.err.eig), breaks=seq(-3, 3, .2),
     col=8, probability=T,
     ylab='Density',
     main='Histogram of Residuals(Spatial Error Model)',
     xlab='Residuals(Spatial Error Model)')
box()
curve(dnorm(x, mean=me2, sd=sd2), from=-3, to=3, add=T,
      col='red', lwd=2)

#checking residuals for spatial autocorrelation using Q1 matrix
names(clev.err.eig)
moran.test(clev.err.eig$residuals, clevspat1.gal, alternative="two.sided")

#estimating error model using KP method
clevKP<-GMerrorsar(lvul~D_LDPNT_2+log(PM25)+log(PTSDF)+log(PRMP)+log(RESP), data=clev,
                   clevspat1.gal, verbose=FALSE)
summary(clevKP)

#comparing autocorrelation in residuals
moran.test(clev3$residuals, clevspat1.gal, alternative="two.sided") #I=0.2952
moran.test(clev.lag.eig$residuals, clevspat1.gal, alternative="two.sided") #I=-0.0385
moran.test(clev.err.eig$residuals, clevspat1.gal, alternative="two.sided") #I=-0.0475
moran.test(clevKP$residuals, clevspat1.gal, alternative="two.sided") #I=0.3542
moran.test(clev2SLS$residuals, clevspat1.gal, alternative="two.sided") #I=-0.0562
moran.test(clev2SLSR$residuals, clevspat1.gal, alternative="two.sided") #I=-0.0562

#comparing AIC values
AIC(clev3) #-1731.703
AIC(clev.lag.eig) #1630.637
AIC(clev.err.eig) #1759.749

#plotting residuals across the different models
windows()
par(mfrow=c(1,2))
par(mar=c(1.2,1.2,1.2,1.2))
colors<-brewer.pal(5, "YlOrBr")
color.cat.lag<-classIntervals(clev.lag.eig$residuals, n=5, style="quantile", dataPrecision=2)
collag<-findColours(color.cat.lag, colors)
plot(clev, col=collag)
title('Map of Lag Model Residuals')
legend('topleft', legend=c(names(attr(collag, 'table'))), fill=c(attr(collag, 'palette')), 
       title='Regression Residuals')
box()

color.cat.err<-classIntervals(clev.err.eig$residuals, n=5, style="quantile", dataPrecision=2)
colerr <- findColours(color.cat.err, colors)
plot(clev, col=colerr)
title('Map of Error Model Residuals')
legend('topleft', legend=c(names(attr(colerr, 'table'))), fill=c(attr(colerr, 'palette')), 
       title='Regression Residuals')
box()

clev$lvul<-logit(clev$VULEOPCT)
writeOGR(clev, getwd(), "EJ_Cleveland2", "ESRI Shapefile")

#Examining source of spatial effects in lag model
effects.SAR<-impacts(clev.lag.eig, listw=clevspat1.gal)
effects.SAR
