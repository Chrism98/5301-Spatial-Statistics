# Setting working directory and loading necessary libraries
setwd("~/OSU/GEOG/5103/NYAIDS")
library(rgdal)
library(spdep)
library(RColorBrewer)
library(classInt)

#Importing shapefile and projecting it 
ny <- readOGR(dsn=getwd(), GDAL1_integer64_policy=T)
proj4string(ny)
proj4string(ny)<-CRS(" +proj=longlat +ellps=WGS84")
ny_NAD<-spTransform(ny, CRS("+init=epsg:3358"))
proj4string(ny_NAD)

#Creating a data frame out of the NYAIDS set
ny.df <- data.frame(ny_NAD)

#Picking a color palate for the map of AIDS rate in NYC
colors <- brewer.pal(5, "RdYlBu") #I chose a gradient from red to yellow to blue.
colcode <- findColours(color.cat.rate, colors)

#Plotting the data onto a map in a separate window, then giving it a title and legend
windows()
plot(ny, col=colcode)
title('AIDS Rate per 1,000 People in New York City')
legend('topleft', legend=c(names(attr(colcode, 'table'))),
       fill=c(attr(colcode, 'palate')), title='AIDS Rate by Quantile')

#Creating queen and rook neighbors
ny_nbq<-poly2nb(ny_NAD)
ny_nbr<-poly2nb(ny_NAD, queen=FALSE)

#Creating multiple nearest neighbor functions
coords<-coordinates(ny_NAD)
IDs<-row.names(as(ny_NAD, "data.frame"))
ny_kn1<-knn2nb(knearneigh(coords, k=1), row.names=IDs)
ny_kn2<-knn2nb(knearneigh(coords, k=2), row.names=IDs)
ny_kn4<-knn2nb(knearneigh(coords, k=4), row.names=IDs)

#Creating distance-based neighbors
dist<-unlist(nbdists(ny_kn1, coords))
summary(dist)
max_k1<-max(dist)

ny_kd1<-dnearneigh(coords, d1=0, d2=0.75*max_k1, row.names=IDs)
ny_kd2<-dnearneigh(coords, d1=0, d2=1*max_k1, row.names=IDs)
ny_kd3<-dnearneigh(coords, d1=0, d2=1.5*max_k1, row.names=IDs)

#Creating weights matracies
ny_nbq_w<- nb2listw(ny_nbq, zero.policy = TRUE)
ny_nbq_wb<-nb2listw(ny_nbq, style="B", zero.policy = TRUE)

ny_nbr_w<- nb2listw(ny_nbr, zero.policy = TRUE)
ny_kn1_w<- nb2listw(ny_kn1, zero.policy = TRUE)
ny_kn2_w<- nb2listw(ny_kn2, zero.policy = TRUE)
ny_kn4_w<- nb2listw(ny_kn4, zero.policy = TRUE)

#Running Global Moran's I tests for each weights matrix
moran.test(ny_NAD$AllCases, listw=ny_nbq_w, zero.policy = TRUE)
moran.test(ny_NAD$AllCases, listw=ny_nbr_w, zero.policy = TRUE)
moran.test(ny_NAD$AllCases, listw=ny_kn1_w, zero.policy = TRUE)
moran.test(ny_NAD$AllCases, listw=ny_kn2_w, zero.policy = TRUE)
moran.test(ny_NAD$AllCases, listw=ny_kn4_w, zero.policy = TRUE)

#Running LISA test for each weights matrix
cases <- order(ny_NAD$POSTAL)
nclocI <- localmoran(ny_NAD$AllCases, ny_nbq_w, alternative="two.sided", zero.policy = TRUE)
printCoefmat(data.frame(nclocI[cases,], row.names=ny_NAD$NAME[cases]), check.names=TRUE)
lmi<-data.frame(nclocI[cases], row.names=ny_NAD$FIPSNO[cases])

#Plotting points on a scatterplot
nci<- moran.plot(ny_NAD$AllCases, ny_kn2_w, labels=ny_NAD$AllCases, xlim=c(-0.5,2700), ylim=c(-0.5,1750), xlab="New York City AIDS Rate", ylab="Spatially Lagged New York City AIDS Rate")
title("Moran Scatterplot of AIDS Rates in New York City")
infl <- apply(nci$is.inf, 1, any)

#Creating a map based on LISA data in a separate window
x <- ny_NAD$AllCases
lhx <- cut(x, breaks=c(min(x), mean(x), max(x)), labels=c("L", "H"), include.lowest=TRUE)
wx <- lag(ny_kn2_w, ny_NAD$AllCases)
lhwx <- cut(wx, breaks=c(min(wx), mean(wx), max(wx)), labels=c("L", "H"), include.lowest=TRUE)
lhlh <- interaction(lhx, lhwx, infl, drop=TRUE)
cols <- rep(1, length(lhlh))
cols[lhlh == "H.L.TRUE"] <- 2
cols[lhlh == "L.H.TRUE"] <- 3
cols[lhlh == "H.H.TRUE"] <- 4
colors <- brewer.pal(5, "RdYlBu")
colcode <- findColours(color.cat.rate, colors)
windows()
plot(ny_NAD, col=colcode)
legend("topright", legend=c("None", "HL", "LH", "HH"), fill=colors, bty="n", cex=0.8, y.intersp=0.8)
title("LISA Map of Spatial Influence on AIDS Rate")

#Creating another scatterplot for outliers
Scrate <- (ny_NAD$AllCases-mean(ny_NAD$AllCases))/sd(ny_NAD$AllCases)
lagScrate <- lag.listw(ny_kn2_w, Scrate, zero.policy=T)
plot(Scrate, lagScrate, xlim=c(-4,4), ylim=c(-2,2),col="blue", xlab="AIDS Rate", 
     ylab="Spatially Lagged AIDS Rate", main="Moran Scatterplot of AIDS Rate in New York City")
abline(h=0, v=0)
abline(lm(lagScrate ~ Scrate), lty=2, lwd=2, col="red")

###Creating a LISA map with significant values###
#Creating standardized variables
ny.df$hh <- (Scrate>= 0 & lagScrate>= 0) & (nclocI[,5]<= 0.05)
ny.df$ll <- (Scrate<= 0 & lagScrate<= 0) & (nclocI[,5]<= 0.05)
ny.df$hl <- (Scrate>= 0 & lagScrate<= 0) & (nclocI[,5]<= 0.05)
ny.df$lh <- (Scrate<= 0 & lagScrate>= 0) & (nclocI[,5]<= 0.05)
ny.df$ns <- nclocI[,5]> 0.05

#Combining significant results into a single variable
ny.df$var <- 0
ny.df$var <- replace(ny.df$var, (ny.df$hh==TRUE), 1)
ny.df$var <- replace(ny.df$var, (ny.df$ll==TRUE), 2)
ny.df$var <- replace(ny.df$var, (ny.df$hl==TRUE), 3)
ny.df$var <- replace(ny.df$var, (ny.df$lh==TRUE), 4)
ny.df$var <- replace(ny.df$var, (ny.df$ns==TRUE), 5)
ny.df$var #Ensuring that all values are accounted for

#Setting breaks, labels, and colors for the classes in the LISA map
breaks <-seq(1,5,1)
labels <- c("High-High", "Low-Low", "High-Low", "Low-High", "Not Signif.")
colors <- c("red", "blue", "pink", "skyblue2", "white")

#Sorting into the correct category for the LISA map
np <- findInterval(ny.df$var, breaks, all.inside=FALSE)

#Creating the LISA map with the specified data and properties in a separate window
windows()
plot(ny)
plot(ny, col=colors[np], add=TRUE)
mtext("LISA Map of AIDS Rates in New York City by Quartile", cex=1.2, side=3, line=1)
legend("topleft", legend=labels, fill=colors, bty="n")

#Creating a graph showing high school graduation rates by quantile in a separate window
int.jenks <- classIntervals(ny$PctHSEd, n=5, style="jenks")
palette <- brewer.pal(5, "Blues")
col.j <- findColours(int.jenks, palette)
windows()
plot(ny,
     col=col.j,
     pch=19)
legend("topleft", fill = attr(col.j, "palette"),
       legend = names(attr(col.j, "table")), bty = "n")
mtext("Percentage of High School Graduates in New York City in Quantiles", cex=1.2, side=3, line=1)

#Creating a graph showing median household income by quantile in a separate window
int.jenks <- classIntervals(ny$MedInc, n=5, style="jenks")
palette <- brewer.pal(5, "Greens")
col.j <- findColours(int.jenks, palette)
windows()
plot(ny,
     col=col.j,
     pch=19)
legend("topleft", fill = attr(col.j, "palette"),
       legend = names(attr(col.j, "table")), bty = "n")
mtext("Median Home Income in New York City in Quantiles", cex=1.2, side=3, line=1)
