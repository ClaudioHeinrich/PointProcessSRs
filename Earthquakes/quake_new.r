### Earthquakes from Entiat, Leavenworth, and Grand Coolee area
### 1969 through Oct 1986

### x has seven columns:
### time (in days from 01011970) x y depth mag yrmodyhrmn sec

x <- read.table("qamar")

####################################################################
### Plot locations on a WA map

library(maps)
positions <- cbind(x[,2],x[,3])
positions[,1] <- -(121+positions[,1]/111.181)
positions[,2] <- 47.5+positions[,2]/75.344

k <- 0.05
x1 <- min(positions[,1])-k
x2 <- max(positions[,1])+k
y1 <- min(positions[,2])-k
y2 <- max(positions[,2])+k

#pdf(file="images/wa_loc_mag.pdf",width=3,height=4.8)
par(mfrow=c(1,1),mex=0.1)
map('state',region=c('washington'),
    xlim=c(x1,x2),
    ylim=c(y1,y2),
    fill=TRUE,col="grey",mar=c(0,0,0,0),lwd=0.5)
points(positions[,1],y=positions[,2],pch=1,cex=0.2*(x[,5]+0.9))
dev.off()

pdf(file="images/wa_loc_g2.pdf",width=3,height=4.8)
par(mfrow=c(1,1),mex=0.1)
map('state',region=c('washington'),
    xlim=c(x1,x2),
    ylim=c(y1,y2),
    fill=TRUE,col="grey",mar=c(0,0,0,0),lwd=0.5)
points(positions[which(x[,5]>1.9),1],y=positions[which(x[,5]>1.9),2],pch=16,cex=0.2)
dev.off()

pdf(file="images/wa.pdf",width=3,height=2)
par(mfrow=c(1,1),mex=0.1)
map('state',region=c('washington'),
    fill=TRUE,col="grey",mar=c(0,0,0,0),lwd=0.5)
lines(c(x1,x1),c(y1,y2),lwd=0.7)
lines(c(x2,x2),c(y1,y2),lwd=0.7)
lines(c(x1,x2),c(y1,y1),lwd=0.7)
lines(c(x1,x2),c(y2,y2),lwd=0.7)
dev.off()
#####################################################################


y <- read.table("cal.data")

####################################################################
### Plot locations on a CA map

positions <- cbind(-y[,8],y[,7])

k <- 0.05
x1 <- min(positions[,1])-k
x2 <- max(positions[,1])+k
y1 <- min(positions[,2])-k
y2 <- max(positions[,2])+k

pdf(file="images/ca_loc_mag.pdf",width=3.3,height=3)
par(mfrow=c(1,1),mex=0.1)
map('state',region=c('california'),
    xlim=c(x1,x2),
    ylim=c(y1,y2),
    fill=TRUE,col="grey",mar=c(0,0,0,0),lwd=0.5)
points(positions[,1],y=positions[,2],pch=1,cex=0.2*y[,9])
dev.off()

pdf(file="images/ca_loc_g3.pdf",width=3.3,height=3)
par(mfrow=c(1,1),mex=0.1)
map('state',region=c('california'),
    xlim=c(x1,x2),
    ylim=c(y1,y2),
    fill=TRUE,col="grey",mar=c(0,0,0,0),lwd=0.5)
points(positions[which(y[,9]>2.9),1],y=positions[which(y[,9]>2.9),2],pch=16,cex=0.2)
dev.off()


pdf(file="images/ca.pdf",width=2,height=2)
par(mfrow=c(1,1),mex=0.1)
map('state',region=c('california'),
    fill=TRUE,col="grey",mar=c(0,0,0,0),lwd=0.5)
lines(c(x1,x1),c(y1,y2),lwd=0.7)
lines(c(x2,x2),c(y1,y2),lwd=0.7)
lines(c(x1,x2),c(y1,y1),lwd=0.7)
lines(c(x1,x2),c(y2,y2),lwd=0.7)
dev.off()
#####################################################################

### Get times for CA earthquakes, in days from 1/1/1970
y.dates <- paste(y[,3],y[,2],sep="/")
y.dates <- paste(y.dates,y[,1],sep="/19")
y.times <- paste(y[,4],y[,5],y[,6],sep=":")
y.time <- paste(y.dates,y.times)
y.time <- strptime(y.time,"%d/%m/%Y %H:%M:%S")
y.julian <- julian(y.time)
ind <- c(1:2049) # which(y[,9]>3.9)
y.julian <- y.julian[ind]

dates <- c("1/1/1962","1/1/1963","1/1/1964","1/1/1965","1/1/1966",
           "1/1/1967","1/1/1968","1/1/1969","1/1/1970","1/1/1971",
           "1/1/1972","1/1/1973","1/1/1974","1/1/1975","1/1/1976",
           "1/1/1977","1/1/1978","1/1/1979","1/1/1980","1/1/1981",
           "1/1/1982")
dates <- strptime(dates,"%d/%m/%Y")
dates.julian <- julian(dates)

pdf(file="images/ca_time_mag.pdf",width=7,height=4,points=12)
par(mfrow=c(4,1),mex=0.6,mar=c(3,2,1,2)+0.01)
k1 <- min(which(y[ind,1]==67))
k2 <- min(which(y[ind,1]==72))
k3 <- min(which(y[ind,1]==77))
k4 <- length(y[ind,1])
plot(y.julian[1:(k1-1)],rep(1,k1-1),type="p",pch="|",yaxs="r",
     cex=y[1:(k1-1),9]/2,xlab="",ylab="",axes=FALSE)
axis(1,at=dates.julian[1:6],labels=dates[1:6])
plot(y.julian[k1:(k2-1)],rep(1,k2-k1),type="p",pch="|",yaxs="r",
     cex=y[k1:(k2-1),9]/2,xlab="",ylab="",axes=FALSE)
axis(1,at=dates.julian[6:11],labels=dates[6:11])
plot(y.julian[k2:(k3-1)],rep(1,k3-k2),type="p",pch="|",yaxs="r",
     cex=y[k2:(k3-1),9]/2,xlab="",ylab="",axes=FALSE)
axis(1,at=dates.julian[11:16],labels=dates[11:16])
plot(y.julian[k3:k4],rep(1,k4-k3+1),type="p",pch="|",yaxs="r",
     cex=y[k3:k4,9]/2,xlab="",ylab="",axes=FALSE)
axis(1,at=dates.julian[16:21],labels=dates[16:21])
dev.off()
#################################################################

### Plot times for WA earthquakes
ind <- c(1:length(x[,1])) #which(x[,5]>1.9)
x.julian <- x[ind,1]

dates <- c("1/1/1969","1/1/1970","1/1/1971",
           "1/1/1972","1/1/1973","1/1/1974","1/1/1975","1/1/1976",
           "1/1/1977","1/1/1978","1/1/1979","1/1/1980","1/1/1981",
           "1/1/1982","1/1/1983","1/1/1984","1/1/1985","1/1/1986",
           "1/1/1987","1/1/1988","1/1/1989")
dates <- strptime(dates,"%d/%m/%Y")
dates.julian <- julian(dates)

x.mag <- x[ind,5] + 0.9

pdf(file="images/wa_time_mag.pdf",width=7,height=4,points=12)
par(mfrow=c(4,1),mex=0.6,mar=c(3,2,1,2)+0.01)
##par(mfrow=c(4,1),mex=0.7)
tmp <- as.numeric(substr(x[,6],start=1,stop=2))
k1 <- min(which(tmp[ind]>=74))
k2 <- min(which(tmp[ind]>=79))
k3 <- min(which(tmp[ind]>=84))
k4 <- 6941
plot(x.julian[1:(k1-1)],rep(1,k1-1),type="p",pch="|",yaxs="r",
     xlab="",ylab="",axes=FALSE,xlim=c(dates.julian[1],dates.julian[6]),
     cex=x.mag[1:(k1-1)]/2)
axis(1,at=dates.julian[1:6],labels=dates[1:6])
plot(x.julian[k1:(k2-1)],rep(1,k2-k1),type="p",pch="|",yaxs="r",
     xlab="",ylab="",axes=FALSE,xlim=c(dates.julian[6],dates.julian[11]),
     cex=x.mag[k1:(k2-1)]/2)
axis(1,at=dates.julian[6:11],labels=dates[6:11])
plot(x.julian[k2:(k3-1)],rep(1,k3-k2),type="p",pch="|",yaxs="r",
     xlab="",ylab="",axes=FALSE,xlim=c(dates.julian[11],dates.julian[16]),
     cex=x.mag[k2:(k3-1)]/2)
axis(1,at=dates.julian[11:16],labels=dates[11:16])
plot(x.julian[k3:k4],rep(1,k4-k3+1),type="p",pch="|",yaxs="r",
     xlab="",ylab="",axes=FALSE,,xlim=c(dates.julian[16],dates.julian[21]),
     cex=x.mag[k3:k4]/2)
axis(1,at=dates.julian[16:21],labels=dates[16:21])
dev.off()
