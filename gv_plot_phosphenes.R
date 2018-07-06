# Plots transparent elipses indicating individual phosphenes
# in each TMS session or average phosphene locations for all
# sessions or both
#
# GV 13/6/16
#
# Works by using a function that DHB wrote to compute finely
# sampled points in the shape of an elipse, the function 
# (getellipse) takes 5 arguments (detailed below). 
#
#

library(R.matlab)
save <- 1
plotting <- 2 # 1=plot individual phosphenes within a session; 2=plot mean phosphenes for all sessions; 3=plot both
tmsmode <- 4 # 1=both, 2=online only, 3=TBS only, 4=spTMS, rTMS(mean), cTBS and iTBS
subjects <- c('1','2','3','4','5','6')
alpha <- 0.5 # determines how transparent the ellipses are
colours <- c(rgb(0,0.4,0.6,alpha), rgb(0.8,0,0,alpha), rgb(0,0.6,0.4,alpha), rgb(0.35,0.25,0.5,alpha), rgb(0.15,0.7,0.7,alpha), rgb(0.8,0.5,0.6,alpha), rgb(0.9,0.6,0,alpha))
# colours: blue, red, green, purple, turquoise, pink, orange
plotwidthx <<- 1024
plotwidthy <<- 1024
res <- plotwidthy/plotwidthx
npixelsperdegree <- 26
crosssize <- 6

stimcircsize <- npixelsperdegree*4

wd <- getwd()

for (sub in 1:length(subjects)){
  
if (tmsmode==1) {
  conditions <- c('spTMS', 'rTMS1', 'rTMS2', 'rTMS3', 'rTMS4', 'coTBS', 'inTBS')
} else if (tmsmode==2){
  conditions <- c('spTMS', 'rTMS1', 'rTMS2', 'rTMS3', 'rTMS4')
} else if (tmsmode==3){
  conditions <- c('coTBS', 'inTBS')
} else if (tmsmode==4){
  conditions <- c('spTMS', 'rTMS', 'coTBS', 'inTBS')
  temp2 <- matrix(nrow=4,ncol=4)
  for (c in 1:4){
    data <- readMat(paste0(wd,"/phospheneLocations/S",subjects[sub],"rTMS",toString(c),".mat"))
    temp1 <- data$markedpositionandsize
    plast <- dim(temp1)[1] # first and last of the phosphenes that were actually used (the last five are used)
    pfirst <- plast - 5
    temp2[c,1:4] <- colMeans(temp1[c(pfirst:plast),])
  }
  rTMScoords <- colMeans(temp2)
}

getellipse <- function(p){
  theta <- -p[5] * (pi / 180)
  alpha <- 360*(0:1000)/1000 * (pi / 180)  
  ellipsematrix <- matrix(data=0,nrow=2,ncol=length(alpha))
  
  for (n in 1:length(alpha)){
    ellipsematrix[1,n] <- (plotwidthx/2) + p[3] + (1.18*p[1] * cos(alpha[n]) * cos(theta) - 1.18*p[2] * sin(alpha[n]) * sin(theta))
    ellipsematrix[2,n] <- (plotwidthy/2) - p[4] + (1.18*p[1] * cos(alpha[n]) * sin(theta) + 1.18*p[2] * sin(alpha[n]) * cos(theta))
  }
  return(ellipsematrix)}

if(save > 0){
  pdf(paste0(wd,"/phos_S",subjects[sub],".pdf"), height = 5, width = 5)
}

plot.new
par(pty='s')
plot(1,1,type='l',lty=1,lwd=2,axes=FALSE, ann=FALSE, xlim=c(0,plotwidthx), ylim=c(0,plotwidthy))

for (cond in 1:length(conditions)){
  if (cond == 2){
    coordinates <- rTMScoords
  } else {
    data <- readMat(paste0(wd,"/phospheneLocations/S",subjects[sub],conditions[cond],".mat"))
    b <- paste0(wd,"/phospheneLocations/S",subjects[sub],conditions[cond],".mat")
    coordinates <- data$markedpositionandsize
  }
  print(coordinates)
  if (plotting!=2){
    for (phosphene in c(1:5)){
      index <- (length(coordinates)/4)+1-phosphene # ridiculous calculation to index from the end of the coords matrix upwards
      coords <- coordinates[index,]
      ellipsematrix <- getellipse(c(coords[3], coords[4], coords[1]-plotwidthx/2, coords[2]-plotwidthy/2, 0))
      polygon(ellipsematrix[1,], ellipsematrix[2,], col=colours[cond], border=NA, lty=1, lwd=2)
    }
    }
  if (plotting!=1){
    indexes <- c((length(coordinates)-4):length(coordinates)) # indexes last five lines in the matrix
    plast <- dim(coordinates)[1] # first and last of the phosphenes that were actually used (the last five are used)
    pfirst <- plast - 5

    if (cond == 2){
      coords <- coordinates
    } else {
      coords <- colMeans(coordinates[c(pfirst:plast),])
    }
    ellipsematrix <- getellipse(c(coords[3]/2, coords[4]/2, coords[1]-plotwidthx/2, (coords[2]-384), 0))
    polygon(ellipsematrix[1,], ellipsematrix[2,], col=colours[cond], border=NA, lty=1, lwd=2)
    
    stimcircposx <- coords[1]
    stimcircposy <- coords[2]
    
    # Since matlab plots (0,0) in the top left corner and R plots (0,0) in the middle,
    # for the x axis we need to do (xpos - plotwidth/2)
    # for the y axis we need to do (ypos - plotheight/2)
    stimcircle <- getellipse(c(stimcircsize/2, stimcircsize/2, stimcircposx-plotwidthx/2, (stimcircposy-384), 0))
    polygon(stimcircle[1,], stimcircle[2,], col=NA, border=rgb(0,0,0,0.8), lty=1, lwd=1)

  }
  
}

locx <- 0
locy <- 0
visang <- getellipse(c(26*5, 26*5, locx, locy, 0))
polygon(visang[1,], visang[2,], col=NA, border=rgb(0.2,0.2,0.2,0.3), lty=1, lwd=1)
visang <- getellipse(c(26*10, 26*10, locx, locy, 0))
polygon(visang[1,], visang[2,], col=NA, border=rgb(0.2,0.2,0.2,0.3), lty=1, lwd=1)
visang <- getellipse(c(26*15, 26*15, locx, locy, 0))
polygon(visang[1,], visang[2,], col=NA, border=rgb(0.2,0.2,0.2,0.3), lty=1, lwd=1)


lines(c(plotwidthx/2-512,plotwidthx/2+512),c(plotwidthy/2,plotwidthy/2), col=rgb(0.2,0.2,0.2,0.3))
lines(c(plotwidthx/2,plotwidthx/2),c(plotwidthy/2-512,plotwidthy/2+512), col=rgb(0.2,0.2,0.2,0.3))

lines(c(plotwidthx/2-crosssize,plotwidthx/2+crosssize),c(plotwidthy/2,plotwidthy/2))
lines(c(plotwidthx/2,plotwidthx/2),c(plotwidthy/2-crosssize,plotwidthy/2+crosssize))

text(plotwidthx/2+20,(26*5+550),'5 deg',cex=1,adj=0,col=rgb(0.2,0.2,0.2,0.5))
text(plotwidthx/2+20,(26*10+575),'10 deg',cex=1,adj=0,col=rgb(0.2,0.2,0.2,0.5))
text(plotwidthx/2+20,(26*15+600),'15 deg',cex=1,adj=0,col=rgb(0.2,0.2,0.2,0.5))

if(save > 0){
  dev.off()
}

}
