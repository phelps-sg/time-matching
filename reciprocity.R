
# Authors: Wing Lon Ng and Steve Phelps
# Released: 23/3/2016
# 
# R source code for analysing the data in Phelps S, Ng WL, Musolesi M, Russell YI. 
# Data from: Chimpanzees exit the market: evidence for immediate but not delayed 
# time matching in social grooming. Dryad Digital Repository. 
# http://dx.doi.org/10.5061/dryad.1rn27
# 
# To the extent possible under law, the authors have dedicated all copyright and 
# related and neighboring rights to this software to the public domain worldwide. 
# This software is distributed without any warranty.
# 
# You should have received a copy of the CC0 Public Domain Dedication along with 
# this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>
# 

library(xtable)
library(sets)
library(lme4)
library(ggplot2)

datadir <- '../../data/'
figsdir <- '../../figs/'

# Raw data
d <- read.table(paste0(datadir,"raw_grooming_data.txt"), 
		  colClasses = "character", sep=",", na.strings = "NA",
		  blank.lines.skip = T, fill=T, col.names=1:4)
d[,2] <- sub(" ", "", d[,2])
d[,3] <- sub(" ", "", d[,3])
d[,4] <- sub(" ", "", d[,4])

idx_schedule <- which(nchar(d[,4])==0)
schedule <- d[idx_schedule,]
nschedules <- nrow(schedule)

session <- matrix(0,nrow(d),1)
session[idx_schedule,1] <- 1
session <- cumsum(session)
d <- cbind(session, d)
d <- d[-idx_schedule,]

chimp <- as.factor(sort(unique(c(d$X1,d$X2))))
nchmp <- length(chimp)
ndyad <- nchmp*(nchmp-1)/2
d$X1  <- factor(d$X1,levels=chimp)
d$X2  <- factor(d$X2,levels=chimp)
adjmat<- matrix(0,nchmp,nchmp)
rownames(adjmat) <- chimp
colnames(adjmat) <- chimp
re1 <- adjmat
re2 <- adjmat

adjmat2 <- matrix(NA,nchmp,nchmp)
rownames(adjmat2)<-chimp
colnames(adjmat2)<-chimp
adjmat2[upper.tri(adjmat2)] <- 0

comov <- read.csv(paste0(datadir, "comovement.csv"))
como2 <- cbind(comov[,26],comov[,-26])
#como2 <- as.matrix(como2)
#como2[is.nan(como2)] <-0
como2[is.na(como2)] <- NaN
#como2 <- como2+t(como2)
rownames(como2) <- rownames(adjmat)
colnames(como2) <- colnames(adjmat)

proxim <- read.csv(paste0(datadir, "proxim.csv"))
proxim <- proxim[,-1]
rownames(proxim) <- rownames(adjmat)
colnames(proxim) <- colnames(adjmat)

kin <- read.csv(paste0(datadir,"kin.csv"))
kin <- kin[,-1]
rownames(kin) <- rownames(kin)

#[1]       AL   BO   CL CR DL  FA  FR HI HO  HP JT KI  KK   KT   KY LA   LI MA NI   PA RO   SA SR   WH WI ZZ
age <- matrix(c(11.5,37.5,32.5, 6,16,27.5,27,31,7.5,28, 3,7.5,19.5,37.5,2 ,12.5,8 ,25,34.5, 6,30.5,15,16.5,10,35,9.5),26,1)

age2 <- age%*%t(rep(1,26))
age2 <- (age2+t(age2))/2
all(age == diag(age2))

for(j in 1:nrow(d)){
  re2[as.numeric(d$X1[j]),as.numeric(d$X2[j])] <- 
      re2[as.numeric(d$X1[j]), as.numeric(d$X2[j])] + 1
}
res <- re2+t(re2)
re1[upper.tri(re1)] <- 1:ndyad
res[lower.tri(res)] <- NA
diag(res) <- NA

window_ws <- round(2 ^ seq(-1, 7, by=0.5) * 60)
res_dyad  <- vector("list", ndyad)
dis_wind  <- vector("list", length(window_ws))
res_wind  <- matrix(NA, length(window_ws), 2)

for(i in 1:nschedules){
  now1 <- Sys.time()
  data <- d[d$session==i,]
  
  t1 <- matrix(as.numeric(unlist(strsplit(data$X3,":"))), nrow(data), 2, byrow=T)
  t1 <- t1[,1]*60 + t1[,2]
  t2 <- matrix(as.numeric(unlist(strsplit(data$X4,":"))), nrow(data), 2, byrow=T)
  t2 <- t2[,1]*60 + t2[,2]
  dur<- t2 - t1
  idx <- which(dur < 0)
  dur[idx] <- dur[idx] + 60*60
  data <- cbind(data, dur)
  
  idx <- which(diff(t1)<0)
  if (length(idx) > 0) {
    for(j in 1:length(idx)) {
      tmp <- (idx[j]+1):min(idx[j+1], nrow(data), na.rm=T)
      t1[tmp] <- t1[tmp]+j*3600
    }
  }
  t3 <- t1+dur
  data<-cbind(data, t1, t3)
  tmax <- max(t3)
  
  x <- cbind(rep(i,tmax), 1:tmax, rep(NA,tmax))
    
  if (i==1) {
    d2 <- data
    timegrids <- x
  } else {
    d2 <- rbind(d2, data)
    timegrids <- rbind(timegrids, x)
  }
}

# Up till this point d2 contains single events with normalised times

JoinEvents<- function(d2) {
  # Pair subsequent events from the same dyad and same session into a single row
  # so that the final result is a data-frame containing rows of dy-events.
  now1 <- Sys.time()
  tmp <- NULL
  for (k in 1:ndyad) {
  #  k=1
    idx_pair<- which(re1==k, arr.ind=T)
    chimp1  <- chimp[idx_pair[1]]
    chimp2  <- chimp[idx_pair[2]]
    idx_pos <- which((d$X1==chimp1 & d$X2==chimp2)|(d$X1==chimp2 & d$X2==chimp1))
    idx_ses <- unique(d$session[idx_pos])
    data    <- d2[idx_pos,]
  
    if (nrow(data) < 2) next
  
    ############################
    idx_session <- (table(data$session)>=2)
    idx_session <- as.integer(names(idx_session)[idx_session])
  
    if (length(idx_session) == 0) next
    
    for (i in 1:length(idx_session)) {
      dta_sel <- data[data$session==idx_session[i],]
      for (j in 2:nrow(dta_sel)) {
        if (dta_sel[j,2] == dta_sel[j-1,3]) {
          # Join reciprocation events
          tmp <- rbind(tmp, cbind(k, dta_sel[j-1,], dta_sel[j,]))
        }
      }
    }
  }
  colnames(tmp) <- paste0(colnames(tmp), ".", rep(c(1,2), each=9))[-18]
  return(tmp)
}

dataset <- JoinEvents(d2)

Y <- dataset$t3.2 - dataset$t1.2
X1<- dataset$t3.1 - dataset$t1.1
X2<- dataset$t1.2 - dataset$t3.1
X3<- dataset$t1.2 - dataset$t1.1
X4<- dataset$t3.2 - dataset$t3.1
X5<- dataset$t3.2 - dataset$t1.1
age_X1  <- age[dataset$X1.1]
age_Y   <- age[dataset$X2.1]
Z <- data.frame(Y,X1,X2,X3,X4)/60
ind<- sign(X2)
kin_idx <- rep(NA,nrow(Z))
com_idx <- kin_idx
proxidx <- kin_idx 

for (i in 1:nrow(Z)) {
  idx1 <- as.numeric(dataset$X1.1[i])
  idx2 <- as.numeric(dataset$X2.1[i])
  kin_idx[i] <- kin[idx1,idx2]
  com_idx[i] <- como2[idx1,idx2]
  proxidx[i] <- proxim[idx1,idx2]
}

Z <- cbind(Z,ind,kin_idx,com_idx,proxidx,age_X1,age_Y)

# Duration of second event
dataset$Y <- dataset$t3.2-dataset$t1.2

# Duration of first event
dataset$X1 <- dataset$t3.1 - dataset$t1.1

# Time between end of first event and start of second
dataset$X2 <- dataset$t1.2 - dataset$t3.1

# Time between start of first event and start of second
#  which is always positive.
dataset$X3 <- dataset$t1.2 - dataset$t1.1

# Time between end of first event and end of second event
dataset$X4 <- dataset$t3.2 - dataset$t3.1

# Time between start of first event and end of second event
dataset$X5 <- dataset$t3.2 - dataset$t1.1

Dyad <- function(i, j) {
  as.list(sort(c(i, j)))
}

Comovement <- function(i) {
  return(como2[dyads[[i,1]], dyads[[i,2]]])
}

DyadNumber <- function(i, j) {
  d <- Dyad(i, j)
  for (i in 1:nrow(dyads)) {
    if (unlist(dyads[i, 1]) == d[[1]] && unlist(dyads[i, 2]) == d[[2]]) {
      return(i)
    }
  }
}

dyads <- c()
for (i in chimp) {
  for (j in chimp) {
    if (i != j) {
      dyads <- rbind(dyads, Dyad(i, j))
    }
  } 
}
dyads <- unique(dyads)

dym <- as.matrix(dataset[,c('X1.1', 'X1.2')])
dataset$dyad <- apply(dym, 1, function(x) { DyadNumber(x[[1]], x[[2]]) })
dataset$comovement <- sapply(dataset$dyad, Comovement)
dataset$reciprocity = dataset$Y / (dataset$X1 + dataset$Y)


CumulativeReciprocity <- function(dataset, ws, treatment='all') {
  dataset$window <- floor(dataset$t1.1 / ws)
  SY <- with(dataset,  aggregate(Y ~ session.1 + dyad + window,  FUN=sum))
  SX1 <- with(dataset,  aggregate(X1 ~ session.1 + dyad + window,  FUN=sum))
  result <- merge(SY, SX1)
  result$R <- log(result$Y / result$X)
#  names(result) = c('session', 'dyad', 'window', 'R')
  result$treatment = replicate(nrow(result), treatment)
  result$ws <- replicate(nrow(result), ws)
  return(result)
}

CumulativeReciprocityByDelay <- function(dataset, ws) {
  reciprocity.all <- CumulativeReciprocity(dataset, ws, 'all')
  reciprocity.delayed <- CumulativeReciprocity(dataset[dataset$X2 >= 0,], ws, 'delayed')
  reciprocity.immediate <- CumulativeReciprocity(dataset[dataset$X3 == 0,], ws, 'immediate') 
  reciprocity.intra <- CumulativeReciprocity(dataset[dataset$X2 < 0,], ws, 'intra')
  return( 
    Reduce(rbind, 
      list(reciprocity.all, reciprocity.delayed, reciprocity.immediate, 
	    reciprocity.intra))
  )
}

CumulativeReciprocityByWindowSize <- 
		      function(dataset, 
				window.sizes = c(20, 40, 60, 240) * 60) {
  
  return(
    Reduce(rbind, 
	    Map(function(ws) CumulativeReciprocityByDelay(dataset, ws), window.sizes)))
}

MakeBoxPlots <- function(data, suffix='', color='blue') {

  FName <- function(name) {
    if (suffix != '') {
      fname.suffix = paste0('-', suffix)
    } else {
      fname.suffix = ''
    }
    return(paste0(figsdir, name, fname.suffix, '.pdf'))
  }
  
  Title <- function(main) {
      return('')
  }
  
  Plot <- function(formula, fname, title, label=expression(rho)) {
    pdf(FName(fname))
    #par(mar=c(0, 0, 0, 0))
    boxplot(formula, outline=F, yaxp=c(0, 1, 10), cex.axis=0.5,
		      main=Title(title), ylab=label, col=color)
    abline(a=0.5, b=0)
    dev.off() 
  }
  
  with(data, {
    Plot(reciprocity ~ X1.1, 'reciprocity-by-initiator', 'Reciprocity by initiator')
    Plot(reciprocity ~ X1.2, 'reciprocity-by-recipricator', 'Reciprocity by reciprocator')
    Plot(reciprocity ~ dyad, 'reciprocity-by-dyad', 'Reciprocity by dyad')
    Plot(reciprocity ~ unlist(dyads[dyad,1]), 'abs-reciprocity-by-chimp', 
	    'Reciprocity by individual')
  })
  
}

MakeWindowedDelayedScatterPlots <- function(windowed) {
  par(mfrow=c(2,2))
  window.sizes <- unique(windowed$ws)
  for(window.size in window.sizes) {
    with(windowed[windowed$treatment == 'delayed' & windowed$ws == window.size, ], {
    	windowed.lm <- lm(Y ~ X1)
    	windowed.lm.summary <- summary(windowed.lm)
	plot(X1, Y, xlim=c(0, 2000), ylim=c(0, 2000), 
	      main=sprintf("ws = %d, R^2 = %f", window.size, windowed.lm.summary$adj.r.squared))
	abline(windowed.lm)
	abline(0, 1)
    })
  }
}

WindowedDelayedRegressions <- function(windowed) {
  window.sizes <- unique(windowed$ws)
  return(unlist(Map(
    function(window.size) {
      with(windowed[windowed$treatment == 'delayed' & windowed$ws == window.size, ], {
	  windowed.lm <- lm(Y ~ X1)
	  windowed.lm.summary <- summary(windowed.lm)
	  return(windowed.lm.summary$adj.r.squared)
      })
    }, window.sizes)))
}

MakeWindowedScatterPlots <- function(windowed) {
  window.sizes <- unique(windowed$ws)
  treatments <- c('all', 'immediate', 'delayed')
  treatment.color <- c('red', 'blue', 'green')
  treatment.label <-
    c(expression(all), expression(Delta < 0), expression(Delta >= 0))
  for (i in 1:length(treatments)) {
    with(windowed[windowed$treatment == treatments[i],], {
      pdf(file=sprintf('%s/windowed-scatter-plot-%s.pdf', figsdir, treatments[i]), 
	    width=9, height=9)
      par(mfrow=c(1,1))
      window.size.factor <- factor(ws, labels=1:4)
      plot(log(X1 / 60), log(Y / 60), 
	    ylim = c(0, log(4000/60)),
	    xlim = c(0, log(4000/60)),
	    xlab = expression(log(sum(X))), ylab=expression(log(sum(Y))),
	    col = treatment.color[i],
# 	    main = treatment.label[i], 
	    main = '',
	    pch = as.integer(window.size.factor))
      legend('topleft', legend=paste0(window.sizes / 60, 'm'), pch=1:4, 
		  col=treatment.color[i])
      abline(0, 1)
      dev.off()
    })
  }
}

MakeNullModelScatterPlot <- function(dataset) {
  pdf(paste0(figsdir,'null-model-scatter-plot.pdf'), width=9, height=9)
  durations.by.dyad <-with(dataset, aggregate(dur.1 ~ dyad, FUN=sum))
  rsum <- 
    function(m, ws) c(ws, sum(runif(ws, max=m))/ws, sum(runif(ws, max=m))/ws)
  window.sizes <- c(20, 40, 60, 240)
  plot(xlim=c(5, 8), ylim=c(5, 8), 
	xlab=expression(log(sum(X))), 
	ylab=expression(log(sum(Y))), 
	x=c(), y=c())
  for(i in 1:length(window.sizes)) {
    null.model <- 
      as.data.frame(t(sapply(durations.by.dyad$dur.1, 
		      function(x) rsum(x, window.sizes[i]))))
    names(null.model) <- c('ws', 'X', 'Y')
    with(null.model, points(log(X), log(Y), pch=i)
    )
  }
  legend('topleft', legend=paste0(window.sizes, ''), pch=1:4)
  abline(0, 1)
  dev.off()
}

MakePDFs <- function(threshold = 0) {
  immediate <<- dataset[dataset$X2 < threshold, ]
  delayed <<- dataset[dataset$X2 >= threshold, ]
  MakeBoxPlots(dataset, color='gray')
  MakeBoxPlots(immediate, 'immediate', color='blue')
  MakeBoxPlots(delayed, 'delayed', color='green')
  
  immediate$treatment <<- 'immediate'
  delayed$treatment <<- 'delayed'
  
  pdf(paste0(figsdir,'reciprocity-by-treatment.pdf')) 
  with(
    rbind(immediate, delayed), 
	  boxplot(reciprocity ~ treatment, outline=F, 
	    col=c('green', 'blue'), yaxp=c(0, 1, 10), ylab=expression(rho),
	    names=c(
	      as.expression(substitute(Delta >= T, list(T=threshold))),
	      as.expression(substitute(Delta < T,  list(T=threshold)))
	    )
	  )
  )
  abline(a=0.5, b=0)
  dev.off()
  
  pdf(paste0(figsdir, 'delay-histogram.pdf'))
  with(dataset, qplot(X2 / (60), colour=I('red'), main='', xlab=expression(Delta)))
  dev.off()
  
  MakeNullModelScatterPlot(dataset)
  windowed <- CumulativeReciprocityByWindowSize(dataset)
  MakeWindowedScatterPlots(windowed)
}