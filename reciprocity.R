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
library(plyr)

datadir <- "../../data/"
figsdir <- "../../figs/"

# Raw data
d <-
  read.table(
    paste0(datadir, "raw_grooming_data.txt"),
    colClasses = "character",
    sep = ",",
    na.strings = "NA",
    blank.lines.skip = T,
    fill = T,
    col.names = 1:4
  )
d[, 2] <- sub(" ", "", d[, 2])
d[, 3] <- sub(" ", "", d[, 3])
d[, 4] <- sub(" ", "", d[, 4])

idx_schedule <- which(nchar(d[, 4]) == 0)
schedule <- d[idx_schedule, ]
nschedules <- nrow(schedule)

session <- matrix(0, nrow(d), 1)
session[idx_schedule, 1] <- 1
session <- cumsum(session)
d <- cbind(session, d)
d <- d[-idx_schedule, ]

chimp <- as.factor(sort(unique(c(d$X1, d$X2))))
nchmp <- length(chimp)
ndyad <- nchmp * (nchmp - 1) / 2
d$X1  <- factor(d$X1, levels = chimp)
d$X2  <- factor(d$X2, levels = chimp)
adjmat <- matrix(0, nchmp, nchmp)
rownames(adjmat) <- chimp
colnames(adjmat) <- chimp
re1 <- adjmat
re2 <- adjmat

adjmat2 <- matrix(NA, nchmp, nchmp)
rownames(adjmat2) <- chimp
colnames(adjmat2) <- chimp
adjmat2[upper.tri(adjmat2)] <- 0

comov <- read.csv(paste0(datadir, "comovement.csv"))
como2 <- cbind(comov[, 26], comov[, -26])
#como2 <- as.matrix(como2)
#como2[is.nan(como2)] <-0
como2[is.na(como2)] <- NaN
#como2 <- como2+t(como2)
rownames(como2) <- rownames(adjmat)
colnames(como2) <- colnames(adjmat)

proxim <- read.csv(paste0(datadir, "proxim.csv"))
#proxim <- as.matrix(proxim)
proxim <- proxim[, -1]
rownames(proxim) <- rownames(adjmat)
colnames(proxim) <- colnames(adjmat)

kin <- read.csv(paste0(datadir, "kin.csv"))
kin <- kin[, -1]
rownames(kin) <- rownames(kin)

#[1]       AL   BO   CL CR DL  FA  FR HI HO  HP JT KI  KK   KT   KY LA   LI MA NI   PA RO   SA SR   WH WI ZZ
age <-
  matrix(
    c(
      11.5,
      37.5,
      32.5,
      6,
      16,
      27.5,
      27,
      31,
      7.5,
      28,
      3,
      7.5,
      19.5,
      37.5,
      2 ,
      12.5,
      8 ,
      25,
      34.5,
      6,
      30.5,
      15,
      16.5,
      10,
      35,
      9.5
    ),
    26,
    1
  )

age2 <- age %*% t(rep(1, 26))
age2 <- (age2 + t(age2)) / 2
all(age == diag(age2))

for (j in 1:nrow(d)) {
  re2[as.numeric(d$X1[j]), as.numeric(d$X2[j])] <-
    re2[as.numeric(d$X1[j]), as.numeric(d$X2[j])] + 1
}
res <- re2 + t(re2)
re1[upper.tri(re1)] <- 1:ndyad
res[lower.tri(res)] <- NA
diag(res) <- NA

window_ws <- round(2 ^ seq(-1, 7, by = 0.5) * 60)
res_dyad  <- vector("list", ndyad)
dis_wind  <- vector("list", length(window_ws))
res_wind  <- matrix(NA, length(window_ws), 2)

group.colours <- c(within.bout= '#0000FF', delayed = '#00FF00')
percentile.colours <- c('#FF0000', '#001889')
window.size.colours <- c( "#7C0607", "#C18989", "#9394C0", "#1F28A2")


for (i in 1:nschedules) {
  now1 <- Sys.time()
  data <- d[d$session == i, ]
  
  t1 <-
    matrix(as.numeric(unlist(strsplit(data$X3, ":"))), nrow(data), 2, byrow =
             T)
  t1 <- t1[, 1] * 60 + t1[, 2]
  t2 <-
    matrix(as.numeric(unlist(strsplit(data$X4, ":"))), nrow(data), 2, byrow =
             T)
  t2 <- t2[, 1] * 60 + t2[, 2]
  dur <- t2 - t1
  idx <- which(dur < 0)
  dur[idx] <- dur[idx] + 60 * 60
  data <- cbind(data, dur)
  
  idx <- which(diff(t1) < 0)
  if (length(idx) > 0) {
    for (j in 1:length(idx)) {
      tmp <- (idx[j] + 1):min(idx[j + 1], nrow(data), na.rm = T)
      t1[tmp] <- t1[tmp] + j * 3600
    }
  }
  t3 <- t1 + dur
  data <- cbind(data, t1, t3)
  tmax <- max(t3)
  
  x <- cbind(rep(i, tmax), 1:tmax, rep(NA, tmax))
  
  if (i == 1)
  {
    d2 <- data
    timegrids <- x
  } else{
    d2 <- rbind(d2, data)
    timegrids <- rbind(timegrids, x)
  }
}

# Up till this point d2 contains single events with nornalised times

JoinEvents <- function(d2, pair.immediate=TRUE) {
  # Pair subsequent events from the same dyad and same session into a single row
  # so that the final result is a data-frame containing rows of dy-events.
  now1 <- Sys.time()
  result <- NULL
  for (k in 1:ndyad)
  {
    #  k=1
    idx_pair <- which(re1 == k, arr.ind = T)
    chimp1  <- chimp[idx_pair[1]]
    chimp2  <- chimp[idx_pair[2]]
    idx_pos <-
      which((d$X1 == chimp1 &
               d$X2 == chimp2) | (d$X1 == chimp2 & d$X2 == chimp1))
    idx_ses <- unique(d$session[idx_pos])
    data    <- d2[idx_pos, ]
    
    if (nrow(data) < 2)
      next
    
    ############################
    idx_session <- (table(data$session) >= 2)
    idx_session <- as.integer(names(idx_session)[idx_session])
    
    if (length(idx_session) == 0)
      next
    
    for (i in 1:length(idx_session)) {
      session <- data[data$session == idx_session[i], ]
      matched.index <- 1
      for (j in 2:nrow(session)) {
        # This event was reciprocation of previous event?
        if (pair.immediate) match <- j-1 else match <- matched.index
        if (session[j, 2] != session[match, 2]) {
          # Join reciprocation events
          result <- rbind(result, cbind(k, session[match, ], session[j, ]))
          matched.index <- j+1
        }
      }
    }
  }
  colnames(result) <-
    paste0(colnames(result), ".", rep(c(1, 2), each = 9))[-18]
  return(result)
}

tmp <- JoinEvents(d2)
# sum(tmp$session!=tmp$session.1)
# sum(tmp$X2.1!=tmp$X1.2)

Y <- tmp$t3.2 - tmp$t1.2
X1 <- tmp$t3.1 - tmp$t1.1
X2 <- tmp$t1.2 - tmp$t3.1
X3 <- tmp$t1.2 - tmp$t1.1
X4 <- tmp$t3.2 - tmp$t3.1
X5 <- tmp$t3.2 - tmp$t1.1
age_X1  <- age[tmp$X1.1]
age_Y   <- age[tmp$X2.1]
Z <- data.frame(Y, X1, X2, X3, X4) / 60
ind <- sign(X2)
kin_idx <- rep(NA, nrow(Z))
com_idx <- kin_idx
proxidx <- kin_idx #  proxim[tmp$X1.1,tmp$X2.1]
for (i in 1:nrow(Z)) {
  idx1 <- as.numeric(tmp$X1.1[i])
  idx2 <- as.numeric(tmp$X2.1[i])
  kin_idx[i] <- kin[idx1, idx2]
  com_idx[i] <- como2[idx1, idx2]
  proxidx[i] <- proxim[idx1, idx2]
}
Z <- cbind(Z, ind, kin_idx, com_idx, proxidx, age_X1, age_Y)

dataset <- tmp

# Duration of second event
dataset$Y <- tmp$t3.2 - tmp$t1.2

# Duration of first event
dataset$X1 <- tmp$t3.1 - tmp$t1.1

# Time between end of first event and start of second
dataset$X2 <- tmp$t1.2 - tmp$t3.1

# Time between start of first event and start of second
#  which is always positive.
dataset$X3 <- tmp$t1.2 - tmp$t1.1

# Time between end of first event and end of second event
dataset$X4 <- tmp$t3.2 - tmp$t3.1

# Time between start of first event and end of second event
dataset$X5 <- tmp$t3.2 - tmp$t1.1

Dyad <- function(i, j) {
  as.list(sort(c(i, j)))
}

Comovement <- function(i) {
  return(como2[dyads[[i, 1]], dyads[[i, 2]]])
}

DyadNumber <- function(i, j) {
  d <- Dyad(i, j)
  for (i in 1:nrow(dyads)) {
    if (unlist(dyads[i, 1]) == d[[1]] &&
        unlist(dyads[i, 2]) == d[[2]]) {
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

dym <- as.matrix(dataset[, c('X1.1', 'X1.2')])
dataset$dyad <-
  apply(dym, 1, function(x) {
    DyadNumber(x[[1]], x[[2]])
  })
dataset$comovement <- sapply(dataset$dyad, Comovement)

Reciprocity <- function(dataset) {
  with(dataset, {
    abs(Y - X1) / (X1 + Y)
  })
}

dataset$reciprocity = Reciprocity(dataset)


CumulativeReciprocity <- function(dataset, ws, treatment = 'all') {
  dataset$window <- floor(dataset$t1.1 / ws)
  SY <-
    with(dataset,  aggregate(Y ~ session.1 + dyad + window,  FUN = sum))
  SX1 <-
    with(dataset,  aggregate(X1 ~ session.1 + dyad + window,  FUN = sum))
  result <- merge(SY, SX1)
  result$R <- log(result$Y / result$X)
  #  names(result) = c('session', 'dyad', 'window', 'R')
  result$treatment = replicate(nrow(result), treatment)
  result$ws <- replicate(nrow(result), ws)
  return(result)
}

CumulativeReciprocityByDelay <- function(dataset, ws) {
  reciprocity.all <- CumulativeReciprocity(dataset, ws, 'all')
  reciprocity.delayed <-
    CumulativeReciprocity(dataset[dataset$X2 >= 0, ], ws, 'delayed')
  reciprocity.immediate <-
    CumulativeReciprocity(dataset[dataset$X3 == 0, ], ws, 'immediate')
  reciprocity.intra <-
    CumulativeReciprocity(dataset[dataset$X2 < 0, ], ws, 'intra')
  return(Reduce(
    rbind,
    list(
      reciprocity.all,
      reciprocity.delayed,
      reciprocity.immediate,
      reciprocity.intra
    )
  ))
}


ReciprocityByChimp <- function(data) {
  
  A <- data[, c('X1.1', 'reciprocity')]
  B <- data[, c('X1.2', 'reciprocity')]
  names(A) <- c('chimpanzee', 'reciprocity')
  names(B) <- names(A)
 
  return(rbind(A, B))
}


MakeBoxPlots <- function(data,
                         suffix = '',
                         color = 'blue') {
  FName <- function(name) {
    if (suffix != '') {
      fname.suffix = paste0('-', suffix)
    } else {
      fname.suffix = ''
    }
    return(paste0(figsdir, name, fname.suffix, '.pdf'))
  }
  
  Title <- function(main) {
    #     if (suffix != '') {
    #       return(paste0(main, ' (', suffix, ')'))
    #     } else {
    #       return(main)
    #     }
    return('')
  }
  
  Plot <- function(formula, fname, title, label = expression('Reciprocity (' ~ rho ~ ')')) {
    pdf(FName(fname))
    #par(mar=c(0, 0, 0, 0))
    boxplot(
      formula,
      outline = F,
      yaxp = c(0, 1, 10),
      cex.axis = 0.5,
      main = Title(title),
      ylab = label,
      col = color
    )
    abline(a = 0., b = 0.)
    dev.off()
  }
  
  with(data, {
    Plot(reciprocity ~ X1.1,
         'reciprocity-by-initiator',
         'Reciprocity by initiator')
    Plot(reciprocity ~ X1.2,
         'reciprocity-by-recipricator',
         'Reciprocity by reciprocator')
    Plot(reciprocity ~ dyad,
         'reciprocity-by-dyad',
         'Reciprocity by dyad')
         })
  
  A <- data[, c('X1.1', 'reciprocity')]
  B <- data[, c('X1.2', 'reciprocity')]
  names(A) <- c('chimpanzee', 'reciprocity')
  names(B) <- names(A)
  
  with(ReciprocityByChimp(data))
    Plot(
      reciprocity ~ chimpanzee,
      'abs-reciprocity-by-chimp',
      'Reciprocity by individual'
    )
  
}


WindowData <- function(dataset, window.sizes=c(20, 40, 60, 240) * 60) {
 
  windowed <-
    Reduce(rbind,
           Map(
             function(ws)
               CumulativeReciprocityByDelay(dataset, ws),
             window.sizes
           ))
           
  windowed$logX1 = log(windowed$X1 / 60)
  windowed$logY = log(windowed$Y / 60)
  windowed$reciprocity = Reciprocity(windowed)
  
  return(windowed)
}


MakeWindowedScatterPlots <- function(dataset, fname='windowed-scatter-plot') {

  window.sizes <- c(20, 40, 60, 240) * 60
  treatments <- c('all', 'immediate', 'delayed')
  windowed <- WindowData(dataset, window.sizes)

  treatment.label <-
    c(expression(all), expression(Delta < 0), expression(Delta >= 0))
  
  adply(1:length(treatments), 1, function(i) {
    with(windowed[windowed$treatment == treatments[i], ], {
        
      if (!is.null(fname)) {
        pdf(
          file = 
            sprintf('%s%s-%s.pdf', 
                      figsdir, fname, treatments[i]),
          width = 8,
          height = 8
        )
        
        par(mfrow = c(1, 1))
      }
          
      window.size.factor <- factor(window.sizes, labels = 1:4)
          
      if (!is.null(fname)) {
        plot(
          logX1,
          logY,
          ylim = c(0, log(4000 / 60)),
          xlim = c(0, log(4000 / 60)),
          xlab = expression(sum(log(X))),
          ylab = expression(sum(log(Y))),
          col = window.size.colours,
          #main = treatment.label[i],
          pch = as.integer(window.size.factor)
        )
      }
          
      regression.results <- adply(1:length(window.sizes), 1, function(ws.index) {
        lm.windowed <- lm(logY ~ logX1, 
                              data=windowed[windowed$ws == window.sizes[ws.index] &
                                windowed$treatment == treatments[i],])
              
        if (!is.null(fname)) {
          abline(lm.windowed$coefficients, 
                    lty=ws.index+1, col=window.size.colours[ws.index]) 
        }
        lm.summary <- summary(lm.windowed)
        print(lm.summary)
        data.frame(coef.1 =
         lm.windowed$coefficients[1], coef.2 = lm.windowed$coefficients[2], adj.r.squared = lm.summary$adj.r.squared, treatment=treatments[i], ws=window.sizes[ws.index])
          c(lm.windowed$coefficients[1], lm.windowed$coefficients[2], adj.r.squared = lm.summary$adj.r.squared)
      })
          
      if (!is.null(fname)) {
        r.squareds <- regression.results$adj.r.squared
        legend(
            'topleft',
            legend = c(paste0(window.sizes / 60, 'm (r-squared = ',
              round(r.squareds, digits=2), ')'), 'Y=X'),
            pch = c(1:4, NA),
            lty = c(2:5, 1),
            col = c(window.size.colours, 1)
        )
          
        abline(0, 1., lty=1)
        
        dev.off()
      }
          
      regression.results
    })
  })
}


MakeGroomingByDyadPlot <- function(dataset) {
  tmp3 <- table(as.factor(dataset$k.1))
  tmp3 <- tmp3/sum(tmp3)

  whi_dyads <- names(which(tmp3>quantile(tmp3,0.95)))
  sel_dyads <- as.factor(dataset$k.1)%in%whi_dyads
  high <- as.numeric(whi_dyads)
  
  dyad.percentage <- data.frame(sort(tmp3, decreasing=T))
  names(dyad.percentage) <- c('Dyad', 'Proportion') 
  dyad.percentage$percentile <- 
    unlist(lapply(dyad.percentage$Dyad, 
            function(x) if (any(high == x)) '1-5' else '6-100'))
  ggplot(dyad.percentage, aes(x=as.numeric(Dyad), y=Proportion, fill=percentile)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values=percentile.colours) + 
    xlab('Dyad') + ylab('Proportion of total grooming')
  ggsave(paste0(figsdir, 'grooming-distribution.pdf'))
  return(dyad.percentage)
}


MakeScatterPlots <- function(dataset) {
  pdf(paste0(figsdir, 'chimptest5a.pdf'))
  tmp3 <- table(as.factor(dataset$k.1))
  tmp3 <- tmp3/sum(tmp3)

  whi_dyads <- names(which(tmp3>quantile(tmp3,0.95)))
  sel_dyads <- as.factor(dataset$k.1)%in%whi_dyads

  out0 <- lm(Y~X1,data=dataset)
  summary(out0)
  out0a <- lm(Y~X1,data=dataset[sel_dyads,])
  summary(out0a)
  out0b <- lm(Y~X1,data=dataset[!sel_dyads,])
  summary(out0b)
  
  MakeLegend <- function() {
    legend("topleft", c("All","Top 5%","Bottom 95%", "Y=X"),
            col=c(1, percentile.colours[1], percentile.colours[2], 1), 
            pch=c(NA, 1, 4, NA), 
            lty=c(1, 1, 1, 3),
            lwd=1)
  }
 
  par(fig=c(0,1,0,1),new=F)
  par(mar=c(4.25,4.25,1,1))
  
  plot(dataset$X1, dataset$Y, type="n",
        xlab=expression(X), ylab="Y",
        asp=1,
        ylim=c(0, ceiling(max(c(dataset$X1, dataset$Y)))))
  
  abline(out0, col=1)
  abline(out0a, col=percentile.colours[1])
  
  points(dataset$X1[sel_dyads], dataset$Y[sel_dyads], col=percentile.colours[1])
  
  abline(out0b, col=percentile.colours[2])
  
  points(dataset$X1[!sel_dyads], dataset$Y[!sel_dyads], col=percentile.colours[2], pch=4)
  
  abline(0, 1, col=1, lty=3)
  
  MakeLegend()
 
  #par(oma=c(0,1,1,0), mar=c(4,0,0,0), fig=c(0.15,0.55,0.6,1),new=T)
#   barplot(sort(tmp3*100,decreasing=T),xaxt="n", yaxs="i", ylab="% of total grooming duration",width=1,xlab="Dyad",
#           col=c(rep(3,length(whi_dyads)),rep(4,length(tmp3)-length(whi_dyads))),
#           space=0,border=NA)
#   legend("topright",c("Top 5%","Bottom 95%"),col=c(3,4),,lwd=10)

  par(fig=c(0,1,0,1),new=F)
  ### Selected Dyads
  idx1 <- dataset$X2>=0
  idx2 <- dataset$X2<0

  tmp2 <- tmp[idx1,]
  tmp3 <- table(as.factor(tmp2$k.1))
  tmp3 <- tmp3/sum(tmp3)

  whi_dyads <- names(which(tmp3>quantile(tmp3,0.95)))
  sel_dyads <- as.factor(tmp2$k.1)%in%whi_dyads

  Z1 <- dataset[idx1,]
  out1 <- lm(Y~X1,data=Z1)
  summary(out1)
  out1a <- lm(Y~X1,data=Z1[sel_dyads,])
  summary(out1a)
  out1b <- lm(Y~X1,data=Z1[!sel_dyads,])
  summary(out1b)

  maxdur <- ceiling(max(c(Z1$X1,Z1$Y)))
  par(fig=c(0,1,0,1),new=F)
  par(mar=c(4.25,4.25,1,1))
  plot(Z1$X1, Z1$Y, type="n", 
        xlab=expression(X), ylab="Y", asp=1,
        ylim=c(0, maxdur), xlim=c(0, maxdur) ) 
  abline(out1,col=1)
  abline(out1a,col=percentile.colours[1])
  points(Z1$X1[sel_dyads], Z1$Y[sel_dyads], col=percentile.colours[1])
  abline(out1b,col=percentile.colours[2])
  points(Z1$X1[!sel_dyads], Z1$Y[!sel_dyads], col=percentile.colours[2], pch=4)
  abline(0, 1, col=1, lty=3)
  
  MakeLegend()
  #legend("topright",c("All","Top 5%","Bottom 95%"),col=c(2,3,4),pch=c(NA,1,4),lwd=1)

#   par(fig=c(0.15,0.65,0.6,1),new=T)
#   barplot(sort(tmp3,decreasing=T),xaxt="n", ylab="% of Grooming",width=1,#xlab="Dyad ID"
#           col=c(rep(3,length(whi_dyads)),rep(4,length(tmp3)-length(whi_dyads))),
#           space=0,border=NA)

  ### 

  tmp2 <- tmp[idx2,]
  tmp3 <- table(as.factor(tmp2$k.1))
  tmp3 <- tmp3/sum(tmp3)

  whi_dyads <- names(which(tmp3>quantile(tmp3,0.95)))
  sel_dyads <- as.factor(tmp2$k.1)%in%whi_dyads

  Z2 <- Z[idx2,]
  out2 <- lm(Y~X1,data=Z2)
  summary(out2)
  out2a <- lm(Y~X1,data=Z2[sel_dyads,])
  summary(out2a)
  out2b <- lm(Y~X1,data=Z2[!sel_dyads,])
  summary(out2b)

  par(fig=c(0,1,0,1),new=F)
  par(mar=c(4.25,4.25,1,1))
  
  plot(Z2$X1,Z2$Y, type="n",
        xlab=expression(X), ylab="Y", asp=1,
        ylim=c(0, ceiling(max(c(Z2$X1,Z2$Y)))))
  
  abline(out2,col=1)
  abline(out2a,col=percentile.colours[1])
  
  points(Z2$X1[sel_dyads], Z2$Y[sel_dyads], col=percentile.colours[1])
  
  abline(out2b,col=percentile.colours[2])
  
  points(Z2$X1[!sel_dyads], Z2$Y[!sel_dyads], col=percentile.colours[2], pch=4)
  abline(0, 1, col=1, lty=3)
  
  MakeLegend()
#   legend("topright",c("All","Top 5%","Bottom 95%"),col=c(2,3,4),pch=c(NA,1,4),lwd=1)

  par(fig=c(0.15,0.65,0.6,1), new=T)
#   barplot(sort(tmp3,decreasing=T),xaxt="n", ylab="% of Grooming",width=1,#xlab="Dyad ID"
#           col=c(rep(3,length(whi_dyads)),rep(4,length(tmp3)-length(whi_dyads))),
#           space=0,border=NA)

  par(fig=c(0,1,0,1), new=F)
  plot(dataset$X1, dataset$Y, xlab=expression(X), ylab="Y",type="n", asp=1,
      ylim=c(0, ceiling(max(c(dataset$X1, dataset$Y)))))
  abline(out0, col=1)
  abline(out1, col=3)
  abline(out2, col=4)
  abline(0, 1, col=1, lty=3)
  points(dataset$X1[idx1], dataset$Y[idx1], col=3)
  points(dataset$X1[idx2], dataset$Y[idx2], col=4, pch=4)
  legend("topleft",c("All",expression(Delta >= 0), expression(Delta < 0), "Y=X"),
            col=c(1, 3, 4, 1), lwd=1, pch=c(NA, 1, 4, NA), lty=c(1, 1, 1, 3))
  # par(fig=c(0,0.5,0.5,1),new=T)
  # pie(c(sum(idx1),sum(idx2)),col=3:4,labels="",radius=0.5)

  dev.off()

  pdf(paste0(figsdir,"chimptest5f.pdf"))

  Z2 <- data.frame(Y,X1,X2)/60
  tgrid <- seq(-30,45,1)
  delta_t <- seq(20,30,5)
  cors <- matrix(NA,length(tgrid),length(delta_t))#rep(NA,length(tgrid))
  regs <- cors
  cors2 <- cors
  regs2 <- cors
  cor.uci<- cors
  cor.lci<- cors
  reg.uci<- cors
  reg.lci<- cors
  cor.uci2<- cors
  cor.lci2<- cors
  reg.uci2<- cors
  reg.lci2<- cors
  rsquared1 <- cors
  rsquared2 <- cors

  for(j in 1:length(delta_t)){
    for(i in 1:length(tgrid)){
      idx <- which(Z2[,3]>(tgrid[i]-delta_t[j]) & Z2[,3]<=(tgrid[i]+delta_t[j]))
  #    if(tgrid[i]<=0){
  #      idx <- which(Z2[,3]>(tgrid[i]-delta_t[j]) & Z2[,3]<=(tgrid[i]))
  #    }else{
  #      idx <- which(Z2[,3]>(tgrid[i]) & Z2[,3]<=(tgrid[i]+delta_t[j]))
  #   }
      cors[i,j] <- cor(Z2[idx,1:2])[1,2]
      cor.lci[i,j] <- cor.test(Z2[idx,1],Z2[idx,2])$conf.int[1]
      cor.uci[i,j] <- cor.test(Z2[idx,1],Z2[idx,2])$conf.int[2]
      fit <- (lm(Z2[idx,1]~Z2[idx,2]))
      regs[i,j] <- summary(fit)$coefficients[2,1]
      rsquared1[i,j] <- summary(fit)$r.squared
      reg.lci[i,j] <- confint(fit)[2,1]
      reg.uci[i,j] <- confint(fit)[2,2]
    }
  }

  tmp2 <- tmp
  tmp3 <- table(as.factor(tmp2$k.1))
  tmp3 <- tmp3/sum(tmp3)

  whi_dyads <- names(which(tmp3>quantile(tmp3,0.95)))
  sel_dyads <- as.factor(tmp2$k.1)%in%whi_dyads

  Z22 <- data.frame(Y[!sel_dyads],X1[!sel_dyads],X2[!sel_dyads])/60

  for(j in 1:length(delta_t)){
    for(i in 1:length(tgrid)){
      idx <- which(Z22[,3]>(tgrid[i]-delta_t[j]) & Z22[,3]<=(tgrid[i]+delta_t[j]))
      cors2[i,j] <- cor(Z22[idx,1:2])[1,2]
      cor.lci2[i,j] <- cor.test(Z22[idx,1],Z22[idx,2])$conf.int[1]
      cor.uci2[i,j] <- cor.test(Z22[idx,1],Z22[idx,2])$conf.int[2]
      fit <- (lm(Z22[idx,1]~Z22[idx,2]))
      regs2[i,j] <- summary(fit)$coefficients[2,1]
      rsquared2[i,j] <- summary(fit)$r.squared
      reg.lci2[i,j] <- confint(fit)[2,1]
      reg.uci2[i,j] <- confint(fit)[2,2]
    }
  }

  #x11()
  #par(mfrow=c(2,1))
  #par(mar=c(4.25,4.25,1,1))

  #hist(X2/60,50,xlim=c(-30,30),xlab=expression(X[2]~"[mins.]"),main="")
  #matplot(tgrid,regs,xlab="X2 [mins.]",ylab=expression(beta),type="l",xlim=c(-30,30))#,ylim=c(-1,1))
  #abline(v=0,col=8)
  #matplot(tgrid,rsquared1,xlab="X2 [mins.]",ylab=expression(R^{2}),type="l",xlim=c(-30,30),ylim=c(0,0.6))
  #abline(v=0,col=8)
  #legend("topright",paste(delta_t," min,"),title="Window Width",lwd=1,col=1:4,lty=1:4)

  #plot(tgrid,cors[,1],xlab=expression(X[2]~"[mins.]"),
  #     ylab=expression("Corr((Y,"~X[1]~")|"~X[2]~")"),type="l",ylim=c(0,0.8))
  #lines(tgrid,cor.uci[,1],col=2,lty=2)
  #lines(tgrid,cor.lci[,1],col=2,lty=2)

  #library("Hmisc")
  #errbar(tgrid,cors[,1],cor.uci[,1],cor.lci[,1], add=T, pch=1, cap=.1)

  par(mfrow=c(1,1))
  par(mar=c(4.25,4.25,1,1))

  hist(X2/60,50,xlim=c(-30,45),xlab=expression('Delay (' ~ Delta ~ ') [mins.]'),main="")

  i<-1
  plot(tgrid,regs[,i],xlab=expression(Delta~"[mins.]"),col=2,lwd=2,
      ylab=expression(beta),type="l",ylim=c(0,2))
  lines(tgrid,reg.uci[,i],col=2,lty=2)
  lines(tgrid,reg.lci[,i],col=2,lty=2)


  lines(tgrid,regs2[,i],col=4,lwd=2)
  lines(tgrid,reg.uci2[,i],col=4,lty=2)
  lines(tgrid,reg.lci2[,i],col=4,lty=2)
  legend("topleft",c(expression(beta~"(All)"),"95% C.I. (All)",
                      expression(beta~"(Bottom 95%)"),"95% C.I. (Bottom 95%)"),
        col=c(2,2,4,4),lwd=c(2,1,2,1),lty=c(1,2,1,2))

  #x11()
  plot(tgrid,rsquared1[,1],ylim=c(0,1),col=2,type="l",ylab=expression(R^{2}),xlab=expression(Delta~"[mins.]"))
  lines(tgrid,rsquared2[,1],col=4,lty=1)
  legend("topright",c(expression(R^{2}~"(All)"),
                      expression(R^{2}~"(Bottom 95%)")),
        col=c(2,4),lwd=1,lty=1)

  #x11()
  #plot(Z2[idx,2],Z2[idx,1])

  dev.off()

}


MakeDeltaVersusRho <- function() {
  matplot(abs(dataset$X2), dataset$reciprocity, pch=4, xlab=expression('Absolute delay (' ~ abs(Delta) ~ ')'), ylab=expression('Reciprocity (' ~ rho ~ ')'))
  with(dataset[dataset$reciprocity == 0,], points(abs(X2), reciprocity, col='red', pch=4))
}


MakeCountBarPlot <- function(dataset) {

  DF <- function(dataset, condition) {
    counts <- aggregate(reciprocity ~ chimpanzee, ReciprocityByChimp(dataset), FUN=length)
    n <- nrow(counts)
    data.frame(x=counts$chimpanzee, y=counts$reciprocity, condition=rep(condition, n))
  }
  
  df <- rbind(DF(dataset[dataset$X2 < 0,], 'within.bout'), DF(dataset[dataset$X2 >= 0,], 'delayed'))
  
  ggplot(df, aes(x=x, y=y, fill=condition)) + 
    geom_bar(stat='identity', position=position_dodge()) +
    scale_fill_manual(values=group.colours) +
    xlab('Chimpanzee') + ylab('n')
  ggsave(paste0(figsdir, 'counts-bar-plot.pdf'))
}


NullModel <- function(dataset) {
# 
#   mean.duration.by.chimp <- with(dataset, aggregate(dur.1 ~ X1.1, FUN=mean))
#   
#   MeanDurationForChimp <- function(i) 
#     mean.duration.by.chimp[mean.duration.by.chimp$X1.1 == i,]$dur.1
#   
#   mean.dur.1 <- sapply(dataset$X1.1, MeanDurationForChimp)
#   mean.dur.2 <- sapply(dataset$X1.2, MeanDurationForChimp)
  
  #RandomDuration <- function(mean_duration) runif(1, min=0, max=2*mean_duration)
  RandomDuration <- function(chimp) {
    durations <- dataset[dataset$X1.1 == chimp, 'dur.1']
    if (length(durations) > 0) {
      return(sample(durations, size=1, replace=T))
    } else {
      return(NA)
    }
  }
  
  null.model <- dataset
  null.model$X1 <- sapply(dataset$X1.1, RandomDuration)
  null.model$Y <-  sapply(dataset$X1.2, RandomDuration)
  
  return(null.model)
}


BootstrapNullModel <- function(dataset, n=100) {
  results <- c()
  for(i in 1:n) {
    results <- 
      rbind(results, 
              aggregate(reciprocity ~ ws + treatment,
                WindowData(NullModel(dataset)), FUN=median))
  }
  return(results)
}


NullModelCriticalValues <- function(bootstrap.results, var, p=0.05) {

  CriticalValues <- function(probability) {
     df <- aggregate(as.formula(paste(var,'~ ws + treatment')), 
                      bootstrap.results, 
                        FUN=function(x) quantile(x, p=probability))
     df
  }
  
  lower <- CriticalValues(p/2.)
  upper <- CriticalValues(1. - (p/2.))
  
  result <- cbind(lower[,1:2], lower[,3], upper[,3])
  names(result)[3] <- paste0('lower.', var)
  names(result)[4] <- paste0('upper.', var)
  
  return(result)
}


MakePDFs <- function(threshold = 0) {
  immediate <<- dataset[dataset$X2 < threshold,]
  delayed <<- dataset[dataset$X2 >= threshold,]
  MakeBoxPlots(dataset, color = 'gray')
  MakeBoxPlots(immediate, 'immediate', color = 'blue')
  MakeBoxPlots(delayed, 'delayed', color = 'green')
  
  immediate$treatment <<- 'immediate'
  delayed$treatment <<- 'delayed'
  
  pdf(paste0(figsdir, 'reciprocity-by-treatment.pdf'))
  with(
    rbind(immediate, delayed),
    boxplot(
      reciprocity ~ treatment,
      outline = F,
      col = c('green', 'blue'),
      yaxp = c(0, 1, 10),
      ylab = expression('Reciprocity (' ~ rho ~ ')'),
      names = c(as.expression(substitute(
        Delta >= T, list(T = threshold)
      )),
      as.expression(substitute(
        Delta < T,  list(T = threshold)
      )))
    )
  )
  abline(a = 0.0, b = 0.)
  dev.off()
  
  dataset$group <- rep(NaN, nrow(dataset))
  dataset[dataset$X2 < 0,]$group = 'within.bout'
  dataset[dataset$X2 >= 0,]$group = 'delayed'
  ggplot(dataset, aes(x=X2/60, group=group, fill=group)) + 
    geom_histogram(closed='left', breaks=seq(-44, 266, by=44/10)) +
    scale_fill_manual(values=group.colours) +
    xlab(expression('Delay (' ~ Delta ~ ') in minutes'))
  ggsave(paste0(figsdir, 'delay-histogram.pdf'))
  
  MakeCountBarPlot(dataset)
  MakeGroomingByDyadPlot(dataset)
  MakeScatterPlots(dataset)
  MakeWindowedScatterPlots(dataset)
}
