require(foreach)
require(doMC)
require(aroma.light)
require(nlme)
require(mgcv)
require(mnormt)

### Wrapper function for estimating the bandwidth for the kernel density estimate at each observation time

bw.est<-function(dat,bw="SJ",record) {
   try(registerDoMC(cores = 6))
   BW<- foreach(i = 1:length(dat), .combine = 'c') %dopar% {
       t.bw<- density(dat[[i]],bw=bw)$bw
       print(paste("BW estimate for t=",record[i],"is",signif(t.bw,3)))
       t.bw
      }
   BW
}

### Estimate the density or its derivative with respect to x
den.fun <- function(dat, x, deriv = 0, trim = TRUE, low.thresh = NULL, hi.thresh = NULL,bw="nrd0") {
# `dat' is a vector of observations sampled from the population at a single observation time.
# `x' is a vector of locations at which the density (or its derivative) is to be estimated
# `deriv' is an indicator:  0 corresponds to the density, 1 corresponds to the derivative with respect to x
# `trim' is an indicator: if TRUE, then estimates from the distribution tails are made to be NA; if FALSE, the estimates for all values of `x' are returned
# `low.thresh' is a scalar specifying the lower bound of the linear dynamic range (used for censored data)
# `hi.thresh' is a scalar specifying the upper bound of the linear dynamic range (used for censored data)
# `bw' bandwidth argument passed to density function

  dat2 <- dat
  prop.thresh <- 0
  if (!is.null(low.thresh)) 
    if (is.na(low.thresh)) 
      low.thresh <- NULL
  if (!is.null(hi.thresh)) 
    if (is.na(hi.thresh)) 
      hi.thresh <- NULL
  if (is.numeric(low.thresh)) {
    # dat2 <- c(dat[dat > low.thresh], 2 * low.thresh - dat[dat > low.thresh])
    prop.thresh <- prop.thresh + mean(dat <= low.thresh)
  }
  if (is.numeric(hi.thresh)) {
#     dat2 <- c(dat2[dat2 < hi.thresh], 2 * hi.thresh - dat2[dat2 < hi.thresh])
    prop.thresh <- prop.thresh + mean(dat >= hi.thresh)
  }
  est <- density(dat2,  n = 5000, bw=bw)
  if (is.numeric(low.thresh)) {
    est$y <- est$y[est$x >= low.thresh]
    est$x <- est$x[est$x >= low.thresh]
  }
  if (is.numeric(hi.thresh)) {
    est$y <- est$y[est$x <= hi.thresh]
    est$x <- est$x[est$x <= hi.thresh]
  }
  est$y <- (1 - prop.thresh) * est$y/(sum(est$y) * mean(diff(est$x)))
  t.fun <- splinefun(x = est$x, y = est$y)
  z <- t.fun(x, deriv)
  if (trim) 
    z[x < min(est$x) | x > max(est$x)] <- NA
  if (deriv == 0) 
    z[z < 0] <- 0
  z
}

### Transform from logistic scale to linear scale
unlogit <- function(x) exp(x)/(1 + exp(x))

### Transform from linear scale to logistic scale
logit <- function(x) log(x/(1 - x))

### Goodness of fit test used when fitting splines to cumulative distribution functions F(x,t)
GOF.test <- function(X, N, fit, logit.p = FALSE, fit.P = NULL) {
  if (is.vector(fit)) 
    p.fit <- fit else p.fit <- predict(fit, fit$x)$y
  if (logit.p) 
    p.fit <- unlogit(p.fit)
  if (!is.null(fit.P)) 
    p.fit <- fit.P + predict(fit, fit$x)$y
  p.fit[p.fit > 1] <- 1
  p.fit[p.fit < 0] <- 0
  alpha <- X + 1
  beta <- N - X + 1
  P.x <- alpha/(alpha + beta)
  t.stat <- 2 * (sum(dbinom(X, N, P.x, log = TRUE)) - sum(dbinom(X, N, p.fit, log = TRUE)))
  if (is.vector(fit)) 
    t.pval <- 1 - pchisq(t.stat, length(X)) else t.pval <- 1 - pchisq(t.stat, length(X) - fit$df)
  if (is.na(t.pval)) 
    stop(paste("GOF.test pval is", t.pval, "for iteration", i))
  t.pval
}

### Make semi-transparent colors.  Code for this function was taken from Nick Sabbe's post on stackoverflow.com

makeTransparent <- function(someColor, alpha = 100) {
  newColor <- col2rgb(someColor)
  apply(newColor, 2, function(curcoldata) {
    rgb(red = curcoldata[1], green = curcoldata[2], blue = curcoldata[3], alpha = alpha, 
      maxColorValue = 255)
  })
}


## Define colors to be used for plots.
COL <- c(makeTransparent(c(2:6)), makeTransparent(c("orange")),1)


### Estimate the density, its derivative with respect to x, and the probability flux.
DnA.ingredients <- function(x.loc, record, dat, Plot = FALSE, hi.thresh = NULL, low.thresh = NULL,bw,n.iter=20){


# `x.loc' is a vector of locations at which the density (or its derivative) is to be estimated

# `record' is a vector specifying the observation times of the population

# `dat' is a list.  Each element of dat is vector of observations sampled from the population at a single observation time.

# `Plot' is an indicator:  if TRUE, the estimated densities will be plotted as will the spline fits to the cumulative distribution functions (F(x,t)); if FALSE, these plots are supressed

# `trim' is an indicator: if TRUE, then estimates from the distribution tails are made to be NA; if FALSE, the estimates for all values of `x' are returned

# `hi.thresh' is a vector specifying the upper bound of the linear dynamic range at each observation time(used for censored data) 

# `low.thresh' is a vector specifying the lower bound of the linear dynamic range at each observation time(used for censored data)

# `bw' bandwidth argument passed to density function

# n.iter is an integer specifying the maximum number of reweighting iterations to use when fitting splines to the cumalative distribution functions


  if (length(hi.thresh == 1)) 
    hi.thresh <- rep(hi.thresh, length(dat))
  if (length(low.thresh == 1)) 
    low.thresh <- rep(low.thresh, length(dat))

  ### Estimate time derivative of P(x,t)
  N <- sapply(dat, length)
  cum.Y<-sapply(dat,function(z){
  zz<-floor(approx(y=1:length(z),x=sort(z),xout=x.loc)$y)
               zz[x.loc<min(z)]<-0
               zz[x.loc>max(z)]<-length(z)
               zz
  })

  ### How many observations fall within each bin formed by adjacent x.loc values
  Y<-rbind(cum.Y[1,],cum.Y[-1,]-cum.Y[-nrow(cum.Y),])
  
  ## Replace counts for x.loc values outside of dynamic range with NA
for(i in 1:ncol(Y)){
 if(!is.null(low.thresh)) if(!is.na(low.thresh[i]))  if(any(x.loc<low.thresh[i])) Y[x.loc<low.thresh[i],i]<-NA
 if(!is.null(hi.thresh))  if(!is.na(hi.thresh[i])) if(any(x.loc>hi.thresh[i])) Y[x.loc>hi.thresh[i],i]<-NA
}

######################################
### This step needs work to handle ###
### left truncation:               ###
    cum.Y<- apply(Y,2,cumsum)
######################################  

  
dP_dt<-NULL

 

  w.mod<-rep(1,length(record))
  p.smooth<-smooth.it<-0
  
  
#while((any(p.smooth<.01)|(!any(p.smooth<.9)))&smooth.it<n.iter){
while(any(p.smooth<.01)&smooth.it<n.iter){ ### Allows overfitting
  plot.fit.P<-fit.P <- dP_dt <- NULL

  for (i in 1:length(x.loc)) {
    spline.dat<-data.frame(y=cum.Y[i,],n=N,x=record)
 ##If applicable, exclude time points for which x.loc[i] is outside of dynamic range (low.thresh,hi.thresh) 
    t.drop<-NULL
    t.w.mod<-w.mod
    if(!is.null(low.thresh)) t.drop<-c(t.drop,which(low.thresh>x.loc[i]))
    if(!is.null(hi.thresh)) t.drop<-c(t.drop,which(hi.thresh<x.loc[i]))
    t.drop<-sort(unique(t.drop))
    t.keep<-(1:length(record))[!(1:length(record))%in%t.drop]
    spline.dat<-spline.dat[t.keep,]
    t.w.mod<-w.mod[t.keep]
    b<- gam(cbind(y,n-y)~s(x,k=min(c(20,nrow(spline.dat)))),data=spline.dat,family="binomial",weights=t.w.mod)
    t.fit.P<-rep(NA,length(record))
    t.fit.P[t.keep]<-unlogit(as.vector(predict(b)))
    fit.P<-rbind(fit.P,t.fit.P)
    if (Plot) {
      t.xx <- seq(min(record), max(record), length.out = 200)
      plot.fit.P<-rbind(plot.fit.P,unlogit(predict(b,newdata=list(x=t.xx))))
    }

    ###Numerically estimate slope of P(x.loc[i],t) w.r.t t 
    P1<-predict(b,newdata=list(x=record[t.keep]+1e-5))
    P2<-predict(b,newdata=list(x=record[t.keep]-1e-5))
    t.dP_dt<-rep(NA,length(record))
    t.dP_dt[t.keep]<-(unlogit(P1)-unlogit(P2))/2e-5
    dP_dt<-rbind(dP_dt,t.dP_dt)
  }


      fit.P[fit.P > 1 &!is.na(fit.P)] <- 1
      fit.P[fit.P < 0 &!is.na(fit.P)] <- 0
      t.stat <- NULL
      t.fit<-fit.P-rbind(0,fit.P[-nrow(fit.P),])
      t.fit<-rbind(t.fit,1-colSums(t.fit,na.rm=TRUE))
      
      t.Y<-rbind(Y,N-colSums(Y,na.rm=TRUE))
      expected<-scale(t.fit,scale=1/N,center=FALSE)
      t.stat<-colSums( (t.Y-expected)^2/expected,na.rm=TRUE)
if(0){
      for (i in 1:ncol(Y)){
        tt.stat<-t.Y[,i]*log(t.fit[,i]*N[i]/t.Y[,i])
        tt.stat[Y[,i]==0]<-0
        t.stat <- c(t.stat,-2*sum(tt.stat,na.rm=TRUE)/(1+(sum(1/t.fit[,i],na.rm=TRUE)-1)/(6*sum(Y[,i],na.rm=TRUE)*(sum(!is.na(Y[,i]))-1))))
      }
    }
      p.smooth <- 1 - pchisq(t.stat, colSums(!is.na(Y)))
#      if(!any(p.smooth<.9)) k.deduct<-k.deduct+1   ### Used when trying to control overfitting
      smooth.it <- smooth.it + 1
      w.mod[which.max(t.stat)] <- w.mod[which.max(t.stat)] * exp(1- smooth.it/n.iter)
}

  ### Plot the emperical CDFs and their fitted splines
  if(Plot){
    dev.new(width = 12, height = 10)
    nf <- layout(matrix(1:2, 1, 2), widths = c(1, 0.4))
    par(mai = c(1, 1, 0.2, 0.1))
    plot(-1e+06, 1e+06, ylab = "Pr(x<X)", xlab = "Time", ylim = c(0, 1), xlim = range(record), 
      main = "", cex.axis = 1.7, cex.lab = 1.7)
    grid()
    t.xx <- seq(min(record), max(record), length.out = 200)
    for (i in 1:length(x.loc)) {
      points(record, cum.Y[i,]/N, col = COL[i%%length(COL)+1])
      lines(t.xx,plot.fit.P[i,], lwd = 2, col = COL[i%%length(COL)+1], lty = 1)
    }
    par(mai = c(1, 0.1, 0.2, 0.1))
    plot(-1e+06, 1e+06, ylab = "", xlab = "", ylim = 0:1, xlim = 0:1, xaxt = "n", 
      yaxt = "n", main = "")
    legend("topright", title = "X", legend = signif(x.loc, 3), ncol = 2, col = COL[c(2:length(COL),1)], 
    lwd = 3, cex = 1.3, bty = "n")
  }

  p <- dp_dx <- NULL
  for (ind in 1:length(dat)) {
    ### Estimate x derivative of log(p(x.loc,record[ind]))
    dp_dx <- cbind(dp_dx, den.fun(dat = dat[[ind]], x = x.loc, deriv = 1,bw=bw[ind],  
      low.thresh = low.thresh[ind], hi.thresh = hi.thresh[ind]))
    p <- cbind(p, den.fun(dat = dat[[ind]], x = x.loc, deriv = 0, bw=bw[ind], 
      low.thresh = low.thresh[ind], hi.thresh = hi.thresh[ind]))
    ### Prevent estimates from tails of distributions
    tails<-which(x.loc<quantile(dat[[ind]],.01)|x.loc>quantile(dat[[ind]],.99))
    dP_dt[tails,ind]<-p[tails,ind]<-dp_dx[tails,ind]<-NA
  }
  rownames(p)<-x.loc
  colnames(p)<-record
  list(dP_dt = dP_dt, dp_dx = dp_dx, p = p)
}

Is.numeric <- function(x) {
  z <- FALSE
  if (!is.null(x)) {
    if (any(!is.na(x))) {
      if (!(any(!is.numeric(x[!is.na(x)])))) 
        z <- TRUE
    }
  }
  z
}


### Estimate the relative entropy between each observed distribution the last observed distribution
rel.entropy <- function(dat,bw="nrd0", hi.thresh = NULL, low.thresh = NULL) {
  p.low.ref <- p.hi.ref <- p.low <- p.hi <- rel.ent <- rep(0, length(dat))
  n <- length(dat)
  
  if (length(hi.thresh == 1)) 
    hi.thresh <- rep(hi.thresh, n)
  if (length(low.thresh == 1)) 
    low.thresh <- rep(low.thresh, n)
  if (length(bw == 1)) 
    bw <- rep(bw, n)
  for (i in 1:n) {
    xx <- seq(min(dat[[i]]), max(dat[[i]]), length.out = 1000)
    t.low.thresh <- t.hi.thresh <- NULL
    if (Is.numeric(low.thresh[c(i, n)])) {
      t.low.thresh <- max(low.thresh[c(i, n)], na.rm = TRUE)
      xx <- xx[xx > t.low.thresh]
      p.low[i] <- mean(dat[[i]] <= t.low.thresh)
      p.low.ref[i] <- mean(dat[[length(dat)]] <= t.low.thresh)
    }
    if (Is.numeric(hi.thresh[c(i, n)])) {
      t.hi.thresh <- min(hi.thresh[c(i, n)], na.rm = TRUE)
      xx <- xx[xx < t.hi.thresh]
      p.hi[i] <- mean(dat[[i]] >= t.hi.thresh)
      p.hi.ref[i] <- mean(dat[[length(dat)]] >= t.hi.thresh)
    }
    t.den1 <- den.fun(dat = dat[[i]], bw=bw[i],x = xx, deriv = 0, low.thresh = t.low.thresh, 
      hi.thresh = t.hi.thresh)
    tot.prob <- sum(t.den1, na.rm = TRUE) * mean(diff(xx)) + p.low[i] + p.hi[i]
    if (abs(tot.prob - 1) > 0.005) 
      print(paste("Sum of estimated probabilities for data list element", i, 
        "is", round(tot.prob, 3)))
    t.den2 <- den.fun(dat = dat[[n]],bw=bw[i], x = xx, deriv = 0, low.thresh = t.low.thresh, 
      hi.thresh = t.hi.thresh)
    t.rel <- log(t.den1/t.den2)
    t.rel[t.den1 == 0] <- 0
    rel.ent[i] <- sum(t.den1 * t.rel * mean(diff(xx)), na.rm = TRUE)
  }
  suppressWarnings(rel.ent[p.low > 0] <- (rel.ent + p.low * log(p.low/p.low.ref))[p.low > 
    0])
  suppressWarnings(rel.ent[p.hi > 0] <- suppressWarnings(rel.ent + p.hi * log(p.hi/p.hi.ref))[p.hi > 
    0])
  rel.ent
}

DnA.analyze <- function(dat, x.loc, record, n.boot = 100, Adj = 1, hi.thresh = NULL, bw.method="nrd0",unnorm=FALSE,
  low.thresh = NULL, dat.mat = NULL,mat.fun=NULL,plot.den=FALSE) {
  n.t <- length(dat)
  if (is.null(n.boot) | n.boot < (2 * length(x.loc))) {
    print(paste("n.boot has been set to the minimum required number of 2*length(x.loc) = ", 
      2 * length(x.loc)))
    n.boot <- 2 * length(x.loc)
  }
  if (length(hi.thresh) == 1) 
    hi.thresh <- rep(hi.thresh, n.t)
  if (length(low.thresh) == 1) 
    low.thresh <- rep(low.thresh, n.t)
  if (!length(hi.thresh) %in% c(0, 1, n.t)) 
    stop(paste("length(hi.thresh)=", length(hi.thresh), ".  hi.thresh must be one of: NULL (no threshold), scalar (constant threshold), or a vector of length(dat)=", 
      n.t, " (providing threshold for each time).", sep = ""))
  if (!length(low.thresh) %in% c(0, 1, n.t)) 
    stop(paste("length(low.thresh)=", length(low.thresh), ".  low.thresh must be one of: NULL (no threshold), scalar (constant threshold), or a vector of length(dat)=", 
      n.t, " (providing threshold for each time).", sep = ""))
  if (length(record) != length(dat)) 
    stop(paste(length(record), " = length(record) != length(dat) = ", length(dat), 
      ". length(record) must equal length(dat).", sep = ""))
  if (!is.numeric(Adj) | Adj < 0) {
    warning(paste("Adj=", Adj, ".  Adj must be a positive real number and has been set to 1.", 
      sep = ""))
    Adj <- 1
  }
  if (!is.vector(x.loc) | !is.numeric(x.loc)) 
    stop("x.loc must be a numeric vector")
  x.loc2 <- seq(min(x.loc, na.rm = TRUE), max(x.loc, na.rm = TRUE), length.out = 200)
  max.den <- 0

### Estimate and plot (if plot.den==TRUE) densities
  for (i in c(which.min(record), which.max(record))) {
    t.den <- den.fun(dat = dat[[i]], x = x.loc2, deriv = 0, low.thresh = low.thresh[i], 
      hi.thresh = hi.thresh[i])
    max.den <- max(c(max.den, t.den[t.den < Inf]), na.rm = TRUE)
  }
  if (n.t > 49) {
    nc <- 7
    nr <- 7
  } else {
    nc <- ceiling(sqrt(n.t))
    nr <- ceiling(n.t/nc)
  }
  keep.den <- x.loc2
  bw<-bw.est(dat,bw=bw.method,record=record)*Adj
  for (i in 1:n.t) {
      keep.den <- cbind(keep.den, den.fun(dat = dat[[i]], bw=bw[i],x = x.loc2, deriv = 0, 
      low.thresh = low.thresh[i], hi.thresh = hi.thresh[i]))
    
    if(plot.den){
      if (i %in% (1 + nr * nc * (0:10))) {
        dev.new(width = 12, height = 12)
        nf <- layout(matrix(1:(nc * nr), nr, nc, byrow = TRUE), widths = c(1.1 + 
          0.05 * nc, rep(1, nc - 1)), heights = c(rep(1, nr - 1), 1.1 + 0.05 * 
          nr))
      }
      b <- l <- t <- r <- 0
      ylab <- xlab <- ""
      yaxt <- xaxt <- "n"
      if (i %in% (1 + nc * (0:nr))) {
        l <- 1
        ylab = "Density"
        yaxt = "s"
      }
      if (i > nc * (nr - 1)) {
        b <- 1
        xlab = "X"
        xaxt = "s"
      }
      par(mai = c(b, l, t, r))
      plot(-10, -10, xlab = xlab, xaxt = xaxt, ylab = ylab, yaxt = yaxt, ylim = c(0, 
        max.den), xlim = range(x.loc), main = "", cex.lab = 1.7, cex.axis = 1.7)
      grid()
      lines(x.loc2, keep.den[, i + 1])
      legend("top", legend = paste("t=", round(record[i], 2), sep = ""), col = "white", 
        pch = 1, bty = "n")
      if (!is.null(low.thresh)) 
        if (is.numeric(low.thresh) & !is.na(low.thresh[i])) {
          if (any(dat[[i]] <= low.thresh[i])) {
            legend("left", title = expression("Pr(x<" * tau[L] * ")"), bty = "n", 
            legend = round(mean(dat[[i]] <= low.thresh[i]), 3), pch = 1, 
            col = "white", cex = 1.7)
          }
        }
      if (!is.null(hi.thresh)) 
        if (is.numeric(hi.thresh) & !is.na(hi.thresh[i])) {
          if (any(dat[[i]] >= hi.thresh[i])) {
            legend("right", title = expression("Pr(x>" * tau[U] * ")"), bty = "n", 
            legend = round(mean(dat[[i]] >= hi.thresh[i]), 3), pch = 1, col = "white", 
            cex = 1.7)
          }
        }
    }
  }
  rel.ent <- rel.entropy(dat, bw=bw, hi.thresh = hi.thresh, low.thresh = low.thresh)

  boot.dat <- dat

## use parallel processing to speed up bootstrapping
  try(registerDoMC(cores = 6))
  hat <- foreach(boot.it = 1:n.boot, .combine = rbind) %dopar% {
#    for(boot.it in 1:20){
   if (boot.it %in% c(20, 50, 100, 200, 500, 1000, 2000, 5000))
      print(paste("Iteration", boot.it))

# set seed for random number generation, which determines bootstrap sample
    set.seed(boot.it)

# Construct bootstrapped sample
    boot.dat <- dat
    if (is.null(dat.mat)) {
      for (i in 1:n.t) boot.dat[[i]] <- sample(dat[[i]], length(dat[[i]]), 
        replace = TRUE)
      if (boot.it == 1) ## Use the original data during the first iteration
        boot.dat <- dat
    }

# Use this code when bootstrapping from image (uses variability across columns in place of assumed binomial variability)
    if (!is.null(dat.mat)) {  
      for (i in 1:n.t) {
        boot.col <- sample(1:ncol(dat.mat[[i]]), ncol(dat.mat[[i]]), replace = TRUE)
        boot.cnt <- round(mat.fun(dat.mat[[i]][, boot.col]))
        if (boot.it == 1) 
          boot.cnt <- mat.fun(dat.mat[[i]])
        boot.cnt[boot.cnt < 0] <- 0
        boot.dat[[i]] <- rep(as.numeric(rownames(dat.mat[[i]])), boot.cnt)
      }
    }

### Estimate required components from bootstrapped sample
    parms <- DnA.ingredients(x.loc, record, boot.dat, Plot = FALSE,
                            hi.thresh = hi.thresh, low.thresh = low.thresh,bw=bw,boot.it=boot.it)
    dP_dt <- parms$dP_dt
    dp_dx <- parms$dp_dx
    p <- parms$p

### To unnormalize estimated distributions (may be useful for truncated data)
   if(unnorm){
     dP_dt<-scale(dP_dt,center=FALSE,scale=1/(sapply(dat,length)))
     dp_dx<-scale(dp_dx,center=FALSE,scale=1/(sapply(dat,length)))
     p<-scale(p,center=FALSE,scale=1/(sapply(dat,length)))
   }

    dlnp_dx <- dp_dx/p
    dP_p <- dP_dt/p

### Estimate (modified) drift and diffusion for each time pair
    nm <- t.a.hat <- t.d.hat <- NULL
    for (t1 in 1:(length(record) - 1)) {  ## Time 1
      for (t2 in (t1 + 1):length(record)) {  ## Time 2
        nm <- c(nm, paste("T1_", t1, "T2_", t2))
        t.d.hat <- cbind(t.d.hat, (p[, t2] * dP_dt[, t1] - p[, t1] * dP_dt[, 
          t2])/(p[, t2] * dp_dx[, t1] - p[, t1] * dp_dx[, t2]))
        t.a.hat <- cbind(t.a.hat, dlnp_dx[, t1] * t.d.hat[, ncol(t.d.hat)] - 
          dP_p[, t1])
      }
    }
    colnames(t.a.hat) <- paste("A", nm)
    colnames(t.d.hat) <- paste("D", nm)
    rownames(t.a.hat) <- rownames(t.d.hat) <- paste("X", x.loc, "X", sep = "")
    signif(cbind(t.d.hat, t.a.hat), 3)
  }  ### end parallel

### Get plots for original data if requested
  if(plot.den){
    parms <- DnA.ingredients(x.loc, record, dat, Plot = TRUE,
                             hi.thresh = hi.thresh, low.thresh = low.thresh,bw=bw,boot.it=100)
  }

### Estimate relative entropy between each pair of times
  rel.ent <- nm <- NULL
  for (t1 in 1:(length(record) - 1)) {
    for (t2 in (t1 + 1):length(record)) {
      nm <- c(nm, paste("T1_", t1, "T2_", t2))
      rel.ent <- c(rel.ent, rel.entropy(dat[c(t1, t2)], bw=bw[c(t1, t2)], hi.thresh = hi.thresh[c(t1, 
        t2)], low.thresh = low.thresh[c(t1, t2)])[1])
    }
  }
  names(rel.ent) <- nm
  return(list(d.hat = hat[, grep("D", colnames(hat))], a.hat = hat[, grep("A", 
    colnames(hat))], rel.ent = rel.ent, x.loc = x.loc, record = record, n.boot = n.boot, 
    Adj = Adj, keep.den = keep.den,bw=bw))
}


























DnA.results <- function(obj, method = "weighted_median", ,smooth="poly", Adj = 1, plot.res = FALSE, df.spline.D = NULL, df.spline.A = NULL){

# 'obj' is the object returned from DnA.analyze
# 'method' describes how to combine estimates across time pairs.  Supported values are "mean", "median", "weighted_median" (default), and "mode"  
#'smooth' describes how to smooth and interpolate drift and diffusion estimates across values of obj$x.loc.  Supported values are 'poly' (default) and 'spline'.  'poly' is recommended as it is the only supported method that incorporates uncertainty from the derivative of the diffusion function into the uncertainty of the drift profile.  
#'plot.res' is an indicator.  If TRUE, the estimated drift and a diffusion profiles are displayed. 
#'df.spline.D' is an optional constraint for the degrees of freedom used when fitting a spline to the pointwise diffusion estimates  (ignored if smooth!='spline') 
#'df.spline.A' is an optional constraint for the degrees of freedom used when fitting a spline to the pointwise drift estimates  (ignored if smooth!='spline') 

  for (i in 1:length(method)) {
    t.method <- grep(method[i], c("mean", "median", "weighted_median", "mode"), 
      value = TRUE)
    if (length(method[i]) == 0) 
      stop(paste("method[", i, "]=", method[i], " does match any of:\n'mean','median','weighted_median', or 'mode'", 
        sep = ""))
    if (length(method[i]) > 1) 
      stop(paste("method[", i, "]=", method[i], " matches", t.method, ". method[", 
        i, "] must uniquely match one of: 'mean','median','weighted_median', or 'mode'", 
        sep = ""))
    if (sum(grepl(method[i], c("mean", "median", "weighted_median", "mode"))) != 
      1) 
      stop(paste("method[", i, "]=", method[i], ".\n Elements of method must uniquely match one of: 'mean','median','weighted_median', or 'mode'.", 
        sep = ""))
    method[i] <- t.method
  }
  if (length(method) == 1) 
    method <- list(local.d = method, local.a = method)
  if (is.null(method$local.d)) 
    stop("method$local.d is not specified")
  if (is.null(method$local.a)) 
    stop("method$local.a is not specified")
  
  rel.ent = obj$rel.ent
  x.loc = obj$x.loc
  record = obj$record
  n.boot = obj$n.boot
  obj$a.hat[obj$d.hat < 0 | is.na(obj$d.hat)] <- NA
  obj$d.hat[obj$d.hat < 0] <- NA
  
  if (!is.null(df.spline.D)) 
    if (!is.numeric(df.spline.D)) 
      stop(paste("If specified, df.spline.D must be a real number between 1 and length(x.loc)=", 
        length(x.loc)))
  if (!is.null(df.spline.A)) 
    if (!is.numeric(df.spline.A)) 
      stop(paste("If specified, df.spline.A must be a real number between 1 and length(x.loc)=", 
        length(x.loc)))
  
  if (is.null(Adj)) 
    Adj = obj$Adj

     dev.new(width = 10, height = 10)
    layout(matrix(1:2, 2, 1))
   
  for (Func in c("D", "A")) {
    if (Func == "D") {
      log.d.BOOT <- log.d.pair.est <- log.d <- se.log.d <- boot.sd.log.d.pair.est <- NULL
      METHOD <- method$local.d
    }
    if (Func == "A") {
      a.BOOT <- a.pair.est <- a <- se.a <- boot.sd.a.pair.est <- NULL
      METHOD <- method$local.a
    }
    
    for (x in x.loc) {
      boot <- NULL
      if (Func == "D"){ 
        t.dat <- suppressWarnings(log(obj$d.hat[grep(paste("X", x, "X", sep = ""), 
          rownames(obj$d.hat)), ]))
      }
      if (Func == "A"){ 
        t.dat <- suppressWarnings(obj$a.hat[grep(paste("X", x, "X", sep = ""), 
          rownames(obj$a.hat)), ])
      }
 
      
      for(i in 1:ncol(t.dat)) 
      if (grepl("median", METHOD)) {
        w <- rep(1, ncol(t.dat))
        if (METHOD == "weighted_median") 
          w <- 1/apply(t.dat, 2, var, na.rm = TRUE)
        boot <- apply(t.dat, 1, weightedMedian, w = w, na.rm = TRUE)
      }
      if (METHOD == "mean") {
        bad <- which(colMeans(t.dat < 0 | t.dat == Inf) > 0.1)
        t.dat <- apply(t.dat, 2, function(x) {
          suppressWarnings(M <- max(x[x < Inf]))
          x[x == Inf] <- M
          x
        })
        for (i in 1:nrow(t.dat)) {
          if (any(!is.na(t.dat[i, ]))) {
          invVAR <- solve(cov(as.matrix(t.dat[, !is.na(t.dat[i, ])]), use = "pairwise.complete"))
          boot <- c(boot, sum(invVAR %*% t.dat[1, !is.na(t.dat[i, ])])/sum(invVAR))
          } else boot <- c(boot, NA)
        }
      }
      if (METHOD == "mode") {
        for (i in 1:nrow(t.dat)) {
          t.d.hat <- t.dat[i, !is.na(t.dat[i, ])]
          dd <- seq(min(t.d.hat), max(t.d.hat), length.out = 500)
          t.den <- den.fun(dat = t.d.hat, x = dd, deriv = 0, Adj = Adj)
          boot <- c(boot, dd[which.max(t.den)])
        }
      }
      if (Func == "D") {
        log.d.pair.est <- rbind(log.d.pair.est, t.dat[1, ])
        colnames(log.d.pair.est) <- colnames(t.dat)
        log.d <- c(log.d, boot[1])
        se.log.d <- c(se.log.d, sd(boot, na.rm = TRUE))
        log.d.BOOT <- cbind(log.d.BOOT, boot)
        boot.sd.log.d.pair.est <- rbind(boot.sd.log.d.pair.est, apply(t.dat, 2, sd, 
          na.rm = TRUE))
      }
      if (Func == "A") {
        a.pair.est <- rbind(a.pair.est, t.dat[1, ])
        colnames(a.pair.est) <- colnames(t.dat)
        a.BOOT <- cbind(a.BOOT, boot)
        boot.sd.a.pair.est <- rbind(boot.sd.a.pair.est, apply(t.dat, 2, sd, 
          na.rm = TRUE))
      }
    }
  }
  rownames(a.pair.est) <- rownames(log.d.pair.est) <- rownames(boot.sd.a.pair.est) <- rownames(boot.sd.log.d.pair.est) <- x.loc
  if(plot.res){
    dev.new(width = 10, height = 10)
    layout(matrix(1:2, 2, 1))
    par(mai = c(0, 1, 1, 0.1))
    plot(x.loc, exp(log.d), type = "l", lwd = 3, log = "y", ylim = exp(range(c(log.d + 
      2.5 * se.log.d, log.d - 2.5 * se.log.d), na.rm = TRUE)), ylab = expression(hat(D)), 
      xlab = "", xaxt = "n")
    grid(equilogs = FALSE)
    for (i in 1:nrow(log.d.pair.est)) points(rep(x.loc[i], ncol(log.d.pair.est)), exp(log.d.pair.est[i, 
      ]), col = COL, pch = rep(1:20, each = length(COL)))

    lines(x.loc, exp(log.d), lwd = 3)
    lines(x.loc, exp(log.d + 2 * se.log.d), lwd = 3)
    lines(x.loc, exp(log.d - 2 * se.log.d), lwd = 3)
  }
  
  for (Func in c("D", "A")) {

    if (Func == "D") {
      use <- which(!is.na(log.d))
      est <- log.d.BOOT[, use]
      given.df <- df.spline.D
    }
    if (Func == "A") {
      use<- which(!is.na(a))
      est <- a.BOOT[,use]
      given.df <- df.spline.A
    }
    full.est<-est
    se.est <- apply(est, 2, sd)
    nn <- ncol(est)
    if (is.null(given.df)) 
      given.df <- 0
    if (given.df < 2) {
      df.try <- 1
      sig <- cov(est,use="pairwise.complete")
      cov.it<-0
      while(any(is.na(diag(sig)))&&cov.it<5){
        use<-use[-which(is.na(diag(sig)))]
        est<-est[,-which(is.na(diag(sig)))]
        sig<-cov(est,use="pairwise.complete")
        if(cov.it==4){
          print(paste("Problem with missing values in bootstrap covariance matrix for", Func))
          print(paste("Remaining column indices are:",use))
        }
        cov.it<-cov.it+1
      }
      invVAR <- NULL
      try(invVAR <- solve(sig), silent = TRUE)
      ### If covariance matrix is singular, slightly increase diagonal elements
      ### before inverting
      if (is.null(invVAR)) {
        sig <- sig + diag(rep(min(diag(sig))/1000, sum(!is.na(log.d))))
        invVAR <- solve(sig)
      }
      t.mn <- sum(invVAR %*% est[1, ])/sum(invVAR)
      t.fit <- list(x = x.loc[!is.na(log.d)], y = rep(t.mn, nn))
      test <- (est[1, ] - t.mn)
      test <- (t(test) %*% invVAR %*% test * (n.boot - nn)/(nn * (n.boot - 
        1)))
      p.constant <- pf(test, nn, n.boot - nn, lower.tail = FALSE)
      print(paste("P-value for test of constant", Func,"=",signif(t.mn,3), "is", signif(p.constant, 3)))


      if(smooth=="LASSO"){
library(lars)

        Y<-est[1,]
        V<-sig
        V1_2<-chol(V)
        Vneg1_2<-solve(V1_2)
  
        Ymod<-as.vector(Vneg1_2%*%Y)
        poly.x<-poly(x.loc,degree=(length(x.loc)-7))
        
        x.matfun<-function(x){
          xx<-t(t(matrix(x-mean(x.loc),length(x),length(x.loc)-2,byrow=FALSE))^(0:(length(x.loc)-3)))
        }
        x.mat<-x.matfun(x.loc)
        Xmod<-Vneg1_2%*%x.mat
        fit<-lars(x=Xmod,y=Ymod)
        predict(fit,s=c(0:10)/10,mode="fraction",type="coef")

        n<-length(Y)
        coef.mat<-predict(fit,type="coef")$coefficients
        pred.mat<-predict(fit,newx=Xmod,type="fit")$fit

        mse<-colMeans((pred.mat-Ymod)^2)

        ### Only allow one model of each d.f. 
        df<-rowSums(coef.mat!=0)
        use<-NULL
        for(DF in unique(df)){
          use<-c(use,max(which(df==DF)))
        }
        bic<-(mse+log(n)/n*df)[use]
        PR<-exp(-.5*bic)
        PR<-PR/sum(PR)

        pred.mat<-pred.mat[,use]
        x.mat_up<-Vneg1_2%*%x.matfun(x.loc+1e-6)
        x.mat_dwn<-Vneg1_2%*%x.matfun(x.loc-1e-6)

        pred.mat_up<-V1_2%*%predict(fit,newx=x.mat_up,type="fit")$fit[,use]
        pred.mat_dwn<-V1_2%*%predict(fit,newx=x.mat_dwn,type="fit")$fit[,use]
        fit$dpred.dx<-t(pred.mat_up-pred.mat_dwn)/2e-6
        t.fit$x<-x.loc
        t.fit$y<-as.vector(V1_2%*%pred.mat%*%PR)      
      }



      if(smooth=="poly"){
        poly.fit<-function(x.loc,deg,BIC=NULL,dpred.dx=NULL,pred=NULL,eps=1e-6){
            xx<-matrix(rep(1,length(x.loc)),length(x.loc))
            poly.x<-1
            if(deg>0){
              poly.x<-NULL
              try(poly.x<-poly(x.loc,degree=deg),silent=TRUE)
              xx<-cbind(xx,poly.x)
            }
            if(!is.null(poly.x)){
            xtxi <- solve(t(xx) %*% invVAR %*% xx)
            beta <- xtxi %*% t(xx) %*% invVAR %*% est[1,]
            t.pred<-as.vector(xx%*%beta)
            pred<-rbind(pred,t.pred)
            bic<-NULL
            try(bic<--2*dmt(est[1,],mean=t.pred,S=sig,df=n.boot-nn,log=TRUE)+(deg+1)*log(nn),silent=TRUE)
            if(is.null(bic)){
              nu<-n.boot-nn
              resid<-est[1,]-t.pred            
              bic<--2*(lgamma(n.boot/2)-lgamma(nu/2)-nn/2*log(nu*pi)-log(det(sig))/2-n.boot/2*log(1+t(resid)%*%invVAR%*%resid/nu))+(deg+1)*log(nn)
}
            BIC<-c(BIC,bic)
            
            if(Func=="D"){
              if(deg==0){dpred.dx<-rbind(dpred.dx,rep(0,length(x.loc)))}
              else{
                xx1<-cbind(rep(1,length(x.loc)),predict(poly.x,newdata=x.loc+eps))
                xx2<-cbind(rep(1,length(x.loc)),predict(poly.x,newdata=x.loc-eps))
 
                z<-as.vector(exp(xx1%*%beta)-exp(xx2%*%beta))/(2*eps)
                dpred.dx<-rbind(dpred.dx,z)
              }
            }
          } else{
            BIC<-c(BIC,Inf)
            pred<-rbind(pred,rep(0,length(x.loc)))
            dpred.dx <- rbind(dpred.dx,rep(0,length(x.loc)))
          }
             return(list(BIC=BIC,dpred.dx=dpred.dx,pred=pred,xx=xx,beta=beta))
        }
        
        ### Use Baysian model averaging across polynomial fits of different degrees
        ### to incorporate uncertainty in model choice in uncertainty in diffusion
        ### and its derivative
          fit<-NULL
#          plot(x.loc[use],est[1,],ylim=c(0,2))
#          for(i in 1:length(use))
#            lines(rep(x.loc[use[i]],2),est[1,i]+c(-2,2)*sqrt(diag(sig)[i]))
        options(show.error.messages=FALSE)
          for(deg in 0:(min(10,length(x.loc[use])-2))){
            fit<-poly.fit(x.loc=x.loc[use],deg=deg,BIC=fit$BIC,dpred.dx=fit$dpred.dx,pred=fit$pred)
#            lines(x.loc[use],as.vector(fit$xx%*%fit$beta),col=deg+1)
          }
          PR<-exp(-.5*fit$BIC)
          PR<-PR/sum(PR)
          round(PR,2)
          while(PR[length(PR)]>max(PR)/1000&deg<(length(x.loc[use])-2)){
            deg<-deg+1
            fit<-poly.fit(x.loc=x.loc[use],deg=deg,BIC=fit$BIC,dpred.dx=fit$dpred.dx,pred=fit$pred)
#            lines(x.loc[use],as.vector(fit$xx%*%fit$beta),col=deg+1)
            PR<-exp(-.5*fit$BIC)
            PR<-PR/sum(PR)
          }
          options(show.error.messages=TRUE)
          t.fit$x<-x.loc[use]
          t.fit$y<-as.vector(t(fit$pred)%*%PR)
           } 
      
      }
    if (Func == "D") {
      smooth.log.d <- t.fit
      d.mod.prob<-if(smooth=="poly") PR else NA
      d.p.const<-p.constant
     if(plot.res) lines(smooth.log.d$x, exp(smooth.log.d$y), col = 3, lwd = 3)
      ### Adjust drift estimates for derivative of diffusion function
      dpred.dx<-matrix(NA,length(PR),length(x.loc))
      dpred.dx[,use]<-fit$dpred.dx
      a.BOOT[1,]<-a.BOOT[1,]+as.vector(t(dpred.dx)%*%PR)
      mod.ind<-sample(1:length(PR),10*n.boot,replace=TRUE,prob=PR)
       a.BOOT <- rbind(a.BOOT[1,],a.BOOT[rep(1:n.boot,10),]+dpred.dx[mod.ind,])
      a<-a.BOOT[1,]
      se.a<-apply(a.BOOT,2,sd)
    }

    if (Func == "A"){ 
      smooth.a <- t.fit
      a.mod.prob<-if(smooth=="poly") PR else NA
      a.p.const<-p.constant
    }
    }

  u.BOOT<--scale(a.BOOT[,use][,-1]+a.BOOT[,use][,-length(use)],center=FALSE,scale=2/diff(x.loc[use]))
  u.BOOT <- cbind(rep(0,nrow(u.BOOT)),t(apply(u.BOOT,1,cumsum)))
  se.u <- apply(u.BOOT, 2, sd)
  u <- u.BOOT[1, ] - median(u.BOOT[1, ], na.rm = TRUE)
if(plot.res){  
  par(mai = c(1, 1, 0, 0.1))
  plot(x.loc, a, type = "l", lwd = 3, ylim = quantile(c(a + 2.5 * se.a, a - 2.5 * 
    se.a),c(.2,.8), na.rm = TRUE), ylab = expression(hat(A)), xlab = "X", main = "")
  grid()
  for (i in 1:nrow(a.pair.est)) points(rep(x.loc[i], ncol(a.pair.est)), a.pair.est[i, 
    ], col = COL, pch = rep(1:20, each = length(COL)))
  lines(x.loc, a, lwd = 3)
  lines(x.loc, a + 2 * se.a, lwd = 3)
  lines(x.loc, a - 2 * se.a, lwd = 3)
  lines(smooth.a$x, smooth.a$y, col = 3, lwd = 3)
 } 

  return(list( x.loc = x.loc, u = u, se.u = se.u,
    a.pair.est = a.pair.est, boot.sd.a.pair.est = boot.sd.a.pair.est, a = a, se.a = se.a, smooth.a =smooth.a, a.mod.prob=a.mod.prob,a.p.const=a.p.const,  
    log.d.pair.est = log.d.pair.est, boot.sd.log.d.pair.est = boot.sd.log.d.pair.est, log.d = log.d, se.log.d = se.log.d, smooth.log.d = smooth.log.d,d.mod.prob=d.mod.prob,d.p.const=d.p.const
    ))
}

