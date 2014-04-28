


### Loop across different choices of diffusion profiles
for(DIFF in c("Constant")){ # subset of "Nonconstant" and "Constant")){

## Function describing diffusion profile
sim.D <- function(x, t) {
    if(DIFF!="Nonconstant") z<-rep(1,length(x))
    if(DIFF=="Nonconstant") z<-0.5 + 2*exp(x)/(1+exp(x))
    z
}

### Loop across different choices of initial distributions
for(INIT in c("3peak")){# subset of "Unif", "1peak", "2peak","3peak", and "Broad"

## How many entities should be in simulated population 
n.part <- 100000

## Factor for scaling bandwidth in density estimation
Adj <- 1

### What times should the evolving distribution be observed
record <- c(0,.01,.02,.05,.1,.2,.5,1); note<-"8times"  ## This is our usual setting
#record<-round(seq(0,1,length.out=8),2);note<-"EvenTime"  ## Used in the Supporting Information

## Discrete time step using during simulation
t.step <- 2e-04

## Total time interval to simulate over
n.days <- max(record)

## Should extra noise be added to data
blur <- 0
blur.fun <- function(x) {
    SD <- exp(-x/(abs(x) + 4))/2
    rnorm(length(x), x, SD)
}

## Names of the different landscape potentials used in simulations
poten<-c("ramp","parabola","bigbarrier","littlebarrier","lowleft","periodicramp")
NAME<-c("Ramp","Parabolic","Large Barrier","Small Barrier","Asymmetric","Periodic Ramp")
names(NAME)<-poten

## Initialize lists where simulation results are stored
sim.res2<-sim.res <- NULL

## Loop across potentials 
for(pot in poten){
print(pot)

### A is the function defining the drift
if(pot=="ramp") A <- function(x) rep(1,length(x))
if(pot=="parabola") A <- function(x) -x
if(pot%in%c("bigbarrier","littlebarrier","lowleft"))
  A <- function(x){
	sig<-c(1,1)
	p<-c( .5, .5 )
        mu<-c(-1.75,1.75)
      if(pot=="lowleft")p <-c(.8,.2) 
      if(pot=="bigbarrier") sig<-c(.75,.75)
        denom<-0;numer<-0
	for(i in 1:length(sig)){ 
		mag<-p[i]*dnorm(x,mu[i],sig[i])
		denom<-denom+mag
		numer<-numer+mag*(mu[i]-x)/sig[i]^2
	}
	z<-numer/denom
	z
  }
if(pot=="periodicramp") A <- function(x) 1+cos(2.5*x)
    if (!"seed" %in% ls()) seed <- 1
    set.seed(seed)
    print(paste("seed=", seed))
    seed <- seed + 1

### Initialize population distribution
if(INIT=="1peak") x <-rgamma(n.part,2,2)-2
if(INIT=="2peak")x <- rnorm(n.part,c(-2,2),.1)
if(INIT=="3peak")x <-rgamma(n.part,2,2)+c(-2,0,2) 
if(INIT=="Unif") x<-runif(n.part,-2,2)
if(INIT=="Broad")x<-rnorm(n.part,-1,3)

### day tracks time during simulation
    day <- 0

### sim.dist is used to record observations at chosen times
    sim.dist <- NULL

### Begin Langevin simulation
    for (i in 1:round(n.days/t.step)) {
       
### If current time is one of the chosen recording times, then save the current values of the population
       if (day %in% record) {
           sim.dist <- cbind(sim.dist, sample(round(x, 3),length(x)))
           colnames(sim.dist)[ncol(sim.dist)] <- day
       }

### Simulate diffusion step size
       diffuse <- sqrt(2 * sim.D(x, day)* t.step) * rnorm(n.part)

### Take step
       x <- x + A(x) * t.step + diffuse

### Update time tracker
       day <- round(day + t.step, 4)

### If simulation is ending, record final population status
       if (day ==max(record)) {
           sim.dist <- cbind(sim.dist,  sample(round(x, 3),length(x)))

           colnames(sim.dist)[ncol(sim.dist)] <- day
       }
   }
  
### If chosen, blur the population distribution by adding noise
if (blur) 
       for (i in 1:ncol(sim.dist)) sim.dist[, i] <- round(blur.fun(sim.dist[, i]), 3)

## Identify where the tails (lower 0-2% and upper 98-100%) of the population begin
quants<-sapply(sim.dist2,quantile,c(.02,.98))

## Estimate the drift and diffusion profiles at 30 locations within the center of the distributions
x.loc<-signif(seq(max(quants[1,]),min(quants[2,]),length.out=30),3)

## No censoring used in the simulations
hi.thresh<-low.thresh<-NULL

### Change observed data to a list instead of a matrix
sim.dist2<-list()
for(i in 1:ncol(sim.dist)) sim.dist2[[i]]<-sim.dist[,i]

### Close all current graphics
while(dev.cur()>1)dev.off()

### Load functions used for analysis
source("DiffusionAnalysisFunctions12.R")

### Analyze simulated distributions
## Estimate (modified) drift and diffusion at each value of x.loc for each pair of observation times
sim.res[[pot]] <- DnA.analyze(dat=sim.dist2, n.boot = 500, record = record, x.loc = x.loc, Adj = 1.75, hi.thresh=hi.thresh,low.thresh=low.thresh,bw.method="SJ")

### Aggregate time pair specific estimates, smooth and interpolate to obtain final estimates of drift and diffusion profiles
sim.res2[[pot]]<-DnA.results(sim.res[[pot]],plot.res=FALSE)
}

### Save results
save(sim.res,file=paste("SimulationPreresultsN",n.part,DIFF,"Diff",INIT,"Init",note,".R",sep=""))
save(sim.res2,file=paste("SimulationResultsN",n.part,DIFF,"Diff",INIT,"Init",note,".R",sep=""))
}
}


#################################################################
## Plot simulation results (e.g. Figure 4 in main manuscript) ##
################################################################


### Loop across different choices of diffusion profiles
for(DIFF in c("Constant")){ # subset of "Nonconstant" and "Constant")){

## Function describing diffusion profile
sim.D <- function(x, t) {
    if(DIFF!="Nonconstant") z<-rep(1,length(x))
    if(DIFF=="Nonconstant") z<-0.5 + 2*exp(x)/(1+exp(x))
    z
}

### Loop across different choices of initial distributions
for(INIT in c("3peak")){# subset of "Unif", "1peak", "2peak","3peak", and "Broad"

n.part <- 100000

### Load simulations results
load(file=paste("SimulationPreresultsN",n.part,DIFF,"Diff",INIT,"Init",note,".R",sep=""))
load(file=paste("SimulationResultsN",n.part,DIFF,"Diff",INIT,"Init",note,".R",sep=""))

### Initialize .eps file to save results
postscript(paste("SimulationPlotsN",n.part,DIFF,"Diff",INIT,"Init",note,"2.eps",sep=""), horizontal=FALSE,onefile=FALSE,paper="special",height=9,width=12)
#dev.new(width = 12, height = 9)

### Format graphics display
nc<-6;nr<-4
                nf <- layout(matrix(1:(nc * nr), nr, nc, byrow = FALSE), widths = c(1.15 + 
                  0.05 * nc, rep(1, nc - 2),1.1), heights = c(1.1+.05*nr,rep(.8, nr - 2), 1.1 + 
                  0.05 * nr))
COL<-c("grey","green","orange")

## Loop across potentials 
for(pot in poten){
if(pot=="ramp") A <- function(x) rep(1,length(x))
if(pot=="parabola") A <- function(x) -x
if(pot%in%c("bigbarrier","littlebarrier","lowleft"))
  A <- function(x){
	sig<-c(1,1)
	p<-c( .5, .5 )
        mu<-c(-1.75,1.75)
      if(pot=="lowleft")p <-c(.8,.2) 
      if(pot=="bigbarrier") sig<-c(.75,.75)
        denom<-0;numer<-0
	for(i in 1:length(sig)){ 
		mag<-p[i]*dnorm(x,mu[i],sig[i])
		denom<-denom+mag
		numer<-numer+mag*(mu[i]-x)/sig[i]^2
	}
	z<-numer/denom
	z
  }
if(pot=="periodicramp") A <- function(x) 1+cos(2.5*x)

### Get components from simulation results required for plotting 
record<-sim.res[[pot]]$record
x.loc<-sim.res2[[pot]]$x.loc
t.x<-seq(min(x.loc),max(x.loc),length.out=500)
true.a<-A(t.x)
true.u<--cumsum(true.a)*mean(diff(t.x))
true.u<-true.u-median(true.u)
a<-sim.res2[[pot]]$a
se.a<-sim.res2[[pot]]$se.a
u<-sim.res2[[pot]]$u
u<-u+(median(true.u)-median(u))#(true.u[1]-u[1])
se.u<-sim.res2[[pot]]$se.u
smooth.a<-sim.res2[[pot]]$smooth.a
pred.a<-smooth.a$y
est.u<--c(0,cumsum((pred.a[-1]+pred.a[-length(pred.a)])/2*diff(smooth.a$x)))
est.u<-est.u+(median(true.u)-median(est.u))#(true.u[1]-est.u[1])
d.p.const<-sim.res2[[pot]]$d.p.const

### Format graphics display
l <- r <- 0
ylab <- xlab <- ""
yaxt <- xaxt <- "n"
if(pot=="ramp"){  ylab="f(x,t)"; yaxt="s";l<-.7;YLIM1<-c(-4,4)}#range(true.u)}
if(pot=="periodicramp"){ r<-.2}
cex.fac=2.7
par(mai = c(0, l, .7, r))
den<-sim.res[[pot]]$keep.den
ylim=c(0,.75)
if(INIT=="3peak") ylim=c(0,.35)

### Plot estimated density from 4 of the observation times
plot(-10,-10,ylab=ylab,xlim=range(den[,1]),xaxt="n",xlab="",yaxt="n",main=NAME[pot],ylim=ylim,cex.main=2.3,cex.lab=cex.fac,cex.axis=cex.fac,type="l",lwd=2)
if(pot=="ramp"){
  axis(2,at=(0:7)/20,labels=rep("",8))
  axis(2,at=c(0,.15,.3),labels=c("0","0.15","0.3"),tick=FALSE,cex.axis=cex.fac)
}
grid()
use<-c(1,3,5,length(record))
den.col<-c("red","orange","green","blue")#makeTransparent(c("red","orange","green","blue"),140)
#for(i in 0:(length(use)-1)) den.col<-c(den.col,rgb(i/(length(use)-1),0,1-i/(length(use)-1)))
for(i in 1:length(use)) lines(den[,1],den[,use[i]+1],col=den.col[i],lwd=2)
if(pot=="ramp") legend("top",legend=paste("t=",record[use],sep=""),ncol=2,col=den.col,lwd=2,cex=1.6,bty="n",seg.len=1)


### Plot true and estimated landscape potential
if(pot=="ramp")   ylab="Potential"
par(mai = c(0, l, 0, r))
plot(t.x,true.u,type="l",col=COL[3],lwd=2,ylim=YLIM1,ylab=ylab,xaxt="n",xlab="",main="",yaxt=yaxt,cex.lab=cex.fac,cex.axis=cex.fac)
grid()
lines(smooth.a$x,est.u,lwd=2,col=COL[2])
lines(smooth.a$x, u, lwd = 2,col=COL[1])
lines(smooth.a$x, u + 2 * se.u, lwd = 1,col=COL[1])
lines(smooth.a$x, u - 2 * se.u, lwd = 1,col=COL[1])
#if(pot=="ramp") legend("bottomleft",legend=c("Truth","Smooth Est","Pt Est"),col=COL[3:1],lwd=2,cex=1.65,bty="n",seg.len=1)

### Plot true and estimated drift profiles
if(nr==4){
if(pot=="ramp")  ylab="Drift"
plot(t.x,true.a,type="l",col=COL[3],lwd=2,ylab=ylab,xaxt="n",xlab="",yaxt=yaxt,main="",ylim=c(-8,8),
     #ylim=c(-2.3,2.3),
     cex.lab=cex.fac,cex.axis=cex.fac)
grid()
lines(x.loc, a, lwd = 2,col=COL[1])
lines(x.loc, a + 2 * se.a, lwd = 1,col=COL[1])
lines(x.loc, a - 2 * se.a, lwd = 1,col=COL[1])
lines(smooth.a$x,smooth.a$y,lwd=2,col=COL[2])
if(pot=="ramp") legend("bottomleft",legend=c("Truth","Bay Mod Ave","Pt Est"),col=COL[3:1],lwd=2,cex=2,bty="n",seg.len=1)
}

### Plot true and estimated diffusion profiles
par(mai = c(.7, l, 0, r))
if(pot=="ramp") ylab<-"Diffusion"
XAXT="s"
if(INIT=="3peak") XAXT<-"n"
plot(t.x,sim.D(t.x,1),ylim=c(.1,10),log="y",ylab=ylab,xlab="x",xaxt=XAXT,yaxt="n",main="",type="l",col=COL[3],lwd=2,cex.lab=cex.fac,cex.axis=cex.fac)
if(INIT=="3peak"){
  axis(1,at=c(-4:4),labels=rep("",9))
  axis(1,at=c(-1,1,3),labels=c(-1,1,3),cex.axis=cex.fac,line=.7,tick=FALSE)
} 
grid(equilogs=FALSE)
if(pot=="ramp"){
  axis(2,at=c(.1,.2,.5,1,2,5,10),labels=rep("",7),cex.axis=cex.fac)
  axis(2,at=c(.3,3),labels=c("0.3","3"),cex.axis=cex.fac)
  axis(2,at=c(1),labels=c("1"),cex.axis=cex.fac)
  axis(2,at=rep(1:9,3)*10^rep(-2:0,each=9),label=rep("",27))
}
log.d<-sim.res2[[pot]]$log.d
se.log.d<-sim.res2[[pot]]$se.log.d
  smooth.log.d<-sim.res2[[pot]]$smooth.log.d
    lines(x.loc, exp(log.d), lwd = 2,col=COL[1])
    lines(x.loc, exp(log.d + 2 * se.log.d), lwd = 1,col=COL[1])
    lines(x.loc, exp(log.d - 2 * se.log.d), lwd = 1,col=COL[1])
lines(smooth.log.d$x,exp(smooth.log.d$y),lwd=2,col=COL[2])
text(mean(range(t.x)),6,paste("p=",signif(d.p.const,2),sep=""),cex=cex.fac)
}

#dev.copy2pdf(file=paste("SimulationPlotsN",n.part,DIFF,"Diff",INIT,"Init",note,".pdf",sep=""))
dev.off()
}
}



### Plot relative entropy (e.g. Figure 5 in main manuscript)
load(file="SimulationPreresultsN1e+05ConstantDiff3peakInit8times.R")
load(file="SimulationResultsN1e+05ConstantDiff3peakInit8times.R")
r<-sim.res[["bigbarrier"]]$rel.ent
r<-c(r[grep(8,names(r))],0)
r<-r/max(r)

dev.new(height=3,width=3)
#postscript("RelEntropy1e+05ConstantDiff3peakInit8timesbigbarrier.eps", horizontal=FALSE,onefile=FALSE,paper="special",height=3,width=3)
par(mai=c(.7,.7,.1,.1))
plot(sim.res[["bigbarrier"]]$record+.1,r,log="x",xlab="",ylab="",cex.lab=1.5,cex.axis=1.5,cex=1.5,lwd=2,xaxt="n")
axis(1,at=(0:10)/10+.1,labels=rep("",11),cex.axis=1.5)
axis(1,at=(0:5)/5+.1,labels=(0:5)/5,cex.axis=1.5)
mtext("time",1,line=2.2,cex=1.5)
mtext("% Relative Entropy",2.2,line=2,cex=1.5)
ob.time<-approxfun(x=r,y=sim.res[["bigbarrier"]]$record)
t.x<-seq(0,1,length.out=200)
lines(ob.time(t.x)+.1,t.x)
for(i in 0:4){
  lines(rep(ob.time(i/4)+.1,2),c(0,i/4),col=i+1,lwd=2)
  lines(c(0,ob.time(i/4))+.1,c(i/4,i/4),col=i+1,lwd=2,lty=2)
}
#dev.copy2pdf(file="RelEntropy1e+05ConstantDiff3peakInit8timesbigbarrier.pdf")
dev.off()

