#setwd("Y:/Diffusion FP/lipids/New Data")
#dir<-"Y:/Diffusion FP/lipids/Raw Data/"; setwd(dir)
dir<-"/home/spl/Diffusion FP/lipids/Raw Data/"
files<-grep(".txt",list.files(path=dir),value=TRUE)
files<-grep("Equil",files,invert=TRUE,value=TRUE)
Rep<-as.numeric(substr(files,4,4))
time<-substr(files,7,8)
Time<-c(0,6+1/6,12+2/3,23,28.25);names(Time)<-unique(time)
Time<-Time[time]

### first row corresponds to left edge of image
### first column corresponds to bottom edge of image
dat<-vector("list",length(files))
for(i in 1:length(dat)){ 
 t.dat<-scan(paste(dir,files[i],sep=""),nlines=1)
 dat[[i]]<-matrix(scan(paste(dir,files[i],sep="")),ncol=length(t.dat),byrow=TRUE)
}


TRIM<-matrix(c(14,10,1,0,
8,16,1,0,
8,16,1,0,
8,16,1,0,
8,16,1,0,
14,10,1,0,
12,12,1,0,
12,12,1,0,
12,12,1,0,
12,12,1,0,
13,11,1,0,
12,12,1,0,
11,13,1,0,
11,13,1,0,
11,13,1,0),15,4,byrow=TRUE)

trim.dat<-dat
for(i in 1:15)trim.dat[[i]]<-dat[[i]][-c(1:TRIM[i,1],nrow(dat[[i]])-0:TRIM[i,2]),-c(1,ncol(dat[[i]]))]

cutoff<-NULL
for(rep in 1:3){
  pval<-j<-1
  while(pval>.05){
    j<-j+1
    pval<-t.test(trim.dat[Rep==rep&time=="05"][[1]][1,]-trim.dat[Rep==rep&time=="05"][[1]][j,],alternative="less")$p.value
  }
  cutoff<-c(cutoff,j-1)
  print(mean(trim.dat[Rep==rep&time=="05"][[1]][1:(j-1),]))
   print(range(trim.dat[Rep==rep&time=="05"][[1]][1:(j-1),]))
  print(range(rowSums(trim.dat[Rep==rep&time=="05"][[1]][1:(j-1),])))
}

if(0){
x11(width=12);
par(mfrow=c(3,5),mai=rep(.1,4))
for(i in 1:15)image(t(trim.dat[[i]]),zlim=range(sapply(trim.dat,range)),xaxt="n",yaxt="n")

x11(width=12)
par(mfrow=c(3,5),mai=rep(.1,4)) 
for(i in 1:15){ plot(rowMeans(trim.dat[[i]]),ylim=c(30,210));grid()}

x11()
 plot(-10,-10,xlim=range(Time),ylim=c(80,170),ylab="Column Means",xlab="Time")
for(i in 11:65){
 lines(Time[1:5],sapply(dat[1:5],function(x)mean(x[-c(1:20,nrow(x)-0:19),i])))
 lines(Time[1:5],sapply(dat[6:10],function(x)mean(x[-c(1:20,nrow(x)-0:19),i])),col=2)
 lines(Time[1:5],sapply(dat[11:15],function(x)mean(x[-c(1:20,nrow(x)-0:19),i])),col=3)
}
x11()
plot(-10,-10,xlim=c(40,220),ylim=c(0,800),ylab="Row Variance",xlab="Row Mean")
for(i in 1:15) points(rowMeans(dat[[i]][-c(1:20,nrow(dat[[i]])-0:19),]),apply(dat[[i]][-c(1:20,nrow(dat[[i]])-0:19),],1,var),col=as.numeric(time[i]),pch=Rep[i])
}


### remove background from dark current
mod.dat<-lapply(trim.dat,function(x)x-38)

total.inten<-sapply(lapply(lapply(trim.dat,function(x)x-38),rowSums),function(x){x[x<0]<-0; sum(x)})
#dev.new(height=4)
postscript("TotalIntensity.eps", horizontal=FALSE,onefile=FALSE,paper="special",height=4,width=7)
par(mai=c(1,1,.1,.1))
plot(Time,total.inten,ylab="Total Intensity",xlab="Time (min)",col=rep(2:4,each=5),lwd=2,cex=1.7,cex.lab=1.7,cex.axis=1.7)
grid()
legend("topright",legend=paste("Region",1:3),pch=1,pt.lwd=2,pt.cex=1.7,cex=1.7,col=2:4,bg="white")
dev.off()

### Find first row where observations go above saturation or quenching threshold
y.thresh<-8850
y.thresh2<-6000
thresh<-suppressWarnings(sapply(mod.dat,function(x)min(which(rowSums(x)>y.thresh))))
thresh2<-suppressWarnings(sapply(mod.dat,function(x)min(which(rowSums(x)>y.thresh2))))
#sapply(mod.dat,function(x)range(rowSums(x)))


### add observations at high values to normalize for saturation or quenching effects (i.e. same number of observations from each time point)
for(i in 0:2) for(j in 2:5) mod.dat[[i*5+j]][ nrow( mod.dat[[i*5+j]]),]<-mod.dat[[i*5+j]][nrow( mod.dat[[i*5+j]]),]+colSums(mod.dat[[i*5+1]])-colSums(mod.dat[[i*5+j]])


### Use row sums to create observation vectors
tab.dat<-vector("list",length(files))
Dat<-NULL
for(i in 1:length(dat)) Dat<-cbind(Dat,rowSums(mod.dat[[i]]))
Dat[Dat<0]<-0

files<-grep("ROI",grep("Equilibrm",list.files(),value=TRUE),value=TRUE)
station.dat<-trim.dat
for(i in 1:3) station.dat[[15+i]]<-as.matrix(read.table(file=files[i]))
station.dat<-station.dat[c(1:5,16,6:10,17,11:15,18)]

dev.new(height=6,width=10.5)
x.convert=2.14
#postscript("LipidImages.eps", horizontal=FALSE,onefile=FALSE,paper="special",height=6,width=10.5)
zlim<-c(30,237) #range(sapply(trim.dat,range))
nf<-matrix(1:18,3,6,byrow=TRUE)
nf<-cbind(nf,rep(19,3))

layout(nf,widths=c(1.375,rep(1,5),.65),heights=c(1.35,1,1.375))
for(i in 1:18){
  b<-l<-t<-r<-0
  ylab=""
   if(i %in% 1:6)  t<-.75
  if(i %in% 13:18)  b<-.75
  if(i %in% c(1,7,13))  l<-.75
  par(mai=c(b,l,t,r))
  image(station.dat[[i]],xaxt="n",yaxt="n",ylab="",xlab="",zlim=zlim,col=heat.colors(40))
  if(i %in% 1:5) mtext(round(Time[i],2),3,line=.35,cex=2)
  if(i ==6) mtext(58,3,line=.35,cex=2)
  if(i==15) mtext(expression("Distance from Left ("*mu*"m)"),1,line=3.65,cex=2,adj=0)
  if(i %in% c(1,7,13))   mtext(paste("Region",(i+5)/6),2,line=1,cex=2)
  if(i %in% c(13:18,6,11)){
     axis(side=1,at=c(0,80,160)/80/x.convert,labels=c(0,"",""),cex.axis=2.3)
     axis(side=1,at=c(40,120)/80/x.convert,labels=c(40,120),cex.axis=2.3)
   }
  if(i==3) mtext("                  Time (min)",3,line=2.8,cex=2)#,adj=70)
#  if(i==11)
#    {
#      lines(c(.05,.05),c(.1,.683),col=1,lwd=3)
#      text(.25,.15,expression(100*mu*"m"),cex=1.75)
#    }
}
par(bty="n");par(mai=c(.75,0,.75,.01))
plot(-10,-10,ylim=c(0,1),xlim=0:1,ylab="",yaxt="n",xlab="",xaxt="n")
mtext("Inten.",1,line=1.3,cex=2)
require(fields)
image.plot(zlim=zlim,legend.only=TRUE,smallplot=c(.1,.4,.125,.9),col=heat.colors(40),axis.args=list(cex.axis=2))
par(bty="o")
#dev.copy2pdf(file="LipidImages.pdf")
dev.off()

norm.fac<-NULL
for(j in 1:5)
norm.fac<-rbind(norm.fac,sapply(trim.dat[5*(0:2)+j],mean)/mean(sapply(trim.dat[5*(0:2)+j],mean)))
norm.fac<-colMeans(norm.fac)

norm.trim.dat<-trim.dat
#for(r in 1:length(trim.dat)) norm.trim.dat[[r]]<-norm.trim.dat[[r]]/norm.fac[Rep[r]] 

combined.norm.trim.dat<-NULL
for(j in 1:5) combined.norm.trim.dat[[j]]<-cbind(norm.trim.dat[[j]],norm.trim.dat[[5+j]],norm.trim.dat[[10+j]])


backsub.norm.trim.dat<-lapply(norm.trim.dat,function(x)x-38)
bb<-backsub.norm.trim.dat
for(i in 0:2){
  for(j in 1:5){
    backsub.norm.trim.dat[[5*i+j]][ nrow( bb[[5*i+j]]),]<-bb[[5*i+j]][nrow( bb[[5*i+j]]),]+colSums(bb[[5*i+1]])-colSums(bb[[5*i+j]])
    rownames(backsub.norm.trim.dat[[5*i+j]])<-2.14*(1:nrow(bb[[5*i+j]]))
  }
}


rowMedians<-function(x,...) apply(x,1,median,...)

backsub.combined.norm.trim.dat<-lapply(combined.norm.trim.dat,function(x)x-38)
bb<-backsub.combined.norm.trim.dat
for(j in 1:5){
   backsub.combined.norm.trim.dat[[j]][ nrow( bb[[j]]),]<-bb[[j]][nrow( bb[[j]]),]+colSums(bb[[1]])-colSums(bb[[j]])
   rownames(backsub.combined.norm.trim.dat[[j]])<-2.14*(1:nrow(bb[[j]]))
 }
mat.fun<-function(x) rowSums(x/10)
mat.fun2<-function(x) rowMeans(x)
DAT<-round(sapply(backsub.combined.norm.trim.dat,mat.fun))
DAT2<-round(sapply(backsub.combined.norm.trim.dat,mat.fun2))
#DAT<-round(sapply(backsub.norm.trim.dat,mat.fun))
 
DAT2[DAT2<0]<-0
DAT[DAT<0]<-0
Dat<-DAT
rownames(DAT)<-2.14*(1:nrow(DAT))

tab.dat<-NULL
for(i in 1:ncol(Dat))tab.dat[[i]]<-rep(as.numeric(rownames(DAT)),DAT[,i])


sim<-FALSE
source("/home/spl/Diffusion FP/DiffusionAnalysisFunctions12.R")
require(foreach)
require(doMC)
require(aroma.light)
require(nlme)


res2<-res<-vector("list",3)
for(r in 1:3){
    thresh<-2.14*suppressWarnings(apply(DAT[,5*(r-1)+1:5],2,function(x)min(which(x>.65*max(DAT[-nrow(DAT),5*(r-1)+1:5])))))
   thresh[1]<-Inf
   t.dat<-tab.dat[Rep==r]
   x.loc<-5*(1:floor(sort(thresh)[4]/5-1))
   x.loc<-x.loc[x.loc<120]
   set.seed(1)
   res[[r]]<-DnA.analyze(dat=t.dat,n.boot=200,record=Time[1:5]*60,Adj=2,x.loc=x.loc,
                         hi.thresh=thresh,
                         unnorm=FALSE,dat.mat=backsub.norm.trim.dat,mat.fun=mat.fun,bw.method="nrd0",plot.den=TRUE)
 }

for(r in 1:2){
  res2[[r]]<-DnA.results(res[[r]],plot.res=FALSE)
}

dev.new(width=5,height=3.5)
postscript("TotalLipid.eps", horizontal=FALSE,onefile=FALSE,paper="special",height=3.5,width=5)
par(mai=c(1,1,.1,.1))
plot(Time[1:5],sapply(bb,mean),ylab="Mean Intensity",xlab="Time (min)",cex.lab=2,cex.axis=2,cex=1.7,lwd=2)
grid()
dev.off()


   thresh<-2.14*suppressWarnings(apply(DAT,2,function(x)min(which(x>.65*max(DAT[-nrow(DAT),])))))
   thresh[1]<-Inf
   x.loc<-5*(1:floor(sort(thresh)[4]/5-1))
   x.loc<-x.loc[x.loc<121]
   thresh[thresh==Inf]<-NA
   set.seed(1)
   res.combined<-DnA.analyze(dat=tab.dat,n.boot=100,record=Time[1:5]*60,Adj=2,x.loc=x.loc,
                         hi.thresh=thresh,
                         unnorm=FALSE,dat.mat=backsub.combined.norm.trim.dat,mat.fun=mat.fun,bw.method="nrd0",plot.den=FALSE)



 dev.new(width = 8, height =6 )

#postscript("UpdatedLipidsResults.eps", horizontal=FALSE,onefile=FALSE,paper="special",height=6,width=8)
 layout(matrix(1:3, 3, 1),heights=c(4,2.5,3.7))
use<-1:nrow(DAT2)
par(bty="o",mai=c(0,.7,.1,.1),mgp=c(3,1,0))
x.convert=2.14
plot(-10,-1,xlim=c(4,2.1*80),ylim=range(0,175),xlab="",ylab="Col. Mean Intensity",main="",cex.axis=2,cex.lab=2,cex.main=1.7,xaxt="n",yaxt="n")
axis(2,at=25*(0:10),labels=rep("",11))
axis(2,at=50*(1:3),labels=50*(1:3),cex.axis=2)
grid()
legend("topleft",legend=c(paste(round(Time[1:5],2),"min"),"threshold"),fill=c(1:5,"orange"),cex=1.7,bg="white",ncol=3)#,xjust=.5,yjust=0)
abline(h=.65*max(DAT2[-nrow(DAT2),]),col="orange",lwd=3)
for(i in 1:ncol(Dat)){
 lines(x.convert*(2:nrow(DAT2)-2),DAT2[-nrow(DAT2),i],col=as.numeric(as.factor(time))[i],lty=Rep[i],lwd=4)
 y<-DAT2[,i]
# if(thresh[i]!=Inf) abline(v=thresh[i],col=as.factor(time)[i],lty=Rep[i],lwd=2)
  # tab.dat[[i]]<-rep(1:length(y),y)
}
#abline(h=38,col="orange")
#legend("topleft",legend=paste(round(Time[1:5],2),"min"),cex=2,bty="n",lwd=4,lty=1:5,ncol=3)
#legend("bottomright",legend=c(paste("Region",1:3),"threshold"),fill=c(2:4,"orange"),cex=2,bg="white")
#dev.copy2pdf(file="LipidProfiles.pdf")
#dev.off()
t1<-c(2,2,2,3,3,4)
t2<-c(3,4,5,4,5,5)
for (Func in c("A", "D")) {
    if (Func == "D") {
      par(bty="o",mai = c(.5, .7, 0, 0.1))
      plot(-1,1e6,xlim=c(4,2.1*80), ylim=
           log(c(0.2,10)),
           #quantile(res.combined$d.hat[,-(1:4)],c(.05,.95),na.rm=TRUE),
           ylab = expression("D ("*mu*"m"^2*"/s)"), yaxt="n",
           xlab = expression("Distance from Left ("*mu*"m)"),cex.axis=2,cex.lab=2)
      axis(2,at=log(rep(1:9,4)*10^rep(-2:1,each=9)),label=rep("",36))
      axis(2,at=log(c(.3,1,3,10)),label=c(.3,1,3,10),cex.axis=2)
      for(ht in log(c(.3,3,.1,1,10))) abline(h=ht,lty=3,col="lightgrey")
      for(vt in 50*(0:10)) abline(v=vt, lty=3,col="lightgrey")
      record=round(Time[1:5],2)
#      legend("topright",title="Time Pair (min)",legend=paste("(",record[c(2,2,3,3,4)],", ",record[c(3,4,4,5,5)],")",sep=""),fill=c(1:6)[-3],ncol=1,bg="white",cex=1.7)
    }
    if (Func == "A") {
      par(mai = c(0, .7, 0, 0.1))
      plot(-1,1e6,xlim=c(4,2.1*80),ylim=  c(-.05,.20),
           ylab = expression("A ("*mu*"m/s)"),
           cex.axis=2,cex.lab=2,xlab = "", xaxt = "n",yaxt="n")
      axis(2,at=.05*(-2:10),label=rep("",13))
      axis(2,at=c(-0.05,.05,0.15),label= c("-.05",".05",".15"),cex.axis=2)
      grid()
    }
    for (x in x.loc) {
      boot <- NULL
      if (Func == "D"){ 
        t.dat <- suppressWarnings(log(res.combined$d.hat[grep(paste("X", x, "X", sep = ""), 
          rownames(res.combined$d.hat)), ]))
      }
      if (Func == "A"){ 
        t.dat <- suppressWarnings(res.combined$a.hat[grep(paste("X", x, "X", sep = ""), 
          rownames(res.combined$a.hat)), ])
      }
      for(i in 5:10){
        if(!is.na(t.dat[1,i])){
#          lines(rep(x+i-7.5,2),t.dat[1,i]+2*c(1,-1)*sd(t.dat[,i]),col=t2[i-4],lwd=6,lty=3)
#          lines(x+i-7.5+c(-.6,.6),rep(t.dat[1,i],2),col=t2[i-4],lwd=6,lty=3)
          lines(rep(x+i-7.5,2),t.dat[1,i]+2*c(1,-1)*sd(t.dat[,i],na.rm=TRUE),col="gray",lwd=2)
          lines(x+i-7.5+c(-.8,.8),rep(t.dat[1,i]-2*sd(t.dat[,i],na.rm=TRUE),2),col=t1[i-4],lwd=4)
          lines(x+i-7.5+c(-.8,.8),rep(t.dat[1,i]+2*sd(t.dat[,i],na.rm=TRUE),2),col=t1[i-4],lwd=4)
          points(x+i-7.5,t.dat[1,i],col=t2[i-4],pch=18,cex=2,lwd=3)
#          lines(x+i-7.5+c(-.6,.6),rep(t.dat[1,i],2),col=t2[i-4],lwd=3)
        }
      }
    }
  }

                                        #dev.copy2pdf(file="UpdatedLipidsResults.pdf")
dev.off()

#   dat=t.dat;n.boot=500;record=Time[Rep==r];Adj=2;hi.thresh=thresh[Rep==r];dat.mat=mod.dat[Rep==r];bw.method="nrd0";low.thresh = NULL;plot.den=FALSE
  res.combined2<-DnA.results(res.combined,plot.res=FALSE)



if(0){
#dev.new()
postscript("LipidThresh.eps", horizontal=FALSE,onefile=FALSE,paper="special",height=7,width=7)
nf<-layout(matrix(1:2,2,1),heights=c(1,.4))
   
par(mai=c(0,1,.2,.1),mgp=c(3,1,0))
x.convert=2.14
plot(-10,-1,xlim=x.convert*c(59,78),ylim=c(9000,13000),xlab="",ylab="Row Intensity",main="",cex.axis=1.7,cex.lab=1.7,cex.main=1.7,xaxt="n")
grid()
for(i in 1:ncol(mn.Dat)){
 lines(x.convert*(2:nrow(Dat)-2),Dat[-nrow(Dat),i],col=as.factor(time)[i],lty=Rep[i],lwd=4)
}
legend("topleft",legend=paste(round(Time[1:5],2),"min"),cex=1.6,bty="n",fill=1:5,ncol=3)
legend(126,178,legend=paste("Reg",1:3),lty=1:3,cex=1.6,lwd=4,bty="n")

move<-c(0,0,-.75,-.25,-.25,0,0,.25,.75,.25,0,0,0,0,0)
par(mai=c(1,1,0,.1),mgp=c(3,1,0))
plot(-10,-1,xlim=x.convert*c(61,80),ylim=c(180,185),xlab=expression("Distance from Bottom ("*mu*"m)"),ylab="",main="",cex.axis=1.7,cex.lab=1.7,cex.main=1.7,yaxt="n")
for(i in 1:ncol(mn.Dat)){
 if(thresh[i]!=Inf) lines(rep(thresh[i]*x.convert+move[i],2),c(170,190),col=as.factor(time)[i],lty=Rep[i],lwd=4)
}

dev.off()
}

thresh[thresh==Inf]<-NA
sim<-FALSE
source("/home/spl/Diffusion FP/DiffusionAnalysisFunctions12.R")
require(foreach)
require(doMC)
require(aroma.light)
require(nlme)

res2<-res<-vector("list",3)
for(r in 1:3){
   t.dat<-tab.dat[Rep==r]
#   x.loc<-seq(quantile(unlist(t.dat),.02),min(thresh[Rep==r],na.rm=TRUE),length.out=49)
     # x.loc<-seq(,min(thresh[Rep==r],na.rm=TRUE),by=1)-.1
   x.loc<-2*(1:floor(sort(thresh[Rep==r]-2)[4]/2))
  
   set.seed(1)
   res[[r]]<-DnA.analyze(dat=t.dat,n.boot=100,record=Time[Rep==r],Adj=2,x.loc=x.loc,
                         hi.thresh=thresh[Rep==r],
                         unnorm=FALSE,dat.mat=mod.dat[Rep==r],bw.method="nrd0")
 }
#   dat=t.dat;n.boot=500;record=Time[Rep==r];Adj=2;hi.thresh=thresh[Rep==r];dat.mat=mod.dat[Rep==r];bw.method="nrd0";low.thresh = NULL;plot.den=FALSE
for(r in 1:2){
  res2[[r]]<-DnA.results(res[[r]],plot.res=FALSE)
}

exp(c(2.95,2.88,3.19))*diff.unit.convert

save(res,file="LipidsPreresultsUnnorm.Rdata")
save(res2,file="LipidsResultsUnnorm.Rdata")

load(file="LipidsPreresultsUnnormBG38.Rdata")
names(res[[1]])
thresh
res[[1]]$d.hat[1:sum(res[[1]]$x.loc<58),-c(1:4,6:7)]

if(0){
COL<-makeTransparent(1:8)
x11(width=12,height=9)
layout(matrix(1:2,2,1))
par(mai=c(0,1,1,.5))
unit.convert<-1/(.468^2*60)
x.convert=.468
plot(-1e6,-1e6,xlim=x.convert*range(res2[[1]]$x.loc),ylim=unit.convert*c(1,500),log="y",type="l",ylab=expression("Diffusion ("*mu*"m"^2*"/sec)"),xlab="",xaxt="n", cex.lab=1.7,cex.axis=1.7,cex.main=1.7,main="Lipid Data Results")
grid(equilogs=FALSE)
axis(2,at=rep(1:9,3)*10^rep(-2:0,each=9),label=rep("",27))
for(i in 1:3){
  lines(x.convert*res2[[i]]$x.loc,unit.convert*exp(res2[[i]]$log.d),col=COL[i+1],lwd=5)
  lines(x.convert*res2[[i]]$x.loc,unit.convert*exp(res2[[i]]$log.d+2*res2[[i]]$se.log.d),col=COL[i+1],lwd=3,lty=2)
  lines(x.convert*res2[[i]]$x.loc,unit.convert*exp(res2[[i]]$log.d-2*res2[[i]]$se.log.d),col=COL[i+1], lwd=3, lty=2)
}
legend("bottomright",legend=paste("Rep",1:3),col=COL[2:4],lwd=2,cex=1.7)
par(mai=c(1,1,0,.5))
unit.convert<-1/(.468*60)
plot(-1e6,-1e6,xlim=x.convert*range(res2[[1]]$x.loc),ylim=c(0,.28), #unit.convert*range(res2[[1]]$a,na.rm=TRUE),
     type="l",ylab=expression("Drift ("*mu*"m/sec)"),xlab=expression("Distance from Top ("*mu*"m)"),cex.lab=1.7,cex.axis=1.7)
grid(equilogs=FALSE)
for(i in 1:3){
  lines(x.convert*res2[[i]]$x.loc,unit.convert*res2[[i]]$a,col=COL[i+1],lwd=5)
  
  lines(x.convert*res2[[i]]$x.loc,unit.convert*(res2[[i]]$a+2*res2[[i]]$se.a),col=COL[i+1],lwd=3,lty=2)
  lines(x.convert*res2[[i]]$x.loc,unit.convert*(res2[[i]]$a-2*res2[[i]]$se.a),col=COL[i+1], lwd=3, lty=2)
}
}

postscript("LipidResults.eps", horizontal=FALSE,onefile=FALSE,paper="special",height=9,width=12)
#dev.new(height=9,width=12)
COL<-makeTransparent(1:8)
COL<-1:8
layout(matrix(1:2,2,1),heights=c(3.5,4.8))
par(mai=c(0,1.5,.2,.5),mgp=c(4,1.65,0))
diff.unit.convert<-2.14^2/60
x.convert=2.14
plot(-1e6,-1e6,xlim=x.convert*range(res2[[1]]$x.loc),ylim=c(.05,4),log="y",type="l",ylab=expression("Diffusion ("*mu*"m"^2*"/s)"),xlab="",xaxt="n", cex.lab=2.7,cex.axis=2.7,cex.main=1.7,main="",yaxt="n")
axis(2,at=rep(1:9,3)*10^rep(-2:0,each=9),label=rep("",27))
axis(2,at=c(.1,.3,1,3),label=c("0.1","0.3","1","3"),cex.axis=2.7)
axis(2,at=c(.3),label=c("0.3"),cex.axis=2.7)
grid(equilogs=FALSE)
for(i in 1:3){
  points(x.convert*res2[[i]]$x.loc,diff.unit.convert*exp(res2[[i]]$log.d),col=COL[i+1],lwd=5)
for(j in 1:length(res2[[i]]$x.loc))  lines(rep(x.convert*res2[[i]]$x.loc[j],2),diff.unit.convert*exp(res2[[i]]$log.d[j]+c(-2,2)*res2[[i]]$se.log.d[j]),col=COL[i+1],lwd=3)
   lines(x.convert*res2[[i]]$smooth.log.d$x,diff.unit.convert*exp(res2[[i]]$smooth.log.d$y),col=COL[i+1],lwd=5)

}
legend("bottomright",legend=paste("Region",1:3),fill=COL[2:4],cex=2.4,bg="white")
par(mai=c(1.5,1.5,0,.5),mgp=c(4,1.65,0))
drift.unit.convert<-2.14/60
plot(-1e6,-1e6,xlim=x.convert*range(res2[[1]]$x.loc),ylim=c(0,.384), #unit.convert*range(res2[[1]]$a,na.rm=TRUE),
     type="l",ylab=expression("Drift ("*mu*"m/s)"),xlab=expression("Distance from Bottom ("*mu*"m)"),cex.lab=2.7,cex.axis=3,yaxt="n")
axis(2,at=(0:4)/10,label=rep("",5))
axis(2,at=c(0,.3),label=c("0","0.3"),cex.axis=2.7)
axis(2,at=c(.1),label=c("0.1"),cex.axis=2.7)
axis(2,at=c(.2),label=c("0.2"),cex.axis=2.7)
grid(equilogs=FALSE)
for(i in 1:3){
  points(x.convert*res2[[i]]$x.loc,drift.unit.convert*res2[[i]]$a,col=COL[i+1],lwd=5)
  for(j in 1:length(res2[[i]]$x.loc))lines(rep(x.convert*res2[[i]]$x.loc[j],2),drift.unit.convert*(res2[[i]]$a[j]+c(-2,2)*res2[[i]]$se.a[j]),col=COL[i+1],lwd=3)
   lines(x.convert*res2[[i]]$smooth.a$x,drift.unit.convert*res2[[i]]$smooth.a$y,col=COL[i+1],lwd=5)
}


#dev.copy2pdf(file="LipidResults.pdf")
dev.off()







plot(1:25,pch=1:25)
par(mfrow=c(4,4),mai=rep(.2,4))
ylim<-range(res2[[r]]$a.pair.est,na.rm=TRUE)
ylim=c(-1,3)
for(i in 1:ncol(res2[[r]]$a.pair.est)){
  plot(res2[[r]]$x.loc,res2[[r]]$a.pair.est[,i],ylim=ylim)
  grid()
  for(j in 1:nrow(res2[[r]]$a.pair.est))lines(rep(res2[[r]]$x.loc[j],2),res2[[r]]$a.pair.est[j,i]+c(-1,1)*res2[[r]]$boot.sd.a.pair.est[j,i])
}

den<-res[[r]]$keep.den
plot(den[,1],den[,2],ylim=range(den[,-1],na.rm=TRUE),type="l")
for(i in 3:6) lines(den[,1],den[,i],col=i-1)

par(mfrow=c(4,4),mai=rep(.2,4))
ylim<-range(res2[[r]]$d.pair.est,na.rm=TRUE)
#ylim=c(-1,3)
for(i in 1:ncol(res2[[r]]$d.pair.est)){
  plot(res2[[r]]$x.loc,res2[[r]]$d.pair.est[,i],ylim=ylim)
  grid()
  for(j in 1:nrow(res2[[r]]$d.pair.est))lines(rep(res2[[r]]$x.loc[j],2),res2[[r]]$d.pair.est[j,i]+c(-1,1)*res2[[r]]$boot.sd.d.pair.est[j,i])
}

den<-res[[r]]$keep.den
plot(den[,1],den[,2],ylim=range(den[,-1],na.rm=TRUE),type="l")
for(i in 3:6) lines(den[,1],den[,i],col=i-1)



sapply(res2,function(x)exp(x$log.constant.D+x$se.log.constant.D*c(-2,2))*2.14^2/60)

FIG<-1; while(!is.null(FIG)) FIG<-dev.off()

x<-res[[r]]$rel.ent
y<-res[[r]]$hat[,grep("A",colnames(res[[r]]$hat))]

t.y<-y[rownames(y)==rownames(y)[i],]
t.sd<-apply(t.y,2,sd)
t.mn<-colMeans(t.y)
YLIM=range(t.mn+.5*c(-1,1)*t.sd)
plot(x,t.mn,ylim=YLIM)
grid(equilogs=FALSE)
for(j in 1:length(x)) lines(rep(x[j],2),t.mn[j]+.5*c(-1,1)*t.sd[j])COL<-makeTransparent(1:8)

X.lim<-range(sapply(res,function(x)x$keep.den[,1]))
plot()


if(0){
plot(-1e6,-1e6,xlim=range(res2[[1]]$x.loc),ylim=range(res2[[1]]$u),type="l",ylab="Potential Landscape",xlab="X",
     cex.lab=1.7,cex.axis=1.7)
grid(equilogs=FALSE)
for(i in 1:3){
  lines(res2[[i]]$x.loc,res2[[i]]$u,col=COL[i+1],lwd=3)
  lines(res2[[i]]$x.loc,res2[[i]]$u+2*res2[[i]]$se.u,col=COL[i+1],lwd=2,lty=2)
  lines(res2[[i]]$x.loc,res2[[i]]$u-2*res2[[i]]$se.u,col=COL[i+1], lwd=2, lty=2)
}
}

