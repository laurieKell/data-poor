library(reshape)
library(plyr)
library(dplyr)
library(ggplot2)
library(FLife)
library(mydas)
library(popbio)

dirMy="/home/laurence-kell/Desktop/papers/data-poor"

load(file.path(dirMy,"data/myers.RData"))

## Fishbase
load(file.path(dirMy,"data/fb.RData"))
fb=ddply(fb, .(species),  function(x)  data.frame(
  k   =mean(x$k,   na.rm=T),
  linf=mean(x$linf,na.rm=T),
  t0  =mean(x$t0,  na.rm=T),
  l50 =mean(x$lmat,   na.rm=T),
  a   =mean(x$a,      na.rm=T),
  b   =mean(x$b,      na.rm=T),
  m   =mean(x$m,      na.rm=T)))

source('~/Desktop/flr/FLife/R/FLife-lhPar.R', echo=TRUE)

popdyn2<-function(x) {
  n<<-n+1
  print(n)
  if (is.na(x["linf"])) return(NULL) 
  
  sx=x[-1]
  x=try(as(x[,names(x)[!is.na(x)]],"FLPar"))
  if ("try.error"%in%is(x)) return(NULL)
  x=lhPar(x)
  
  x=try(popdyn(x))
  if ("try.error"%in%is(x)) return(NULL) else return(model.frame(x))}

n=0
pd=subset(fb,species%in%unique(myers$species))
pd=adply(pd[,c("species","linf","k","t0","a","b","l50")],1,popdyn2)
save(pd,file=file.path(dirMy,"data/pd.RData"))