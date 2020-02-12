library(reshape)
library(plyr)
library(dplyr)
library(ggplot2)

library(doParallel)
library(foreach)

cl= makeCluster(5)
registerDoParallel(cl)

dirMy="/home/laurence-kell/Desktop/papers/data-poor"

## FAO Taxa
load(file.path(dirMy,"data/fao_taxa.rda"))

## ICES
load(file.path(dirMy,"data/sag_2019.RData"))

## Myers database
load(file.path(dirMy,"data/myersDB.RData"))
myers=subset(timeseries,tsid%in%c("BdivBmsytouse-dimensionless","TCbest-MT"))
myers=cast(myers,assessid+stockid+tsyear~tsid,value="tsvalue")
names(myers)[3:5]=c("year","biomass","catch")
myers=subset(myers,!(is.na(biomass)|is.na(catch)))
myers=merge(myers,stock[,c("stockid","scientificname","region")],by="stockid")
names(myers)[6]="species"

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

fnBD<-function(assessid,stockid,biomass,catch,year,spp,dir,cv=NULL){
  
  stockid =stockid[1]
  assessid=assessid[1]
  
  print(paste(spp[1],assessid))
  
  cat(paste(assessid,"\n"),file=file.path(dir,stockid),append=TRUE)

  if (is.numeric(cv)){
    set.seed(1233)
    biomass=biomass*exp(rlnorm(length(biomass),0,cv))}
  biomass=biomass*mean(catch)*4
  
  cdp<-try(sraplus::format_driors(
    b_ref_type ="b",
    index      = biomass,
    index_year = year,
    catch      = catch,
    years      = year,
    initial_state=biomass[1],
    terminal_state=biomass[length(biomass)],
    initial_state_cv=0.05,
    terminal_state_cv=0.05,
    taxa       = spp[1],
    use_heuristics=TRUE))
  
  if ("try-error"%in%is(cdp)) return(NULL)
  
  fit<-try(sraplus::fit_sraplus(
    driors = cdp,      
    engine = "tmb",
    model  = "sraplus_tmb",
    adapt_delta = 0.9,
    max_treedepth = 10,
    n_keep = 4000,
    chains = 1,
    cores = 1,
    estimate_qslope = FALSE,
    estimate_proc_error = TRUE))
  
  if (!("try-error"%in%is(fit)))
    save(fit,file=file.path(dir,paste(stockid,assessid,"RData",sep=".")))
  
  save(biomass,file=file.path(dir,paste(stockid,assessid,"u",sep=".")))
  
  "try-error"%in%is(fit)}

control=myers[!duplicated(myers[,c("assessid","stockid")]),c("assessid","stockid")]

##Benchmark
foreach(i=seq(dim(control)[1]),
        .combine="c",
        .packages=c("sraplus"),
        .export  =c("myers","control")
        ) %do% 
     with(subset(myers,assessid==control[i,"assessid"]&stockid==control[i,"stockid"]),  
          fnBD(assessid,stockid,biomass,catch,year,species,dir=file.path(dirMy,"results/bd")))

##CV of 40%
foreach(i=rev(seq(dim(control)[1])),
        .combine="c",
        .packages=c("sraplus"),
        .export  =c("myers","control")
) %do% 
  with(subset(myers,assessid==control[i,"assessid"]&stockid==control[i,"stockid"]),  
       fnBD(assessid,stockid,biomass,catch,year,species,dir=file.path(dirMy,"results/bd4"),cv=0.4))

##SRA
fnSRA<-function(assessid,stockid,biomass,catch,year,spp,dir,cv=NULL){
  
  stockid =stockid[1]
  assessid=assessid[1]
  
  print(paste(spp[1],assessid))
  
  cat(paste(assessid,"\n"),file=file.path(dir,stockid),append=TRUE)
  
  if (is.numeric(cv)){
    set.seed(1233)
    biomass=biomass*exp(rlnorm(length(biomass),0,cv))}
  
  cdp<-try(sraplus::format_driors(
    taxa           =spp[1],
    use_heuristics =TRUE,
    catch          =catch,
    years          =year))
  
  if ("try-error"%in%is(cdp)) return(NULL)
  
  fit<-try(sraplus::fit_sraplus(
    driors = cdp,
    engine = "sir",
    adapt_delta = 0.9,
    max_treedepth = 10,
    n_keep = 4000,
    chains = 1,
    cores = 1,
    estimate_qslope = FALSE,
    estimate_proc_error = TRUE))
  
  if (!("try-error"%in%is(fit)))
    save(fit,file=file.path(dir,paste(stockid,assessid,"RData",sep=".")))
  
  "try-error"%in%is(fit)}

foreach(i=rev(seq(dim(control)[1])),
        .combine="c",
        .packages=c("sraplus"),
        .export  =c("myers","control")
) %do% 
  with(subset(myers,assessid==control[i,"assessid"]&stockid==control[i,"stockid"]),  
       fnSRA(assessid,stockid,biomass,catch,year,species,dir=file.path(dirMy,"results/sra")))

runs=as.data.frame(rbind(
           cbind(run="bd",   assessid=list.files(file.path(dirMy,"results/bd"))),
           cbind(run="bd4",  assessid=list.files(file.path(dirMy,"results/bd4"))),
           cbind(run="sra",  assessid=list.files(file.path(dirMy,"results/sra")))))
runs=runs[grep(".RData",runs$assessid),]

runs=mdply(runs, function(run,assessid) {
       print(paste(run,assessid))  
       load(file.path(file.path(dirMy,"results"),run,assessid))
       
       if ("try-error"%in%is(fit)) return (NULL)
       
       fit$results})
runs$assessid=gsub(".RData","",runs$assessid)
runs$assessid=laply(strsplit(runs[,2],"\\."),function(x) x[[2]])

save(runs,file=file.path(dirMy,"results/myersRuns.RData"))



##################################
    driors,
    ~ sfs(
      driors = .x,
      engine = "tmb",
      model = "sraplus_tmb",
      adapt_delta = 0.9,
      max_treedepth = 10,
      n_keep = 4000,
      chains = 1,
      cores = 1,
      estimate_qslope = FALSE,
      estimate_proc_error = TRUE)

driors = map2(
  taxa,
  data,
  ~
    format_driors(
      taxa = .x,shape_prior=1.01,      #use_heuristics = T,shape_prior=2,
      catch = .y$capture,
      years = .y$year,
      initial_state = 0.9,initial_state_cv = 0.2,b_ref_type = "k",
      index = .y$E1[!is.na(.y$E1)],index_years=.y$year[!is.na(.y$E1)]) 