library(reshape)
library(plyr)
library(dplyr)
library(ggplot2)

library(sraplus)

load("/home/laurence-kell/Desktop/rfmo/iccat/sc-eco/iccat/jabba/jabbaiccat.RData")

jabbaiccat=jabbaiccat[sort(names(jabbaiccat))]

spp=c("Thunnus obesus",
      "Kajikia albida",
      "Xiphias gladius",
      "Makaira nigricans",
      "Thunnus albacares")
names(spp)=names(jabbaiccat)
jb=ldply(jabbaiccat,function(x)
       data.frame(year=x$timeseries$year,b=x$timeseries$ssb,b_div_bmsy=x$timeseries$ssb/x$refpts$bmsy[1]))
names(jb)[1]="sp"


bds=mlply(data.frame(sp=names(spp)), function(sp){

  x=jabbaiccat[[sp]]
  
  cdp<-try(sraplus::format_driors(
    use_heuristics=!TRUE,
    taxa          =spp[sp],
    
    index      = x$timeseries$ssb,
    index_year = x$timeseries$year,
    catch      = x$timeseries$catch,
    years      = x$timeseries$year
     
    # k_prior	            =x$pfunc$k,
    # growth_rate_prior   =x$pfunc$r,
    # shape_prior         =x$pfunc$p+1,
    # initial_state       =x$timeseries$ssb[1]/x$refpts[1,"k"],
    # terminal_state      =x$timeseries$ssb[length(x$timeseries$ssb)]/x$refpts[1,"k"],

    # k_prior_cv          =0.1,
    # growth_rate_prior_cv=0.1,
    # shape_prior_cv      =0.1,
    # initial_state_cv    =0.1,
    # terminal_state_cv   =0.1
    ))
  
  fit<-try(sraplus::fit_sraplus(
    driors = cdp,      
    engine = "tmb",
    model  = "sraplus_tmb",
    estimate_qslope = FALSE,
    estimate_proc_error = TRUE))

  if ("try-error"%in%is(fit)) return(NULL)
  
  fit})
bd=ldply(bds,function(x){
  if (is.null(x)) return(NULL)
  subset(x$results,variable%in%c("b_div_bmsy","b"))})


#### SRA
sras=mlply(data.frame(sp=names(spp)), function(sp){
  
  x=jabbaiccat[[sp]]
  
  cdp<-try(sraplus::format_driors(
    use_heuristics=TRUE,
    taxa          =spp[sp],
    
    catch      = x$timeseries$catch,
    years      = x$timeseries$year
  ))
  
  
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
  
  if ("try-error"%in%is(fit)) return(NULL)
  
  fit})
sra=ldply(sras,function(x){
  if (is.null(x)) return(NULL)
  subset(x$results,variable%in%c("b_div_bmsy","b"))})

ggplot(rbind(cbind(Method="Catch & Index",subset(bd, variable=="b")),
             cbind(Method="Catch",        subset(sra,variable=="b"))))+
  geom_ribbon(aes(year,ymax=upper,ymin=lower),fill="red",alpha=0.25)+
  geom_line(aes(year,mean),col="red")+
  geom_line(aes(year,b),data=jb)+
  facet_grid(sp~Method,scale="free")

ggplot(rbind(cbind(Method="Catch & Index",subset(bd, variable=="b_div_bmsy")),
             cbind(Method="Catch",        subset(sra,variable=="b_div_bmsy"))))+
  geom_ribbon(aes(year,ymax=upper,ymin=lower),fill="red",alpha=0.25)+
  geom_line(aes(year,mean),col="red")+
  geom_line(aes(year,b_div_bmsy),data=jb)+
  facet_grid(sp~Method,scale="free")
