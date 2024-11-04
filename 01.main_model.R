library(purrr)
library(tidyverse)
library(gnm)
library(tidyverse)
library(dlnm)
library(pbs)
library(splines)
library(mixmeta)
library(FluMoDL)
library(data.table)
rm(list=ls())
setwd("")
df <- read_rds('analysis.rds')

df %>% 
  dplyr::mutate(dow=weekdays(date),
                year=year(date),
                month=month(date),
                doy=yday(date)) %>% 
  dplyr::arrange(code) %>% 
  dplyr::mutate(stratum=factor(paste0(code,'_',year,'_',month)))-> df

ylist <- c('count.all','count.female','count.male',
           'count.age0_64','count.age65_',
           'count.cardio','count.resp',
           'count.mental','count.nervous',
           'count.injury','count.symptom')

outcome=c('All cause',
          'Female','Male',
          'Younger(Age<65)','Older(Age≥65)',
          'Cardiovascular','Respiratory',
          'Mental','Nervous',
          'Injury','Symptom')

code <- unique(df$code)
dlist <- lapply(code, function(x) df[df$code==x,])

# predict temperature percentiles
predper <- c(seq(0,1,0.1),2:98,seq(99,100,0.1))
set.seed(13041975)
tmeancountry <- round(rowMeans(sapply(dlist,function(x) quantile(jitter(x$tmean),
                                                                 predper/100,na.rm=T))),2)
# basis function
lag=14
lagnk=2
argvar <- list(fun='ns',df=3)
cb <- crossbasis(df$tmean,lag=lag,argvar=argvar,
                 arglag=list(fun='ns',knots=logknots(lag,lagnk)),
                 group=df$code)


reslist <- list()
RRlist <- list()

# time function
spldoy <- onebasis(df$doy,'ns',df=6)
# formula
fmod <- y ~ cb + ns(rh_mv02,3) + holiday+ factor(dow) + spldoy:factor(year)

for (i in 1:length(ylist)) {
  
  # i=1
  print(i)
  df$y = df[[ylist[i]]]
  
  df2 <- df
  df2 %>% 
    group_by(stratum) %>%
    summarise(all=sum(y)) %>% 
    mutate(keep=ifelse(all>0,T,F)) %>% 
    right_join(df2,by='stratum') -> df2
  
  mod <- gnm(fmod,eliminate=stratum,
             na.action="na.exclude",subset=keep, # keep only strata with at least one death, otherwise the estimate would have bias
             data=df2,family=quasipoisson)
  
  pred <- crosspred(cb,mod,
                    at=tmeancountry,
                    cen=15)
  
  mmt <- pred$predvar[11:109][which.min(pred$allRRfit[11:109])] # limit the mmt in 1%-99% range

  pred <- crosspred(cb,mod,
                    at=tmeancountry,
                    cen=mmt)
  
  ##############################################
  # interaction
  intval <- c(0,1)
  cbint1 <- cb * (df$po - intval[1])
  cbint2 <- cb * (df$po - intval[2])
  
  # RUN THE MODELS
  modint1 <- update(mod, .~. + cbint1 + I(po-intval[1]))
  modint2 <- update(mod, .~. + cbint2 + I(po-intval[2]))
  
  
  
  # PREDICT FOR EACH OF THE TWO IMD VALUES
  ind <- grep('cb',names(coef(mod)))
  cpint1 <- crosspred(cb,
                      coef=coef(modint1)[ind],vcov=vcov(modint1)[ind,ind],
                      cen=mmt, model.link='log',
                      at=tmeancountry,
                      cumul=T)
  
  cpint2 <- crosspred(cb,
                      coef=coef(modint2)[ind],vcov=vcov(modint2)[ind,ind],
                      cen=mmt, model.link='log',
                      at=tmeancountry,
                      cumul=T)
  
  cpint <- crosspred(cbint2,
                     modint2,
                     cen=mmt,
                     at=tmeancountry,
                     cumul=T)

  
  data.frame(outcome=outcome[i],
             mmt=mmt,
             perc=names(tmeancountry),
             predvar=tmeancountry,
             nonintRR=pred$allRRfit,
             nonintRRlow=pred$allRRlow,
             nonintRRhigh=pred$allRRhigh,
             intRR_0=cpint1$allRRfit,
             intRRlow_0=cpint1$allRRlow,
             intRRhigh_0=cpint1$allRRhigh,
             intRR_1=cpint2$allRRfit,
             intRRlow_1=cpint2$allRRlow,
             intRRhigh_1=cpint2$allRRhigh,
             intRR=cpint$allRRfit,
             intRRlow=cpint$allRRlow,
             intRRhigh=cpint$allRRhigh
  ) -> reslist[[i]]
  
  logrr <- c(
    cpint1$allfit[pred$predvar==tmeancountry['90.0%']],
    cpint2$allfit[pred$predvar==tmeancountry['90.0%']]
  )
  logvar <- c(
    (cpint1$allse[pred$predvar==tmeancountry['90.0%']])^2,
    (cpint2$allse[pred$predvar==tmeancountry['90.0%']])^2
  )
  p90=summary(mixmeta(logrr,logvar))$qstat$pvalue
  
  rbind(
    data.frame(outcome=outcome[i],
               mmt=mmt,
               perc='90.0%',
               type='Non-outage',
               RR=cpint1$allRRfit[pred$predvar==tmeancountry['90.0%']],
               RRlow=cpint1$allRRlow[pred$predvar==tmeancountry['90.0%']],
               RRhigh=cpint1$allRRhigh[pred$predvar==tmeancountry['90.0%']]),
    data.frame(outcome=outcome[i],
               mmt=mmt,
               perc='90.0%',
               type='Outage',
               RR=cpint2$allRRfit[pred$predvar==tmeancountry['90.0%']],
               RRlow=cpint2$allRRlow[pred$predvar==tmeancountry['90.0%']],
               RRhigh=cpint2$allRRhigh[pred$predvar==tmeancountry['90.0%']])
  ) -> RR90
  RR90$p.inter=p90
  
  logrr <- c(
    cpint1$allfit[pred$predvar==tmeancountry['95.0%']],
    cpint2$allfit[pred$predvar==tmeancountry['95.0%']]
  )
  logvar <- c(
    (cpint1$allse[pred$predvar==tmeancountry['95.0%']])^2,
    (cpint2$allse[pred$predvar==tmeancountry['95.0%']])^2
  )
  p95=summary(mixmeta(logrr,logvar))$qstat$pvalue
  
  rbind(
    data.frame(outcome=outcome[i],
               mmt=mmt,
               perc='95.0%',
               type='Non-outage',
               RR=cpint1$allRRfit[pred$predvar==tmeancountry['95.0%']],
               RRlow=cpint1$allRRlow[pred$predvar==tmeancountry['95.0%']],
               RRhigh=cpint1$allRRhigh[pred$predvar==tmeancountry['95.0%']]),
    data.frame(outcome=outcome[i],
               mmt=mmt,
               perc='95.0%',
               type='Outage',
               RR=cpint2$allRRfit[pred$predvar==tmeancountry['95.0%']],
               RRlow=cpint2$allRRlow[pred$predvar==tmeancountry['95.0%']],
               RRhigh=cpint2$allRRhigh[pred$predvar==tmeancountry['95.0%']])
  ) -> RR95
  RR95$p.inter=p95
  
  logrr <- c(
    cpint1$allfit[pred$predvar==tmeancountry['99.0%']],
    cpint2$allfit[pred$predvar==tmeancountry['99.0%']]
  )
  logvar <- c(
    (cpint1$allse[pred$predvar==tmeancountry['99.0%']])^2,
    (cpint2$allse[pred$predvar==tmeancountry['99.0%']])^2
  )
  p99=summary(mixmeta(logrr,logvar))$qstat$pvalue
  
  rbind(
    data.frame(outcome=outcome[i],
               mmt=mmt,
               perc='99.0%',
               type='Non-outage',
               RR=cpint1$allRRfit[pred$predvar==tmeancountry['99.0%']],
               RRlow=cpint1$allRRlow[pred$predvar==tmeancountry['99.0%']],
               RRhigh=cpint1$allRRhigh[pred$predvar==tmeancountry['99.0%']]),
    data.frame(outcome=outcome[i],
               mmt=mmt,
               perc='99.0%',
               type='Outage',
               RR=cpint2$allRRfit[pred$predvar==tmeancountry['99.0%']],
               RRlow=cpint2$allRRlow[pred$predvar==tmeancountry['99.0%']],
               RRhigh=cpint2$allRRhigh[pred$predvar==tmeancountry['99.0%']])
  ) -> RR99
  
  RR99$p.inter=p99
  
  logrr <- c(
    cpint1$allfit[pred$predvar==tmeancountry['10.0%']],
    cpint2$allfit[pred$predvar==tmeancountry['10.0%']]
  )
  logvar <- c(
    (cpint1$allse[pred$predvar==tmeancountry['10.0%']])^2,
    (cpint2$allse[pred$predvar==tmeancountry['10.0%']])^2
  )
  p10=summary(mixmeta(logrr,logvar))$qstat$pvalue
  
  rbind(
    data.frame(outcome=outcome[i],
               mmt=mmt,
               perc='10.0%',
               type='Non-outage',
               RR=cpint1$allRRfit[pred$predvar==tmeancountry['10.0%']],
               RRlow=cpint1$allRRlow[pred$predvar==tmeancountry['10.0%']],
               RRhigh=cpint1$allRRhigh[pred$predvar==tmeancountry['10.0%']]),
    data.frame(outcome=outcome[i],
               mmt=mmt,
               perc='10.0%',
               type='Outage',
               RR=cpint2$allRRfit[pred$predvar==tmeancountry['10.0%']],
               RRlow=cpint2$allRRlow[pred$predvar==tmeancountry['10.0%']],
               RRhigh=cpint2$allRRhigh[pred$predvar==tmeancountry['10.0%']])
  ) -> RR10
  RR10$p.inter=p10
  
  logrr <- c(
    cpint1$allfit[pred$predvar==tmeancountry['5.0%']],
    cpint2$allfit[pred$predvar==tmeancountry['5.0%']]
  )
  logvar <- c(
    (cpint1$allse[pred$predvar==tmeancountry['5.0%']])^2,
    (cpint2$allse[pred$predvar==tmeancountry['5.0%']])^2
  )
  p5=summary(mixmeta(logrr,logvar))$qstat$pvalue
  
  rbind(
    data.frame(outcome=outcome[i],
               mmt=mmt,
               perc='5.0%',
               type='Non-outage',
               RR=cpint1$allRRfit[pred$predvar==tmeancountry['5.0%']],
               RRlow=cpint1$allRRlow[pred$predvar==tmeancountry['5.0%']],
               RRhigh=cpint1$allRRhigh[pred$predvar==tmeancountry['5.0%']]),
    data.frame(outcome=outcome[i],
               mmt=mmt,
               perc='5.0%',
               type='Outage',
               RR=cpint2$allRRfit[pred$predvar==tmeancountry['5.0%']],
               RRlow=cpint2$allRRlow[pred$predvar==tmeancountry['5.0%']],
               RRhigh=cpint2$allRRhigh[pred$predvar==tmeancountry['5.0%']])
  ) -> RR5
  RR5$p.inter=p5
  
  logrr <- c(
    cpint1$allfit[pred$predvar==tmeancountry['1.0%']],
    cpint2$allfit[pred$predvar==tmeancountry['1.0%']]
  )
  logvar <- c(
    (cpint1$allse[pred$predvar==tmeancountry['1.0%']])^2,
    (cpint2$allse[pred$predvar==tmeancountry['1.0%']])^2
  )
  p1=summary(mixmeta(logrr,logvar))$qstat$pvalue
  
  rbind(
    data.frame(outcome=outcome[i],
               mmt=mmt,
               perc='1.0%',
               type='Non-outage',
               RR=cpint1$allRRfit[pred$predvar==tmeancountry['1.0%']],
               RRlow=cpint1$allRRlow[pred$predvar==tmeancountry['1.0%']],
               RRhigh=cpint1$allRRhigh[pred$predvar==tmeancountry['1.0%']]),
    data.frame(outcome=outcome[i],
               mmt=mmt,
               perc='1.0%',
               type='Outage',
               RR=cpint2$allRRfit[pred$predvar==tmeancountry['1.0%']],
               RRlow=cpint2$allRRlow[pred$predvar==tmeancountry['1.0%']],
               RRhigh=cpint2$allRRhigh[pred$predvar==tmeancountry['1.0%']])
  ) -> RR1
  
  RR1$p.inter=p1
  
  rbind(RR90,RR95,RR99,RR10,RR5,RR1) -> RRlist[[i]]
  
}


## result display ---------------------------------------------------------------
### Figure 2 ER plot ----------------------------------------------------------------------
source('plot_theme.R')
plotf <- function(outcome){
  # outcome='All cause'
  rbindlist(reslist) %>%
    dplyr::filter(outcome==!!outcome) -> tmp
  ggplot(tmp) +
    geom_ribbon(aes(x=predvar,ymin=nonintRRlow,ymax=nonintRRhigh),fill='#E8E9EB')+
    geom_line(aes(x=predvar,y=nonintRR),
              color=ifelse(tmp$predvar>tmp$mmt,2,4),linewidth=1) +
    scale_x_continuous(n.breaks=10) +
    geom_hline(yintercept=1,linetype='dashed')+
    labs(title=outcome,x=expression("Temperature ("~degree~"C)"),y='Relative risks(95% CI)') +
    theme_Publication(base_size=16) 
  
}
plotfx <- function(outcome){
  # outcome='All cause'
  rbindlist(reslist) %>%
    dplyr::filter(outcome==!!outcome) -> tmp
  ggplot(tmp) +
    geom_ribbon(aes(x=predvar,ymin=nonintRRlow,ymax=nonintRRhigh),fill='#E8E9EB')+
    geom_line(aes(x=predvar,y=nonintRR),
              color=ifelse(tmp$predvar>tmp$mmt,2,4),linewidth=1) +
    scale_x_continuous(n.breaks=10) +
    geom_hline(yintercept=1,linetype='dashed')+
    labs(title=outcome,x=expression("Temperature ("~degree~"C)"),y='Relative risks(95% CI)') +
    theme_Publication_x(base_size=16) 
  
}
plotfy <- function(outcome){
  # outcome='All cause'
  rbindlist(reslist) %>%
    dplyr::filter(outcome==!!outcome) -> tmp
  ggplot(tmp) +
    geom_ribbon(aes(x=predvar,ymin=nonintRRlow,ymax=nonintRRhigh),fill='#E8E9EB')+
    geom_line(aes(x=predvar,y=nonintRR),
              color=ifelse(tmp$predvar>tmp$mmt,2,4),linewidth=1) +
    scale_x_continuous(n.breaks=10) +
    geom_hline(yintercept=1,linetype='dashed')+
    labs(title=outcome,x=expression("Temperature ("~degree~"C)"),y='Relative risks(95% CI)') +
    theme_Publication_y(base_size=16) 
  
}


p1 <- plotf('All cause')
p2 <- plotf('Male')
p3 <- plotf('Female')
p4 <- plotf('Younger(Age<65)')
p5 <- plotf('Older(Age≥65)')
p6 <- plotfy('Cardiovascular')
p7 <- plotf('Respiratory')
p8 <- plotf('Mental')
p9 <- plotf('Nervous')
p10 <- plotfx('Injury')
p11 <- plotf('Symptom')
library(patchwork)
(p1+ (p2+p3)/(p4+p5))/((p6+p7+p8)/(p9+p10+p11)) -> p
p
ggsave_pdf2png(p,'fig/ER_noninteract',width=18,height=15)

### Figure 3 Interaction plot -----------------------------------------------------------
library(ggpattern)
source('plot_theme.R')
rockthemes::muse_pal()(2)
rockthemes::deelite_pal()(2)
rockthemes::californication_pal()(2)
rockthemes::husker_pal()(2)
rockthemes::miles_pal()(2)
col <- c("#972D15","#08306B")
# col <- c("#7a6dbc","#dec2b0")
# col <- c("#88b1c1","#e0dcc2")
# col <- c("#739C9C","#EFDB15")
# col <- c("#739C9C","#e0dcc2")
# col <- c("#48448e","#fc4d97")
# col <- c("#cf80d5","#e4d4d2")
plotf <- function(outcome){
  # outcome='All cause'
  rbindlist(reslist) %>%
    dplyr::filter(outcome==!!outcome) -> tmp
  tmp %>%
    dplyr::select(predvar,rr=intRR_1,rrlow=intRRlow_1,rrhigh=intRRhigh_1) %>%
    mutate(type='Outage') -> tmp1
  tmp %>%
    dplyr::select(predvar,rr=intRR_0,rrlow=intRRlow_0,rrhigh=intRRhigh_0) %>%
    mutate(type='No outage') -> tmp2
  tmp <- rbind(tmp1,tmp2)
  tmp$type <- factor(tmp$type,levels=c('Outage','No outage'))
  ggplot(tmp) +
    geom_ribbon(aes(x=predvar,ymin=rrlow,ymax=rrhigh,fill=type),alpha=0.25) +
    scale_fill_manual(values = c("Outage"=col[1], "No outage"=col[2])) +
    geom_line(aes(x=predvar,y=rr,color=type),linewidth=1) +
    scale_color_manual(values = c("Outage"=col[1], "No outage"=col[2])) +
    scale_x_continuous(n.breaks=10) +
    geom_hline(yintercept=1,linetype='dashed')+
    labs(title=outcome,y='Relative risks(95% CI)',x=expression("Temperature ("~degree~"C)")) +
    theme_Publication()
}
plotfx <- function(outcome){
  # outcome='All cause'
  rbindlist(reslist) %>%
    dplyr::filter(outcome==!!outcome) -> tmp
  tmp %>%
    dplyr::select(predvar,rr=intRR_1,rrlow=intRRlow_1,rrhigh=intRRhigh_1) %>%
    mutate(type='Outage') -> tmp1
  tmp %>%
    dplyr::select(predvar,rr=intRR_0,rrlow=intRRlow_0,rrhigh=intRRhigh_0) %>%
    mutate(type='No outage') -> tmp2
  tmp <- rbind(tmp1,tmp2)
  tmp$type <- factor(tmp$type,levels=c('Outage','No outage'))
  ggplot(tmp) +
    geom_ribbon(aes(x=predvar,ymin=rrlow,ymax=rrhigh,fill=type),alpha=0.25) +
    scale_fill_manual(values = c("Outage"=col[1], "No outage"=col[2])) +
    geom_line(aes(x=predvar,y=rr,color=type),linewidth=1) +
    scale_color_manual(values = c("Outage"=col[1], "No outage"=col[2])) +
    scale_x_continuous(n.breaks=10) +
    geom_hline(yintercept=1,linetype='dashed')+
    labs(title=outcome,y='Relative risks(95% CI)',x=expression("Temperature ("~degree~"C)")) +
    theme_Publication_x()
}
plotfy <- function(outcome){
  # outcome='All cause'
  rbindlist(reslist) %>%
    dplyr::filter(outcome==!!outcome) -> tmp
  tmp %>%
    dplyr::select(predvar,rr=intRR_1,rrlow=intRRlow_1,rrhigh=intRRhigh_1) %>%
    mutate(type='Outage') -> tmp1
  tmp %>%
    dplyr::select(predvar,rr=intRR_0,rrlow=intRRlow_0,rrhigh=intRRhigh_0) %>%
    mutate(type='No outage') -> tmp2
  tmp <- rbind(tmp1,tmp2)
  tmp$type <- factor(tmp$type,levels=c('Outage','No outage'))
  ggplot(tmp) +
    geom_ribbon(aes(x=predvar,ymin=rrlow,ymax=rrhigh,fill=type),alpha=0.25) +
    scale_fill_manual(values = c("Outage"=col[1], "No outage"=col[2])) +
    geom_line(aes(x=predvar,y=rr,color=type),linewidth=1) +
    scale_color_manual(values = c("Outage"=col[1], "No outage"=col[2])) +
    scale_x_continuous(n.breaks=10) +
    geom_hline(yintercept=1,linetype='dashed')+
    labs(title=outcome,y='Relative risks(95% CI)',x=expression("Temperature ("~degree~"C)")) +
    theme_Publication_y()
}

p1 <- plotf('All cause')
p2 <- plotf('Male')
p3 <- plotf('Female')
p4 <- plotf('Younger(Age<65)')
p5 <- plotf('Older(Age≥65)')
p6 <- plotfy('Cardiovascular')
p7 <- plotf('Respiratory')
p8 <- plotf('Mental')
p9 <- plotf('Nervous')
p10 <- plotfx('Injury')
p11 <- plotf('Symptom')
library(patchwork)
(p1+ (p2+p3)/(p4+p5))/((p6+p7+p8)/(p9+p10+p11)) +
  plot_layout(guides = 'collect') &
  theme(legend.position = 'bottom')-> p
p
ggsave_pdf2png(p,'fig/ER_inter_0.7',width=18,height=15)

### Figure S RRR plot --------------------------------------------------------------------
source('plot_theme.R')
plotf <- function(outcome){
  # outcome='All cause'
  rbindlist(reslist) %>% 
    dplyr::filter(outcome==!!outcome) -> tmp
  ggplot(tmp) +
    geom_ribbon(aes(x=predvar,ymin=intRRlow,ymax=intRRhigh),fill='#E8E9EB')+
    geom_line(aes(x=predvar,y=intRR),
              color=ifelse(tmp$predvar>tmp$mmt,2,4),linewidth=1) +
    scale_x_continuous(n.breaks=10) +
    geom_hline(yintercept=1,linetype='dashed')+
    labs(title=outcome,y='Ratio of Relative risks(95% CI)',x=expression("Temperature ("~degree~"C)")) +
    theme_Publication(base_size=16) 
}
plotfx <- function(outcome){
  # outcome='All cause'
  rbindlist(reslist) %>% 
    dplyr::filter(outcome==!!outcome) -> tmp
  ggplot(tmp) +
    geom_ribbon(aes(x=predvar,ymin=intRRlow,ymax=intRRhigh),fill='#E8E9EB')+
    geom_line(aes(x=predvar,y=intRR),
              color=ifelse(tmp$predvar>tmp$mmt,2,4),linewidth=1) +
    scale_x_continuous(n.breaks=10) +
    geom_hline(yintercept=1,linetype='dashed')+
    labs(title=outcome,y='Ratio of Relative risks(95% CI)',x=expression("Temperature ("~degree~"C)")) +
    theme_Publication_x(base_size=16) 
}
plotfy <- function(outcome){
  # outcome='All cause'
  rbindlist(reslist) %>% 
    dplyr::filter(outcome==!!outcome) -> tmp
  ggplot(tmp) +
    geom_ribbon(aes(x=predvar,ymin=intRRlow,ymax=intRRhigh),fill='#E8E9EB')+
    geom_line(aes(x=predvar,y=intRR),
              color=ifelse(tmp$predvar>tmp$mmt,2,4),linewidth=1) +
    scale_x_continuous(n.breaks=10) +
    geom_hline(yintercept=1,linetype='dashed')+
    labs(title=outcome,y='Ratio of Relative risks(95% CI)',x=expression("Temperature ("~degree~"C)")) +
    theme_Publication_y(base_size=16) 
}
p1 <- plotf('All cause')
p2 <- plotf('Male')
p3 <- plotf('Female')
p4 <- plotf('Younger(Age<65)')
p5 <- plotf('Older(Age≥65)')
p6 <- plotfy('Cardiovascular')
p7 <- plotf('Respiratory')
p8 <- plotf('Mental')
p9 <- plotf('Nervous')
p10 <- plotfx('Injury')
p11 <- plotf('Symptom')
library(patchwork)
(p1+ (p2+p3)/(p4+p5))/((p6+p7+p8)/(p9+p10+p11)) -> p
p
ggsave_pdf2png(p,'fig/RRR_0.7',width=18,height=15)
 
### Figure 4 RR plot ---------------------------------------------------------------------
plotf <- function(outcome){
  # outcome='All cause'
  rbindlist(RRlist) %>%
    dplyr::filter(outcome==!!outcome) %>%
    mutate(perc=factor(perc,levels=c('1.0%','5.0%','10.0%','90.0%','95.0%','99.0%')),
           type=factor(type,levels=c('Outage','Non-outage'),labels=c('Outage','No outage'))) %>%
    ggplot() +
    geom_point(aes(x=perc,y=RR,shape=type),position=position_dodge(width=0.5),size=2.5) +
    geom_errorbar(aes(x=perc,y=RR,ymin=RRlow,ymax=RRhigh,group=type),
                  width=0.2,position=position_dodge(width=0.5)) +
    geom_hline(yintercept=1,linetype='dashed') +
    labs(title=outcome,y='Relative risks(95% CI)',x='Percentile of Temperature') +
    theme_Publication()
}
plotfx <- function(outcome){
  # outcome='All cause'
  rbindlist(RRlist) %>%
    dplyr::filter(outcome==!!outcome) %>%
    mutate(perc=factor(perc,levels=c('1.0%','5.0%','10.0%','90.0%','95.0%','99.0%')),
           type=factor(type,levels=c('Outage','Non-outage'),labels=c('Outage','No outage'))) %>%
    ggplot() +
    geom_point(aes(x=perc,y=RR,shape=type),position=position_dodge(width=0.5),size=2.5) +
    geom_errorbar(aes(x=perc,y=RR,ymin=RRlow,ymax=RRhigh,group=type),
                  width=0.2,position=position_dodge(width=0.5)) +
    geom_hline(yintercept=1,linetype='dashed') +
    labs(title=outcome,y='Relative risks(95% CI)',x='Percentile of Temperature') +
    theme_Publication_x()
}
plotfy <- function(outcome){
  # outcome='All cause'
  rbindlist(RRlist) %>%
    dplyr::filter(outcome==!!outcome) %>%
    mutate(perc=factor(perc,levels=c('1.0%','5.0%','10.0%','90.0%','95.0%','99.0%')),
           type=factor(type,levels=c('Outage','Non-outage'),labels=c('Outage','No outage'))) %>%
    ggplot() +
    geom_point(aes(x=perc,y=RR,shape=type),position=position_dodge(width=0.5),size=2.5) +
    geom_errorbar(aes(x=perc,y=RR,ymin=RRlow,ymax=RRhigh,group=type),
                  width=0.2,position=position_dodge(width=0.5)) +
    geom_hline(yintercept=1,linetype='dashed') +
    labs(title=outcome,y='Relative risks(95% CI)',x='Percentile of Temperature') +
    theme_Publication_y()
}
rbindlist(RRlist) %>%
  dplyr::filter(p.inter<0.05)

rbindlist(RRlist) %>% 
  mutate(RR=sprintf("%.2f (%.2f, %.2f)",RR,RRlow,RRhigh)) %>% 
  dplyr::select(outcome,perc,type,RR,p.inter) %>% 
  pivot_wider(names_from = type,values_from = RR) %>% 
  relocate(p.inter,.after=Outage) %>% 
  mutate(perc=factor(perc,levels=c('1.0%','5.0%','10.0%','90.0%','95.0%','99.0%')),
         outcome=factor(outcome,levels=c('All cause',
                                         'Female','Male',
                                         'Younger(Age<65)','Older(Age≥65)',
                                         'Cardiovascular','Respiratory',
                                         'Mental','Nervous',
                                         'Injury','Symptom'))) %>% 
  arrange(outcome,perc) %>% 
  openxlsx::write.xlsx('tab/RR.xlsx')


# 
p1 <- plotf('All cause')
p2 <- plotf('Male')
p3 <- plotf('Female')
p4 <- plotf('Younger(Age<65)')
p5 <- plotf('Older(Age≥65)')
p6 <- plotfy('Cardiovascular')
p7 <- plotf('Respiratory')
p8 <- plotf('Mental')
p9 <- plotf('Nervous')
p10 <- plotfx('Injury')
p11 <- plotf('Symptom')

p1 + geom_text(aes(x=4,y=1.05,label='*'),size=10,vjust=0,hjust=0.5) +
  geom_text(aes(x=5,y=1.1,label='*'),size=10,vjust=0,hjust=0.5) +
  geom_text(aes(x=6,y=1.3,label='*'),size=10,vjust=0,hjust=0.5) -> p1;p1
p2
p3 + geom_text(aes(x=4,y=1.05,label='*'),size=10,vjust=0,hjust=0.5) +
  geom_text(aes(x=5,y=1.1,label='*'),size=10,vjust=0,hjust=0.5) +
  geom_text(aes(x=6,y=1.3,label='*'),size=10,vjust=0,hjust=0.5) -> p3
p4
p5 + geom_text(aes(x=4,y=1.05,label='*'),size=10,vjust=0,hjust=0.5) +
  geom_text(aes(x=5,y=1.1,label='*'),size=10,vjust=0,hjust=0.5) +
  geom_text(aes(x=6,y=1.32,label='*'),size=10,vjust=0,hjust=0.5) -> p5
p6
p7
p8
p9
p11 + geom_text(aes(x=4,y=1.05,label='*'),size=10,vjust=0,hjust=0.5) +
  geom_text(aes(x=5,y=1.1,label='**'),size=10,vjust=0,hjust=0.5) +
  geom_text(aes(x=6,y=1.3,label='**'),size=10,vjust=0,hjust=0.5) -> p11
library(patchwork)
(p1+ (p2+p3)/(p4+p5))/((p6+p7+p8)/(p9+p10+p11)) +
  plot_layout(guides = 'collect') &
  theme(legend.position = 'bottom') -> p
p
ggsave_pdf2png(p,'fig/RR_0.7_percentile',width=18,height=15)