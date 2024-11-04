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
predper <- c(seq(0,1,0.1),2:98,seq(99,100,0.1))
set.seed(13041975)
tmeancountry <- round(rowMeans(sapply(dlist,function(x) quantile(jitter(x$tmean),
                                                                 predper/100,na.rm=T))),2)
# basis
lag=14
lagnk=2
argvar <- list(fun='ns',df=3)
cb <- crossbasis(df$tmean,lag=lag,argvar=argvar,
                 arglag=list(fun='ns',knots=logknots(lag,lagnk)),
                 group=df$code)

and1 <- list()
and2 <- list()
anlist <- list()
spldoy <- onebasis(df$doy,'ns',df=6)
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
             na.action="na.exclude",subset=keep,
             data=df2,family=quasipoisson)
  
  pred <- crosspred(cb,mod,
                    at=tmeancountry,
                    cen=15)
  
  red <- crossreduce(cb,mod,type='overall')
  coef <- coef(red)
  vcov <- vcov(red)
  
  mmt <- pred$predvar[11:109][which.min(pred$allRRfit[11:109])] # 1%-99%
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
  
  # grep('cb',names(coef(modint1)))
  ind <- 1:12
  redint1 <- crossreduce(cb,coef=coef(modint1)[ind],
                         vcov=vcov(modint1)[ind,ind],
                         model.link='log')
  coefint1 <- coef(redint1)
  vcovint1 <- vcov(redint1)
  
  redint2 <- crossreduce(cb,coef=coef(modint2)[ind],
                         vcov=vcov(modint2)[ind,ind],
                         model.link='log')
  coefint2 <- coef(redint2)
  vcovint2 <- vcov(redint2)
  
  # simulation
  nsim <- 1000
  library(MASS)
  coefsim <- mvrnorm(nsim,coef,vcov)
  coefsimint1 <- mvrnorm(nsim,coefint1,vcovint1)
  coefsimint2 <- mvrnorm(nsim,coefint2,vcovint2)
  
  # ATTRDL -----------------------------
  # in fact
  codes <- unique(df$code)
  attrtmp <- list()
  for(m in 1:length(codes)){
    cat(m,'')
    # m=1
    df %>%
      dplyr::filter(code==codes[m]) -> tmp
    
    death <- rowMeans(as.matrix(tsModel::Lag(tmp$y,-seq(0, 14))))
    
    argvar <- list(fun='ns',knots=attr(cb,'argvar')$knots,
                   Boundary.knots=attr(cb,'argvar')$Boundary.knots) #use fulltime knots
    bvar <- do.call(onebasis,c(list(x=tmp$tmean),argvar))
    cenvec <- do.call(onebasis,c(list(x=mmt),argvar))
    bvarcen <- scale(bvar,center=cenvec,scale=F)
    
    an1 <- (1-exp(-bvarcen%*%coefint1))*death
    an2 <- (1-exp(-bvarcen%*%coefint2))*death
    ind <- ifelse(tmp$tmean<mmt,T,F)
    
    po <- tmp$po
    poinflu <- rowMeans(tsModel::Lag(po,seq(0,14)))
    poind <- ifelse(poinflu>0,T,F)
    poind[is.na(poind)] <- F

    pocold <- ind & poind
    pohot <- (!ind) & poind
    poan1cold <- sum(an1[pocold],na.rm=T)
    poan1hot <- sum(an1[pohot],na.rm=T)
    poan2cold <- sum(an2[pocold],na.rm=T)
    poan2hot <- sum(an2[pohot],na.rm=T)
    
    an1cold <- sum(an1[ind],na.rm=T)
    an2cold <- sum(an2[ind],na.rm=T)
    an1hot <- sum(an1[!ind],na.rm=T)
    an2hot <- sum(an2[!ind],na.rm=T)
    
    apply(coefsimint1,1,function(coef){
      an <- (1-exp(-bvarcen%*%coef))*death
    }) -> an1sim
    apply(coefsimint2,1,function(coef){
      an <- (1-exp(-bvarcen%*%coef))*death
    }) -> an2sim
    
    
    if(sum(pohot)>1){
      cbind(total = colSums(an1sim,na.rm=T),
            cold = colSums(an1sim[ind,],na.rm=T),
            heat = colSums(an1sim[!ind,],na.rm=T),
            pocold = colSums(an1sim[pocold,],na.rm=T),
            pohot = colSums(an1sim[pohot,],na.rm=T)) -> an1sim
      cbind(total = colSums(an2sim,na.rm=T),
            cold = colSums(an2sim[ind,],na.rm=T),
            heat = colSums(an2sim[!ind,],na.rm=T),
            pocold = colSums(an2sim[pocold,],na.rm=T),
            pohot = colSums(an2sim[pohot,],na.rm=T)) -> an2sim
    }else{
      cbind(total = colSums(an1sim,na.rm=T),
            cold = colSums(an1sim[ind,],na.rm=T),
            heat = colSums(an1sim[!ind,],na.rm=T),
            pocold = colSums(an1sim[pocold,],na.rm=T),
            pohot = an1sim[pohot,]) -> an1sim
      cbind(total = colSums(an2sim,na.rm=T),
            cold = colSums(an2sim[ind,],na.rm=T),
            heat = colSums(an2sim[!ind,],na.rm=T),
            pocold = colSums(an2sim[pocold,],na.rm=T),
            pohot = an2sim[pohot,]) -> an2sim
    }

    
    ncount <- sum(death,na.rm=T)
    pocount <- sum(death[poind],na.rm=T)
    
    rbind(data.frame(city=codes[m],
                     type='total',
                     sim=c('est',paste0('sim',1:1000)),
                     an=c(an1hot+an1cold,an1sim[,1]),
                     ndeath=rep(ncount,1000+1)),
          data.frame(city=codes[m],
                     type='heat',
                     sim=c('est',paste0('sim',1:1000)),
                     an=c(an1hot,an1sim[,3]),
                     ndeath=rep(ncount,1000+1)),
          data.frame(city=codes[m],
                     type='cold',
                     sim=c('est',paste0('sim',1:1000)),
                     an=c(an1cold,an1sim[,2]),
                     ndeath=rep(ncount,1000+1)),
          data.frame(city=codes[m],
                     type='pocold',
                     sim=c('est',paste0('sim',1:1000)),
                     an=c(poan1cold,an1sim[,4]),
                     ndeath=rep(pocount,1000+1)),
          data.frame(city=codes[m],
                     type='pohot',
                     sim=c('est',paste0('sim',1:1000)),
                     an=c(poan1hot,an1sim[,5]),
                     ndeath=rep(pocount,1000+1))) -> an0
    rbind(data.frame(city=codes[m],
                     type='total',
                     sim=c('est',paste0('sim',1:1000)),
                     an=c(an2hot+an2cold,an2sim[,1]),
                     ndeath=rep(ncount,1000+1)),
          data.frame(city=codes[m],
                     type='heat',
                     sim=c('est',paste0('sim',1:1000)),
                     an=c(an2hot,an2sim[,3]),
                     ndeath=rep(ncount,1000+1)),
          data.frame(city=codes[m],
                     type='cold',
                     sim=c('est',paste0('sim',1:1000)),
                     an=c(an2cold,an2sim[,2]),
                     ndeath=rep(ncount,1000+1)),
          data.frame(city=codes[m],
                     type='pocold',
                     sim=c('est',paste0('sim',1:1000)),
                     an=c(poan2cold,an2sim[,4]),
                     ndeath=rep(pocount,1000+1)),
          data.frame(city=codes[m],
                     type='pohot',
                     sim=c('est',paste0('sim',1:1000)),
                     an=c(poan2hot,an2sim[,5]),
                     ndeath=rep(pocount,1000+1))) -> an1
    
    an0$kind <- 'cf0'
    an1$kind <- 'cf1'
    
    rbind(an0,an1) -> attrtmp[[m]]    
  }
  
  rbindlist(attrtmp) -> attrtmp
  attrtmpf <- attrtmp[, list(an=sum(an,na.rm=T), ndeath=sum(ndeath)),by=c('sim','type','kind')]
  attrtmpf %>%
    dplyr::filter(sim=='est') -> est
  
  attrtmpf %>%
    dplyr::filter(sim!='est') %>%
    group_by(type,kind) %>%
    summarise(anlow=quantile(an,0.025,na.rm=T),
              anhigh=quantile(an,0.975,na.rm=T)) -> sim
  est %>%
    left_join(sim,by=c('type','kind')) %>% 
    relocate(an,.after=ndeath) %>%
    mutate(outcome=outcome[i]) %>%
    relocate(outcome) %>%
    dplyr::select(-sim) -> an
  
  calculate_difference_and_ci <- function(type1, attrtmpf) {
    # type1='heat'
    x0 <- attrtmpf[type == type1 & (sim != 'est') & (kind=='cf0')]$an
    x1 <- attrtmpf[type == type1 & (sim != 'est') & (kind=='cf1')]$an
    
    if(type1=='heat'){
      potype='pohot'
    }else{
      potype='pocold'
    }
    
    x0po <- attrtmpf[type == potype & (sim != 'est') & (kind=='cf0')]$an
    x1po <- attrtmpf[type == potype & (sim != 'est') & (kind=='cf1')]$an
    
    xfact <- x0-x0po+x1po
    
    diff <- as.vector(outer(x1po, x0po, '-'))
    ci <- quantile(diff,c(0.025,0.975))
    
    af <- as.vector(outer(diff,xfact,'/')*100) 
    afci <- quantile(af,c(0.025,0.975))
    
    x0_e <- attrtmpf[type == type1 & (sim == 'est') & (kind=='cf0')]$an
    x1_e <- attrtmpf[type == type1 & (sim == 'est') & (kind=='cf1')]$an
    x0po_e <- attrtmpf[type == potype & (sim == 'est') & (kind=='cf0')]$an
    x1po_e <- attrtmpf[type == potype & (sim == 'est') & (kind=='cf1')]$an
    
    xfact_e <- x0_e - x0po_e + x1po_e
    diff_e <- x1po_e - x0po_e
    af_e <- diff_e/xfact_e*100
    
    data.frame(type=type1,
               fact=xfact_e,
               fact_low=xfact_e-ci[2],
               fact_high=xfact_e-ci[1],
               diff=diff_e,
               diff_low=ci[1],
               diff_high=ci[2],
               af=af_e,
               af_low=afci[1],
               af_high=afci[2]) -> tmp
    return(tmp)
  }
  
  rbind(
    calculate_difference_and_ci('heat', attrtmpf),
    calculate_difference_and_ci('cold', attrtmpf)
  ) %>% 
    mutate(outcome=outcome[i]) %>% 
    relocate(outcome) -> an_d
  
  calculate_difference_and_ci <- function(type1, attrtmpf) {
    # type1='pocold'
    
    x0po <- attrtmpf[type == type1 & (sim != 'est') & (kind=='cf0')]$an
    x1po <- attrtmpf[type == type1 & (sim != 'est') & (kind=='cf1')]$an
    
    diff <- as.vector(outer(x1po, x0po, '-'))
    ci <- quantile(diff,c(0.025,0.975))
    
    x0po_e <- attrtmpf[type == type1 & (sim == 'est') & (kind=='cf0')]$an
    x1po_e <- attrtmpf[type == type1 & (sim == 'est') & (kind=='cf1')]$an
    
    diff_e <- x1po_e - x0po_e
    af_e <- diff_e/x1po_e*100
    afci <- quantile(as.vector(outer(diff,x1po_e,'/')*100),c(0.025,0.975))
    
    data.frame(type=type1,
               diff=diff_e,
               diff_low=ci[1],
               diff_high=ci[2],
               af=af_e,
               af_low=afci[1],
               af_high=afci[2]) -> tmp
    return(tmp)
  }
  
  rbind(
    calculate_difference_and_ci('pohot', attrtmpf),
    calculate_difference_and_ci('pocold', attrtmpf)
  ) %>% 
    mutate(outcome=outcome[i]) %>% 
    relocate(outcome) -> an_d2
  
  # AN_diff
  anlist[[i]] <- an
  and1[[i]] <- an_d
  and2[[i]] <- an_d2
}

# limit the range and show the error bar
rbindlist(and1) -> and1
rbindlist(and2) -> and2
rbindlist(anlist) -> anlist
and1
and2


## Figure 5---------------------------------------------------------------------
source('plot_theme.R')
lower_limit <- -10
upper_limit <- 20
and1 %>% 
  as_tibble() %>% 
  dplyr::filter(type %in% c('cold','heat')) %>%
  dplyr::select(outcome,type,af,af_low,af_high) %>% 
  mutate(ll=ifelse(af_low<lower_limit,af-lower_limit,NA),
         ul=ifelse(af_high>upper_limit,upper_limit-af,NA)) %>% 
  mutate(outcome=factor(outcome,levels=(c('All cause','Male','Female',
                                          'Younger(Age<65)','Older(Age≥65)',
                                          'Cardiovascular','Respiratory',
                                          'Mental','Nervous','Injury','Symptom')))) %>% 
  ggplot(aes(x=outcome,y=af,color=type))+
  geom_point(position=position_dodge(width=0.5))+
  geom_errorbar(aes(ymin=af_low,ymax=af_high,color=type),position='dodge',width=0.5)+
  coord_cartesian(ylim=c(lower_limit,upper_limit))+
  labs(title='Full period',y='Attributable fraction (%)')+
  geom_hline(yintercept=0,linetype='dashed',color='black')+
  scale_color_manual(values=rev(c("#972D15","#08306B")))+
  scale_y_continuous(expand=c(0,0)) +
  # scale_x_discrete(expand=c(0,0)) +
  labs(x='Outcome',y='Attributable fraction (%)',color='')+
  theme_Publication_y() +
  theme(
    axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=12),
    axis.ticks.x = element_blank(),
    legend.position = 'none',
    legend.direction = 'vertical',
    legend.background = element_rect(fill = "transparent")) -> p1

source('plot_theme.R')
lower_limit <- -75
upper_limit <- 100
and2 %>% 
  as_tibble() %>% 
  dplyr::filter(type %in% c('pocold','pohot')) %>%
  mutate(type=ifelse(type=='pocold','Cold','Heat')) %>%
  dplyr::select(outcome,type,af,af_low,af_high) %>% 
  mutate(ll=ifelse(af_low<lower_limit,af-lower_limit,NA),
         ul=ifelse(af_high>upper_limit,upper_limit-af,NA)) %>% 
  mutate(outcome=factor(outcome,levels=(c('All cause','Male','Female',
                                          'Younger(Age<65)','Older(Age≥65)',
                                          'Cardiovascular','Respiratory',
                                          'Mental','Nervous','Injury','Symptom')))) %>% 
  ggplot(aes(x=outcome,y=af,color=type))+
  geom_point(position=position_dodge(width=0.5))+
  geom_errorbar(aes(ymin=af_low,ymax=af_high,color=type),position='dodge',width=0.5)+
  coord_cartesian(ylim=c(lower_limit,upper_limit))+
  labs(title='Power outage period',y='Attributable fraction (%)')+
  geom_hline(yintercept=0,linetype='dashed',color='black')+
  scale_color_manual(values=rev(c("#972D15","#08306B")))+
  scale_y_continuous(expand=c(0,0)) +
  # scale_x_discrete(expand=c(0,0)) +
  labs(x='Outcome',y='Attributable fraction (%)',color='')+
  theme_Publication_y() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = c(0.2,0.18),
    legend.direction = 'vertical',
    legend.background = element_rect(fill = "transparent")) -> p2

library(patchwork)
p <- p2/p1
ggsave_pdf2png(p,'fig/AF_ambulance',width=10, height=8)

################################################################################
data.frame(outcome=c('All cause','Female','Male',
                     'Younger(Age<65)','Older(Age≥65)',
                     'Cardiovascular','Respiratory',
                     'Mental','Nervous','Injury','Symptom'),
           ndeath=c(sum(df$count.all),sum(df$count.female),sum(df$count.male),
                    sum(df$count.age0_64),sum(df$count.age65_),
                    sum(df$count.cardio),sum(df$count.resp),
                    sum(df$count.mental),sum(df$count.nervous),
                    sum(df$count.injury),sum(df$count.symptom))) -> ndeath
and1 %>% 
  left_join(ndeath,by='outcome') %>% 
  mutate(AN_fact=sprintf('%.0f (%.0f, %.0f)',fact,fact_low,fact_high),
         AF_fact=sprintf('%.2f (%.2f, %.2f)',fact/ndeath*100,fact_low/ndeath*100,diff_high/ndeath*100),
         Po_addition=sprintf('%.0f (%.0f, %.0f)',diff,diff_low,diff_high),
         Po_addiAF=sprintf('%.2f (%.2f, %.2f)',af,af_low,af_high)) %>% 
  dplyr::select(outcome,type,AN_fact,AF_fact,Po_addition,Po_addiAF) -> aftable

library(openxlsx)
write.xlsx(aftable,'tab/AF_ambulance.xlsx')
