library(data.table)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(viridis)
library(mgcv)
library(ggforce)
source('scripts/utils.R')
nbs <- function(x,...) nbss(x,...)$growth_rate
US_counties <- read.csv('data/nbss_us_counties.csv') %>% as.data.table()
US_counties[,date:=as.Date(date)]
lockdown_date <- as.Date('2020-03-22')

NYC <- US_counties[county=='New York City']
nyc <- NYC[date<as.Date('2020-05-01'),c('date','new_confirmed','new_deaths','deaths','population')]
nyc[,lagged_dpc:=shift(deaths,11,type='lead')/population]


### Cases & Deaths pre-lockdown
g_cd <- ggplot(nyc[date<=lockdown_date],aes(date,new_confirmed))+
  geom_bar(stat='identity',col='red',fill='red')+
  geom_bar(aes(y=new_deaths),stat='identity',col='black',fill='black')+
  geom_vline(xintercept = as.Date('2020-03-22'))+
  theme_bw(base_size=15)+
  scale_y_continuous('N(t)',trans='log',breaks=10^(0:4))+
  ggtitle('Cases (red), Deaths (black)')




# Exposed projections -----------------------------------------------------

lockdown_deaths <- NYC[date==lockdown_date,new_deaths]
ifr <- 0.014
death_grs <- cbind(NYC[,c('date')],nbss(NYC$new_deaths))
r <- death_grs[date<=lockdown_date,mean(growth_rate,na.rm=T)]
r2.5 <- death_grs[date<=lockdown_date,mean(p2.5_growth_rate,na.rm=T)]
r97.5 <- death_grs[date<=lockdown_date,mean(p97.5_growth_rate,na.rm=T)]


HypothNYC <- data.table('date'=seq(as.Date('2020-01-15'),as.Date('2020-03-22'),by='day'))
HypothNYC[,death_date:=date+20]
HypothNYC[,I:=(lockdown_deaths/ifr)*exp(r*(as.numeric(death_date-lockdown_date)))]
HypothNYC[,I2.5:=(lockdown_deaths/ifr)*exp(r2.5*(as.numeric(death_date-lockdown_date)))]
HypothNYC[,I97.5:=(lockdown_deaths/ifr)*exp(r97.5*(as.numeric(death_date-lockdown_date)))]

HypothNYC[,new_confirmed:=I]
HypothNYC <- HypothNYC[I>1]

HypothNYC[,est_deaths:=predict(deaths_fit,newdata=.SD,type='response')]


g_hypoth <- ggplot(nyc[date<=lockdown_date],aes(date,new_confirmed))+
  geom_ribbon(data=HypothNYC,aes(xmin=min(HypothNYC$date),xmax=max(HypothNYC$date),
                  ymin=nyc$population[1]*.6,ymax=max(HypothNYC$I97.5)),fill='grey')+
  geom_hline(yintercept=nyc$population[1]*0.6,lty=2,lwd=2)+
  geom_ribbon(data=HypothNYC[date>=lockdown_date-20],aes(ymin=I2.5,ymax=I97.5),fill='orange',alpha=0.5)+
  geom_line(data=HypothNYC,aes(y=I),col='orange',lwd=2)+
  geom_bar(stat='identity',col='red',fill='red')+
  geom_bar(aes(y=new_deaths),stat='identity',col='black',fill='black')+
  theme_bw(base_size=15)+
  # geom_hline(yintercept = nyc$population[1],lwd=2)+
  # geom_ribbon(data=HypothNYC[I97.5>nyc$population[1]*0.5],aes(ymin=nyc$population[1]*0.6,ymax=I97.5),fill='orange')+
  scale_y_continuous('N(t)',trans='log',breaks=10^(0:11),limits=c(1,max(HypothNYC$I97.5)))+
  ggtitle('NYC COVID Epidemic, pre-lockdown')+
  geom_segment(aes(x = lockdown_date, y = lockdown_deaths, 
                   xend = lockdown_date-20, yend = lockdown_deaths),
               arrow = arrow(length = unit(0.5, "cm")),lwd=2)+
  geom_segment(aes(x = lockdown_date-20, y = lockdown_deaths,
                   xend = lockdown_date-20, yend = HypothNYC[date==(lockdown_date-20),I]),
               arrow = arrow(length = unit(0.5, "cm")),lwd=2)+
  annotate(geom='text',x=as.Date('2020-03-09'),y=lockdown_deaths*2,
           label=expression(paste(tau,'=20',sep='')),size=8)+
  annotate(geom='text',x=lockdown_date-15.5,y=320,
           label='IFR=0.01',size=8)+
  annotate(geom='text',x=as.Date('2020-03-01'),y=nyc$population[1]*1.3,
           label='HIT',size=8)
g_hypoth


# IFR - r and HIT boundary plot -------------------------------------------

r_mean <- death_grs[date<=lockdown_date,mean(growth_rate,na.rm=T)]
r_sd <- death_grs[date<=lockdown_date,mean(p75_growth_rate-p25_growth_rate,na.rm=T)/(qnorm(.75)-qnorm(.25))]
rset <- seq(r_mean-3*r_sd,r_mean+3*r_sd,length.out = 100)
rd <- dnorm(rset,r_mean,r_sd)

R_df <- data.table('growth_rate'=rset,'density'=rd)


nyc_hit <- function(r,lockdown_deaths.=lockdown_deaths,N=nyc$population[1],HIT=0.6,lag_exposure_to_deaths=20){
  (lockdown_deaths/(0.6*N))*exp(r*lag_exposure_to_deaths)
}

R_df[,HIT:=nyc_hit(growth_rate)]

####### IFR quantile matching with log-norml
## https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(20)30769-6/fulltext
### IFR 95% CI: 1.04 - 1.77, median 1.39
ql <- function(x,qs=c(0.0104,0.0139,0.0177)) sum((qlnorm(c(0.025,0.5,0.975),meanlog=x[1],sdlog=x[2])-qs)^2)
res <- optim(par = c(0.0139,0.002),fn=ql)
ifr_median <- 0.0139
ifr_sd <- res$par[2]
IFR_df <- data.table('IFR'=seq(min(R_df$HIT),max(R_df$HIT),length.out=1000))
IFR_df[,density:=dlnorm(IFR,meanlog=res$par[1],sdlog=res$par[2])]


blnk <- ggplot()+theme_void()

g_r_density <- ggplot(R_df,aes(growth_rate,density))+
  geom_line()+
  geom_polygon(fill='steelblue',alpha=0.4)+
  theme_bw(base_size=12)+
  scale_x_continuous(name=NULL)+
  ggtitle("Growth rate density")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(1,0.1,-1,1),'cm'))

n=1e4
dum <- data.table('growth_rate'=rnorm(n,mean=r_mean,sd=r_sd),
                  'IFR1'=rlnorm(n,meanlog = res$par[1],sdlog = res$par[2]),
                  'IFR2'=rlnorm(n,meanlog = res$par[1],sdlog=res$par[2])/2)

g_HIT <- ggplot(R_df,aes(growth_rate,HIT))+
  geom_line(lwd=2,lty=2)+
  geom_ribbon(aes(ymin=1.1e-3,ymax=HIT),fill='grey')+
  theme_bw(base_size=12)+
  stat_ellipse(data=dum,aes(y=IFR1),geom='polygon',fill='purple',alpha=0.4,lwd=2)+
  stat_ellipse(data=dum,aes(y=IFR2),geom='polygon',fill='green',alpha=0.4,lwd=2)+
  scale_x_continuous('Growth Rate')+
  geom_point(aes(x=r_mean,y=ifr_median),fill='purple',pch=21)+
  geom_point(aes(x=r_mean,y=ifr_median/2),fill='green',pch=21)+
  geom_line()+
  # geom_vline(xintercept = r_mean,col='steelblue',lwd=2)+
  # geom_hline(yintercept = ifr_median,col='red',lwd=2)+
  # geom_hline(yintercept = ifr_median/2,col='green',lwd=2)+
  annotate(geom='text',x = 0.45,y=0.003,label='HIT',size=10)+
  annotate(geom='text',x=0.3,y=0.1,label='Avoid HIT',size=10)+
  scale_y_continuous('IFR',trans='log',breaks=10^seq(-4,-.5,by=0.5),
                     labels = signif(10^seq(-4,-.5,by=0.5),2),
                     limits=c(1e-3,0.3))+
  theme(plot.margin = unit(c(-2,0.1,1,1),'cm'))

ifr_density <- ggplot(IFR_df,aes(IFR,density))+
  geom_line()+
  geom_polygon(fill='purple',alpha=0.4)+
  geom_line(aes(x=IFR/2))+
  geom_polygon(fill='green',aes(x=IFR/2),alpha=0.4)+
  theme_bw(base_size=12)+
  scale_x_continuous(name=NULL,trans='log',breaks=10^seq(-4,-.5,by=0.5),
                     labels = signif(10^seq(-4,-.5,by=0.5),2),
                     limits=c(1e-3,0.3))+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  coord_flip()+
  ggtitle("IFR density")

d=4
g_NY_HIT <- ggarrange(ggarrange(blnk,g_r_density,blnk,ncol=3,widths=c(0.1,4,1),align='h'),
          ggarrange(g_HIT,ifr_density,ncol=2,widths=c(d,1),align='h'),
          nrow=2,heights=c(1,d),align='v')


ggarrange(g_hypoth,g_NY_HIT,ncol=2,labels = c("A",'B') )
ggsave('figures/nyc_hit_plausibility.png',height=8,width=17,units='in')



## percent HIT
n=1e5
dum <- data.table('growth_rate'=rnorm(n,mean=r_mean,sd=r_sd),
                  'IFR1'=rlnorm(n,meanlog = res$par[1],sdlog = res$par[2]),
                  'IFR2'=rlnorm(n,meanlog = res$par[1],sdlog=res$par[2])/2,
                  'IFR3'=rlnorm(n,meanlog = log(0.01),sdlog=res$par[2])/2)

dum[,sum((lockdown_deaths/IFR1)*exp(growth_rate*20)>0.6*nyc$population[1])]/n
dum[,sum((lockdown_deaths/IFR2)*exp(growth_rate*20)>0.6*nyc$population[1])]/n
dum[,sum((lockdown_deaths/IFR3)*exp(growth_rate*20)>0.6*nyc$population[1])]/n


# Other early epidemics ---------------------------------------------------

US_counties[,lagged_dpc:=shift(deaths,11,type='lead')/population,by=c('state','county')]
NYC <- US_counties[county=='New York City']
NYC[,region:=county]
EssexNJ <- US_counties[state=='New Jersey' & county=='Essex']
EssexNJ[,region:="Essex, NJ"]
Boston <- US_counties[state=='Massachusetts' & county=='Suffolk']
Boston[,region:='Suffolk, MA']
Detroit <- US_counties[state=='Michigan' & county=='Wayne']
Detroit[,region:="Detroit, MI"]
Nawlins <- US_counties[state=='Louisiana' & county=='Orleans']
Nawlins[,region:='New Orleans']

load('data/fits.RData')
fits <- as.data.table(fits)
Lombardy <- fits[administrative_area_level==2 & administrative_area_level_2=='Lombardia']
Lombardy[,lagged_dpc:=shift(deaths,11,type='lead')/population]
Lombardy[,region:='Lombardy, Italy']

cols <- c('date','growth_rate','lagged_dpc','region','new_confirmed','new_deaths','p2.5_growth_rate','p97.5_growth_rate')
X <- rbind(NYC[,cols,with=F],
           EssexNJ[,cols,with=F],
           Boston[,cols,with=F],
           Detroit[,cols,with=F],
           Nawlins[,cols,with=F],
           Lombardy[,cols,with=F])
X[,region:=factor(region,levels = c(setdiff(unique(region),'New York City'),'New York City'))]


g_rD <- ggplot(X[date<as.Date('2020-05-01')],aes(lagged_dpc,growth_rate,col=region))+
  geom_hline(yintercept = 0)+
  geom_line(lwd=2)+
  geom_vline(xintercept = 0.006,lty=2,lwd=2)+
  geom_ribbon(data=X[date<as.Date('2020-05-01') & region=='New York City'],
              aes(ymin=p2.5_growth_rate,ymax=p97.5_growth_rate),fill='red',alpha=0.2,
              show.legend=FALSE)+
  scale_x_continuous('D(t+11)',trans='log',limits=c(1e-5,6e-3),breaks=10^(-5:-2))+
  scale_y_continuous('r(t)')+
  theme_bw(base_size=15)+
  theme(legend.position=c(0.75,0.7))+
  scale_color_manual(values=c(viridis(5),'red'))+
  ggtitle("Timescale of deaths")


ggarrange(ggarrange(g_hypoth,g_NY_HIT,ncol=2,labels = c("A",'B')),
                    g_rD,labels=c(NA,'C'),nrow=2,heights=c(1.5,1))
ggsave('figures/nyc_hit_plausibility_2.png',height=12,width=16,units='in')
