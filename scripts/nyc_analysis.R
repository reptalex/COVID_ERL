library(data.table)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(viridis)
source('scripts/utils.R')

load('data/fits.RData')


fits <- as.data.table(fits)

US_counties <- fits[administrative_area_level==3 & administrative_area_level_1=='United States']
US_counties[,state:=administrative_area_level_2]
US_counties[,county:=administrative_area_level_3]
US_counties[,lagged_dpc:=shift(deaths,11,type='lead')/population,by=c('state','county')]
US_counties[,polyname:=tolower(paste(state,county,sep=','))]



# A view of early epidemics -----------------------------------------------
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


Lombardy <- fits[administrative_area_level==2 & administrative_area_level_2=='Lombardia']
Lombardy[,lagged_dpc:=shift(deaths,11,type='lead')/population]
Lombardy[,region:='Lombardy, Italy']

cols <- c('date','growth_rate','lagged_dpc','region','new_confirmed','new_deaths')
X <- rbind(NYC[,cols,with=F],
           EssexNJ[,cols,with=F],
           Boston[,cols,with=F],
           Detroit[,cols,with=F],
           Nawlins[,cols,with=F],
           Lombardy[,cols,with=F])
X[,region:=factor(region,levels = c(setdiff(unique(region),'New York City'),'New York City'))]


rD <- ggplot(X[date<as.Date('2020-05-01')],aes(lagged_dpc,growth_rate,col=region))+
  geom_hline(yintercept = 0)+
  geom_line(lwd=2)+
  scale_x_continuous(trans='log',limits=c(1e-5,2e-3),breaks=10^(-5:-3))+
  theme_bw(base_size=15)+
  theme(legend.position=c(0.8,0.8))+
  scale_color_manual(values=c(viridis(5),'red'))

Nt <- ggplot(X[date<as.Date('2020-05-01')],aes(date,new_confirmed,col=region,fill=region))+
  geom_bar(stat='identity')+
  geom_bar(aes(y=new_deaths),stat='identity',col='black',fill='black')+
  theme_bw(base_size=15)+
  facet_wrap(.~region,nrow=3,ncol=2,scales = 'free_y')+
  theme(legend.position='none')+
  scale_y_continuous(trans='log',breaks=round(10^(seq(0,3,by=0.5))))+
  scale_color_manual(values=c(viridis(5),'red'))+
  scale_fill_manual(values=c(viridis(5),'red'))

ggarrange(Nt,rD,ncol=2,widths = c(1,1))


# NYC analysis ------------------------------------------------------------

### using information prior to May 1, 2020

nyc <- NYC[date<as.Date('2020-05-01'),c('date','new_confirmed','new_deaths','deaths','population')]


nyc[,lockdown:=c('pre-lockdown','lockdown')[as.numeric(date>=as.Date('2020-03-23'))+1]]
nyc[,lockdown:=factor(lockdown,levels=c('pre-lockdown','lockdown'))]


nyc[,growth_rate:=nbs(new_confirmed)]
nyc[,growth_rate_deaths:=nbs(new_deaths)]
nyc[,lagged_dpc:=shift(deaths,11,type='lead')/population]
nyc[,deaths_pc:=deaths/population]

g_cd <- ggplot(nyc[lockdown=='pre-lockdown'],aes(date,new_confirmed))+
  geom_bar(stat='identity',col='red',fill='red')+
  geom_bar(aes(y=new_deaths),stat='identity',col='black',fill='black')+
  geom_vline(xintercept = as.Date('2020-03-22'))+
  theme_bw(base_size=15)+
  scale_y_continuous('N(t)',trans='log',breaks=10^(0:4))+
  ggtitle('Cases (red), Deaths (black)')


# New York City projections -----------------------------------------------
### get credible intervals etc
x <- nyc[lockdown=='pre-lockdown']

x$q2.5=nbs(x$new_deaths,'p2.5_growth_rate')
x$q97.5=nbs(x$new_confirmed,'p97.5_growth_rate')
r_pre_ld <- x[,list('q2.5'=mean(q2.5,na.rm=T),
                    'r'=mean(growth_rate_deaths,na.rm=T),
                    'q97.5'=mean(q97.5,na.rm=T))]
r_pre_ld[,sd:=(q97.5-q2.5)/(qnorm(.975)-qnorm(.025))]

project_exp <- function(deaths,r,IFR=0.01,time_exposure_to_death=20) (deaths/IFR)*exp(time_exposure_to_death*r)

### assume r is gamma distributed, 
### IFR beta distributed mean 0.06 ...

### simpler plot: show exponential projection to "HIT=0.6" as a boundary line in plot of r vs. IFR, and then show distribution of estimates of r, IFR.
### and/or include IFR, r, tau_D in three pairwise plots, fixing other at its median estimate.
