rm(list=ls())
gc()
library(data.table)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(COVID19)
source('scripts/utils.R')


US <- read.csv('data/nbss_us_counties.csv') %>% as.data.table()
US[,date:=as.Date(date)]

mn <- US[,list(growth_rate=mean(growth_rate,na.rm=T)),by=date]
mn[,growing:=factor(sign(growth_rate),levels=c(-1,1))]
US[,polyname:=tolower(paste(state,county,sep=','))]

setkey(US,polyname)

#### Every County - time intensive plot.
g_c <- ggplot(US[date>=as.Date('2020-03-01')],aes(date,growth_rate))+
  geom_line(aes(group=polyname),alpha=0.01)+
  geom_line(data=mn,lwd=2)+
  theme_bw()+
  scale_y_continuous('Growth Rate, r(t)',limits=c(-0.12,0.25))+
  theme(legend.position='none')+
  geom_hline(yintercept=0)+
  scale_x_continuous(breaks=seq(as.Date('2020-03-01'),as.Date('2021-02-01'),by='month'),
                     labels=c('Mar','Apr','May','June','July','Aug','Sep','Oct','Nov','Dec','Jan','Feb'))+
  ggtitle('Growth rates across US Counties')


# ILI ---------------------------------------------------------------------



nat <- cdcfluview::ilinet('national') %>% as.data.table()
nat_avg <- nat[year!=2020,list(age_65=sum(age_65),
                               total_patients=sum(total_patients)),by=week]
ili_mn=nat_avg[,age_65/total_patients]
ili_mx=nat[year==2020,age_65/total_patients]


ili <- cdcfluview::ilinet('state') %>% as.data.table()

ili[,ryear:=paste(region,year,sep='_')]
avg_ili <- ili[,list(unweighted_ili=mean(unweighted_ili,na.rm=T)),by=week]
approx_date <- ili[year==2019 & region=='New York',list(date=as.Date(week_start)+368),by=week]

ix <- round(seq(1,52,length.out = 10))
lims <- approx_date[date>=as.Date('2020-02-27') & date<=max(US$date),list(m=c(min(week),max(week)))]
g_ili <- ggplot(ili,aes(week,unweighted_ili))+
  geom_line(aes(group=ryear),col='black',lwd=2,alpha=0.01)+
  geom_line(data=avg_ili,lwd=2,col='black')+
  theme_bw()+
  scale_y_continuous('Percent ILI',trans='log',breaks=c(0.005,0.01,0.02,0.04,0.08,0.16)*100,limits=c(0.25,20))+
  scale_x_continuous('Approximate 2020 date',breaks=approx_date$week[ix],labels = approx_date$date[ix],
                     limits=lims$m)+
  ggtitle('ILI over all regions since 2010')



# Combine -----------------------------------------------------------------



ggarrange(g_c+geom_vline(xintercept=as.Date('2020-08-13'),col='red'),
          g_ili+geom_vline(xintercept=33,col='red'),nrow=2,align='v')
ggsave('figures/rt_vs_ILI_seasonality.png',height=12,width=10,units='in')
