library(data.table)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(mgcv)
library(qgam)
library(visreg)

states <- c('California','Florida','Arizona','Illinois','Texas','Louisiana')
counties <- c('Los Angeles','Miami-Dade','Maricopa','Cook','Harris','East Baton Rouge')

excluded_states <- c('Guam','Northern Mariana Islands','Virgin Islands','Puerto Rico','Rhode Island')
### The first four states above are too small with data that isn't reported with the same reliability as the other US states.
## The trajectory of Rhode Island is too irregular, with massive Monday spikes in the second wave but not the first.
## This time-varying intensity of day-of-week effects produces poor fits of position & growth rates from the NBSS model.

find_peaks <- function(r,include_valleys=FALSE){
  y <- sign(r)
  sy <- y*shift(y)
  peak <- rep(0,length(y))
  if (include_valleys){
    peak[sy==-1] <- 1
  } else {
    peak[sy==-1 & y==(-1)] <- 1
  }
  return(peak)
}

prominence <- function(peak,z,log_scale=TRUE){
  ix <- which(peak==1)
  output <- rep(NA,length(z))
  if (length(ix)>0){
    mnz <- split(z,cumsum(peak)) %>% sapply(min)
    if (log_scale){
      output[ix] <- log(z[ix]/mnz[1:length(ix)])
    } else {
      output[ix] <- z[ix]-mnz[1:length(ix)]
    }
    return(output)
  } else {
    return(NA)
  }
}

# counties ----------------------------------------------------------------
US <- read.csv('data/nbss_us_counties.csv') %>% as.data.table()
US[,date:=as.Date(date)]


US[,peak:=find_peaks(growth_rate),by=c('state','county')]
US[,prominence_log:=prominence(peak,exp(mean_position)),by=c('state','county')]
US[,prominence:=prominence(peak,exp(mean_position),log_scale = FALSE),by=c('state','county')]
US[,lagged_dpc:=shift(deaths_pc,11,type='lead'),by=c('state','county')]
US[,polyname:=tolower(paste(state,county,sep=','))]

NYC <- US[county=='New York City' & date<as.Date('2020-05-01')]
largest_counties <- US[!state %in% excluded_states & county !='New York City',
                       list(county=unique(county[which.max(population)]),
                             population=max(population),
                            polyname=unique(polyname[which.max(population)])),by=state][order(population,decreasing=T)]



# Plotting cases, growth rates, and deaths pc at peak ---------------------

focal_polynames <- largest_counties$polyname[1:5]
cols <- viridis::viridis(length(focal_polynames))

X <- US[polyname %in% focal_polynames]
X[,nyc_growth_rate:=approx(NYC$lagged_dpc,NYC$growth_rate,xout=lagged_dpc)$y]
nyc_max <- US[peak==1 & county=='New York City'][prominence_log==max(prominence_log),lagged_dpc]




g_pk <- ggplot(X,aes(date,new_confirmed))+
  geom_bar(stat='identity',col='grey',fill='grey')+
  geom_line(aes(date,exp(mean_position)),col='black',lwd=2)+
  geom_point(data=X[peak==1],aes(color=county,shape=county,y=exp(mean_position),size=prominence_log))+
  theme_bw(base_size=14)+
  facet_wrap(.~county,scales = 'free_y',ncol=5)+
  theme(legend.position='none')+
  ggtitle('Peaks of large US counties')+
  scale_size_continuous(range = c(3,8))+
  scale_shape_manual(values=15:19)+
  scale_y_continuous('Cases')+
  scale_x_continuous(breaks=as.Date(c('2020-04-01','2020-07-01','2020-10-01','2021-01-01')),
                     labels=c('Apr','July','Oct','Jan'))

g_rD <- ggplot(X,aes(lagged_dpc,growth_rate))+
  geom_line(aes(y=nyc_growth_rate),col='red',lwd=2)+
  geom_line(aes(color=county),lwd=2)+
  theme_bw(base_size=14)+
  geom_hline(yintercept=0)+
  ggtitle('Growth Rate vs. Cumulative Incidence')+
  facet_wrap(.~county,ncol=5)+
  theme(legend.position='none')+
  scale_y_continuous('r(t)')+
  scale_x_continuous('D(t+11)',trans='log',limits=c(1e-5,4e-3),breaks=10^(-6:-3))+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

gdpc <- US[peak==1 & date<as.Date('2021-01-20')] %>% ggplot(aes(date,lagged_dpc))+
  geom_hline(yintercept=nyc_max,lwd=2,col='red')+
  geom_hline(yintercept=3.6e-3,lwd=2,lty=2)+
  geom_point(aes(size=prominence_log),pch=16,alpha=0.015)+
  geom_smooth(col='black')+
  geom_smooth(data=X[peak==1],lwd=2,alpha=0.5)+
  geom_point(data=X[peak==1],aes(color=county,shape=county,size=prominence_log))+
  theme_bw(base_size=14)+
  geom_point(data=US[county=='New York City' & prominence_log==max(prominence_log[county=='New York City'],na.rm=T)],pch=16,cex=10,col='red')+
  scale_y_continuous('D(t+11)',trans='log',breaks=10^(-5:-3),limits=c(1e-5,4e-3))+
  theme(legend.position='none')+
  ggtitle('Deaths per-capita at peak')+
  scale_shape_manual(values=15:19)+
  scale_size_continuous(range = c(5,10))
  
ggarrange(g_pk,g_rD,gdpc,nrow=3,heights = c(1,1,2),align='v')

ggsave('figures/large_counties_and_peak_D_trends.png',height=12,width=14,units='in')



# Regulation level --------------------------------------------------------

US[,month:=factor(paste(year(date),month(date),sep='-'),
                  levels=c(paste('2020',1:12,sep='-'),paste('2021',1:12,sep='-')))]


US[,week:=week(date)]
US[week<23,season:='Spring']
US[week>=23 & week<36,season:='Summer']
US[week>=36 & week<45,season:='Fall']
US[week>=45 | week<7,season:='Winter']
US[,season:=factor(season,levels=c('Spring','Summer','Fall','Winter'))]
US[,week:=factor(week(date),levels=c(7:52,1:6))]

US[,less_regulated_state:=state %in% c('Tennessee',"South Dakota",'South Carolina','Idaho','Georgia',
                                       'Mississippi','Oklahoma','Nebraska','Arizona','Texas')]
g_reg <- US[peak==1  & date>=as.Date('2020-04-01') & date<as.Date('2021-01-28')] %>% 
  ggplot(aes(month,log10(lagged_dpc)))+
  geom_hline(yintercept=log10(3.6e-3),lwd=2,lty=2)+
  geom_hline(yintercept=log10(nyc_max),col='red',lwd=2)+
  geom_boxplot(aes(color=less_regulated_state),lwd=1.2)+
  theme_bw(base_size=14)+
  scale_y_continuous('D(t+11)',breaks=-5:-3,labels=c('1e-5','1e-4','1e-3'),limits=c(-5.1,-2.4))+
  ggtitle('Deaths per-capita at peak, US counties')+
  guides(color=guide_legend(ncol=2))+
  theme(legend.position=c(0.8,0.2))

# season ------------------------------------------------------------------


g_season <- US[peak==1 & prominence_log>log(3)& date>=as.Date('2020-04-01') & !is.na(week) & date<as.Date('2021-01-28')] %>% 
  ggplot(aes(week,log10(lagged_dpc)))+
  # geom_jitter(aes(cex=population),alpha=0.3)+
  geom_hline(yintercept=log10(nyc_max),lwd=2,col='red')+
  geom_hline(yintercept=log10(3.6e-3),lwd=2,lty=2)+  
  geom_boxplot(aes(color=season),lwd=1.2)+
  theme_bw(base_size=14)+
  scale_y_continuous('D(t+11)',breaks=-5:-3,labels=c('1e-5','1e-4','1e-3'),limits=c(-5.1,-2.4))+
  # theme(legend.position='none')+
  ggtitle('Deaths per-capita at peak, all US counties')+
  guides(color=guide_legend(ncol=4))+
  theme(legend.position=c(0.8,0.18),legend.key = element_rect(fill = NA))+
  scale_color_manual(values=rev(viridis::plasma(4)))

ggarrange(g_reg,g_season,nrow=2,align='v')

ggsave('figures/peak_deaths_over_regulation_intensity_and_season.png',height=6,width=11,units='in')


# deaths per-capita at final peak -----------------------------------------

pops <- US[,list(pop=unique(population)),by=c('state','county')][order(pop,decreasing=T)]
pops[,pop_rank:=rank(1/pop,na.last = T)]
setkey(pops,state,county)
setkey(US,state,county)
US <- pops[US]



final_stop <-US[peak==1 & season=='Winter',list(growth_rate=growth_rate[.N],
                             lagged_dpc=lagged_dpc[.N],
                             polyname=unique(polyname),
                             less_regulated_state=unique(less_regulated_state),
                             population=population[1],
                             date=date[.N]),by=c('state','county')]
nrow(final_stop[!is.na(lagged_dpc)]) # 2864
length(unique(US$polyname)) #3146 total ---> 2/3146 = 
# 91% of US counties have peaked >=11 days before our final observation on February 8, 2021.
mds <- final_stop[,list(lagged_dpc=median(lagged_dpc,na.rm=T)),by=less_regulated_state]

g_rD_stop=US[polyname %in% largest_counties$polyname] %>% ggplot(aes(lagged_dpc,growth_rate))+
  geom_line(aes(group=county),lwd=2,alpha=0.08)+
  theme_bw(base_size=12)+
  geom_point(data=final_stop[polyname %in% largest_counties$polyname],col='black',cex=3)+
  geom_vline(data=final_stop[polyname %in% largest_counties$polyname],aes(xintercept=lagged_dpc),alpha=0.5)+
  geom_line(data=US[county=='New York City'],lwd=2,col='red')+
  geom_point(data=US[peak==1 & county=='New York City' & date<as.Date('2020-05-01')],col='red',cex=4)+
  scale_x_continuous('Deaths per-capita',trans='log',limits=c(1e-6,4e-3),breaks=10^(-6:-3))+
  theme(legend.position='none')+
  ggtitle('Largest County in Each State')
g_final_hist <- ggplot(final_stop,aes(lagged_dpc,group=less_regulated_state))+
  geom_histogram(bins=100,alpha=0.4,aes(fill=less_regulated_state),position='identity')+
  theme_bw(base_size=12)+
  scale_x_continuous('Deaths per-capita',trans='log',limits=c(4e-5,4e-3),breaks=10^(-6:-3))+
  geom_vline(data=mds,aes(xintercept=lagged_dpc,color=less_regulated_state),lwd=2)+
  geom_vline(xintercept = nyc_max,lwd=2,col='red')+
  geom_vline(xintercept=3.6e-3,lwd=2,lty=2)+
  theme(legend.position=c(0.28,0.72))+
  ggtitle('Final peak D*, histogram')



g_final_cdf <- ggplot(final_stop,aes(lagged_dpc))+
  stat_ecdf(aes(color=less_regulated_state),lwd=2,alpha=0.5)+
  scale_x_continuous('D*',trans='log',breaks=10^(-7:-3),limits=c(4e-5,4e-3))+
  geom_vline(xintercept = 3.6e-3,lwd=2,lty=2)+
  geom_vline(xintercept=nyc_max,lwd=2,col='red')+
  scale_y_continuous("P(D<D*)")+
  ggtitle('Final peak D*, CDF')+
  theme_bw(base_size=12)+
  geom_vline(data=mds,aes(xintercept=lagged_dpc,color=less_regulated_state),lwd=2)+
  theme(legend.position=c(0.28,0.72))



g_dynamics <- ggarrange(g_reg,g_season,nrow=2,align='v',labels = c("A","B"))
g_endpoints <- ggarrange(g_rD_stop,g_final_hist,g_final_cdf,nrow=3,ncol=1,labels=c('C',"D",'E'))
  
ggarrange(g_dynamics,g_endpoints,widths=c(3,1),ncol=2)
ggsave('figures/deaths_at_peaks_reg_season_largest_counties.png',height=8,width=16,units='in')



# Statistical analysis of regulations, populatoin size --------------------
final_stop[,deaths:=lagged_dpc*population]
final_stop[,survivors:=population-deaths]
final_stop[,datenum:=as.numeric(date)]
fit <- gam(cbind('Successes'=deaths,'Failures'=survivors)~s(datenum)+log10(population)+less_regulated_state,
           data=final_stop[!state %in% c('New York','New Jersey')],family=binomial)

summary(fit)

ilogit <- function(x) 1/(1+exp(-x))

ests <- cumsum(fit$coefficients[c('(Intercept)','less_regulated_stateTRUE')]) + 4*fit$coefficients['log10(population)']
ilogit(ests)

ilogit(ests)[2]/ilogit(ests[1]) ## +13% increase in deaths from COVID at fall peak in a 10,000 person county of a less-regulated state.
### This is a likely a limited view of the determinants of COVID mortality - state-specific effects, and more detailed county-level
### data may improve our understanding of the importance of mask-wearing and many other NPIs in place.

# Parametric coefficients:
#                             Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)              -7.542121   0.014026 -537.74   <2e-16 ***
#   log10(population)         0.124742   0.002474   50.42   <2e-16 ***
#   less_regulated_stateTRUE -0.143259   0.004258  -33.64   <2e-16 ***