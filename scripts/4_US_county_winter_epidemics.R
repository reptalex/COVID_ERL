library(data.table)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(mgcv)
library(qgam)
library(visreg)
library(usmap)
library(parallel)
ncores=7
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

# Load and preformat dataset  ----------------------------------------------------------------
US_counties <- read.csv('data/nbss_us_counties.csv') %>% as.data.table()
US_counties[,date:=as.Date(date)]

### filter to pre-February 2021 to remove effect of more transmissible variants

US_counties <- US_counties[date<as.Date('2021-02-01')]

US_counties[,peak:=find_peaks(growth_rate),by=c('state','county')]
US_counties[,prominence_log:=prominence(peak,exp(mean_position)),by=c('state','county')]
US_counties[,prominence:=prominence(peak,exp(mean_position),log_scale = FALSE),by=c('state','county')]
US_counties[,lagged_dpc:=shift(deaths_pc,11,type='lead'),by=c('state','county')]
US_counties[,polyname:=tolower(paste(state,county,sep=','))]

NYC <- US_counties[county=='New York City' & date<as.Date('2020-05-01')]
largest_counties <- US_counties[!state %in% excluded_states & county !='New York City',
                       list(county=unique(county[which.max(population)]),
                             population=max(population),
                            polyname=unique(polyname[which.max(population)])),by=state][order(population,decreasing=T)]



# Plotting cases, growth rates, and deaths pc at peak ---------------------

focal_polynames <- largest_counties$polyname[1:5]
cols <- viridis::viridis(length(focal_polynames))

X <- US_counties[polyname %in% focal_polynames]
X[,nyc_growth_rate:=approx(NYC$lagged_dpc,NYC$growth_rate,xout=lagged_dpc)$y]
nyc_max <- US_counties[peak==1 & county=='New York City' & date<as.Date('2020-07-01'),lagged_dpc]

X[date>as.Date('2020-08-01'),county_max:=max(prominence_log,na.rm=T),by=polyname]
X[peak==1 & prominence_log==county_max,c('polyname','date','lagged_dpc')]


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

gdpc <- US_counties[peak==1 & date<as.Date('2021-01-20')] %>% ggplot(aes(date,lagged_dpc))+
  geom_hline(yintercept=nyc_max,lwd=2,col='red')+
  geom_hline(yintercept=3.6e-3,lwd=2,lty=2)+
  geom_point(aes(size=prominence_log),pch=16,alpha=0.015)+
  geom_smooth(col='black')+
  geom_smooth(data=X[peak==1],lwd=2,alpha=0.5)+
  geom_point(data=X[peak==1],aes(color=county,shape=county,size=prominence_log))+
  theme_bw(base_size=14)+
  geom_point(data=US_counties[peak==1 & county=='New York City' & date<as.Date('2020-07-01')],pch=16,cex=10,col='red')+
  # geom_line(data=US_counties[county=='New York City'],col='red')+
  # geom_point(data=US_counties[county=='New York City' & peak==1 & prominence > 2000],col='red',cex=4)+ ## shows the full NYC trajectory & fall peak
  scale_y_continuous('D(t+11)',trans='log',breaks=10^(-5:-3),limits=c(1e-5,4e-3))+
  theme(legend.position='none')+
  ggtitle('Deaths per-capita at peak')+
  scale_shape_manual(values=15:19)+
  scale_size_continuous(range = c(5,10))
  
ggarrange(g_pk,g_rD,gdpc,nrow=3,heights = c(1,1,2),align='v',labels = c('A','B','C'))

ggsave('figures/large_counties_and_peak_D_trends.png',height=12,width=14,units='in')



# Season & Mask Wearing --------------------------------------------------------

US_counties[,month:=factor(paste(year(date),month(date),sep='-'),
                  levels=c(paste('2020',1:12,sep='-'),paste('2021',1:12,sep='-')))]


US_counties[,week:=week(date)]
US_counties[week<23,season:='Spring']
US_counties[week>=23 & week<36,season:='Summer']
US_counties[week>=36 & week<45,season:='Fall']
US_counties[week>=45 | week<7,season:='Winter']
US_counties[,season:=factor(season,levels=c('Spring','Summer','Fall','Winter'))]
US_counties[,week:=factor(week(date),levels=c(7:52,1:6))]

# https://abcnews.go.com/Health/states-dropped-mask-mandates/story?id=76249857 ### states below had no mask-wearing mandates as of march 4
US_counties[,mask_wearing:= ! state %in% c('Arizona','Florida','Georgia','Idaho','Missouri','Nebraska','Oklahoma',
                                           'South Carolina','South Dakota','Tennessee')]


US_counties[peak==1  & date>=as.Date('2020-04-01') & prominence_log>log(3)] %>% 
  ggplot(aes(month,log10(lagged_dpc)))+
  geom_hline(yintercept=log10(3.6e-3),lwd=2,lty=2)+
  geom_hline(yintercept=log10(swe$lagged_dpc),lwd=2,col=swy)+
  geom_hline(yintercept=log10(nyc_max),col='red',lwd=2)+
  geom_boxplot(aes(color=mask_wearing),lwd=1.2)+
  theme_bw(base_size=14)+
  scale_y_continuous('D(t+11)',breaks=-5:-3,labels=c('1e-5','1e-4','1e-3'),limits=c(-5.1,-2.4))+
  ggtitle('D* across US counties by Mask Wearing Mandate')+
  guides(color=guide_legend(ncol=2))+
  theme(legend.position=c(0.8,0.2))

final_stop[,month:=month(date)]

ggplot(final_stop,aes(date,log10(lagged_dpc)))+
  geom_hline(yintercept=log10(3.6e-3),lwd=2,lty=2)+
  geom_hline(yintercept=log10(swe$lagged_dpc),lwd=2,col=swy)+
  geom_hline(yintercept=log10(nyc_max),col='red',lwd=2)+
  geom_point(aes(color=mask_wearing),cex=2,alpha=0.2)+
  geom_smooth(aes(color=mask_wearing))+
  theme_bw(base_size=14)+
  scale_y_continuous('D(t+11)',breaks=-4:-3,labels=c('1e-4','1e-3'),limits=c(-4.1,-2.4))+
  ggtitle('D* across US counties by Mask Wearing Mandate')+
  guides(color=guide_legend(ncol=2))+
  theme(legend.position=c(0.8,0.2))

SWE <- read.csv('data/nbss_countries.csv') %>% as.data.table
SWE <- SWE[country=='Sweden' & date<as.Date('2020-09-01')]
SWE[,lagged_dpc:=shift(deaths/population,11,type='lead')]
SWE[,peak:=find_peaks(growth_rate)]
SWE[,prominence:=prominence(peak,exp(mean_position))]

swe <- SWE[peak==1 & prominence==max(prominence,na.rm=T)]
swe[,date]  ## June 17th peak
swe[,lagged_dpc] ### 0.0005345342 deaths per-capita
swy <- rgb(1,205/255,0) ## swedish yellow

g_reg <- US_counties[peak==1  & date>=as.Date('2020-08-01') & prominence_log>log(3) & !week %in% c('4','5',NA)] %>% 
  ggplot(aes(week,log10(lagged_dpc)))+
  geom_hline(yintercept=log10(3.6e-3),lwd=2,lty=2)+
  geom_hline(yintercept=log10(swe$lagged_dpc),lwd=2,col=swy)+
  geom_hline(yintercept=log10(nyc_max),col='red',lwd=2)+
  geom_boxplot(aes(color=mask_wearing),lwd=1.2)+
  theme_bw(base_size=14)+
  scale_y_continuous('D(t+11)',breaks=-5:-3,labels=c('1e-5','1e-4','1e-3'),limits=c(-5.1,-2.4))+
  ggtitle('D* across US counties by Mask Wearing Mandate')+
  guides(color=guide_legend(ncol=2))+
  theme(legend.position=c(0.8,0.2))


# season ------------------------------------------------------------------


g_season <- US_counties[peak==1 & prominence_log>log(3)& date>=as.Date('2020-04-01') & !is.na(week) & !week %in% c('4','5')] %>% 
  ggplot(aes(week,log10(lagged_dpc)))+
  # geom_jitter(aes(cex=population),alpha=0.3)+
  geom_hline(yintercept=log10(nyc_max),lwd=2,col='red')+
  geom_hline(yintercept=log10(swe$lagged_dpc),lwd=2,col=swy)+
  geom_hline(yintercept=log10(3.6e-3),lwd=2,lty=2)+  
  geom_boxplot(aes(color=season),lwd=1.2)+
  theme_bw(base_size=14)+
  scale_y_continuous('D(t+11)',breaks=-5:-3,labels=c('1e-5','1e-4','1e-3'),limits=c(-5.1,-2.4))+
  # theme(legend.position='none')+
  ggtitle('Shifting endpoints under seasonal forcing')+
  guides(color=guide_legend(ncol=4))+
  theme(legend.position=c(0.8,0.18),legend.key = element_rect(fill = NA))+
  scale_color_manual(values=rev(viridis::magma(4)))+
  geom_text(aes(x=factor(15,levels=c(7:52,1:6)),y=log10(nyc_max*1.5),label='NYC'),cex=8,col='red')+
  geom_point(aes(x=factor(week(NYC[peak==1]$date),levels=c(7:52,1:6)),y=log10(nyc_max)),col='red',cex=6)+
  geom_point(aes(x=factor(week(swe$date),levels=c(7:52,1:6)),y=log10(swe$lagged_dpc)),cex=6,col=swy)+
  geom_text(aes(x=factor(week(swe$date),levels=c(7:52,1:6)),y=log10(swe$lagged_dpc*1.5),label='SWE'),cex=6,col=swy)


# deaths per-capita at final peak -----------------------------------------

pops <- US_counties[,list(pop=unique(population)),by=c('state','county')][order(pop,decreasing=T)]
pops[,pop_rank:=rank(1/pop,na.last = T)]
setkey(pops,state,county)
setkey(US_counties,state,county)
US_counties <- pops[US_counties]



final_stop <-US_counties[peak==1 & season=='Winter' & prominence_log>log(3),list(growth_rate=growth_rate[.N],
                             lagged_dpc=lagged_dpc[.N],
                             polyname=unique(polyname),
                             mask_wearing=unique(mask_wearing),
                             population=population[1],
                             date=date[.N]),by=c('state','county')]
nrow(final_stop[!is.na(lagged_dpc)]) # 2498
length(unique(US_counties$polyname)) #3146 total
# 91% of US counties have peaked >=11 days before our final observation on February 8, 2021.
mds <- final_stop[,list(lagged_dpc=mean(lagged_dpc,na.rm=T)),by=mask_wearing]

g_rD_stop=US_counties[polyname %in% largest_counties$polyname] %>% ggplot(aes(lagged_dpc,growth_rate))+
  geom_line(aes(group=county),lwd=2,alpha=0.08)+
  theme_bw(base_size=12)+
  geom_point(data=final_stop[polyname %in% largest_counties$polyname],col='black',cex=3)+
  geom_vline(data=final_stop[polyname %in% largest_counties$polyname],aes(xintercept=lagged_dpc),alpha=0.5)+
  geom_line(data=US_counties[county=='New York City'],lwd=2,col='red')+
  geom_point(data=US_counties[peak==1 & county=='New York City' & date<as.Date('2020-05-01')],col='red',cex=4)+
  scale_x_continuous('D(t+11)',trans='log',limits=c(1e-6,4e-3),breaks=10^(-6:-3))+
  scale_y_continuous('r(t)')+
  theme(legend.position='none')+
  ggtitle('Largest County in Each State')
g_final_hist <- ggplot(final_stop,aes(lagged_dpc,group=mask_wearing))+
  geom_histogram(bins=100,alpha=0.4,aes(fill=mask_wearing),position='identity')+
  theme_bw(base_size=12)+
  scale_x_continuous('D*',trans='log',limits=c(4e-5,4e-3),breaks=10^(-6:-3))+
  geom_vline(data=mds,aes(xintercept=lagged_dpc,color=mask_wearing),lwd=2)+
  geom_vline(xintercept = nyc_max,lwd=2,col='red')+
  geom_vline(xintercept=3.6e-3,lwd=2,lty=2)+
  theme(legend.position=c(0.28,0.72))+
  ggtitle('Deaths per-capita, final peak')

# Vector field inference --------------------------------------------------
US_counties[,dr:=shift(growth_rate,type='lead')-growth_rate,by=c('state','county')]
US_counties[,dD:=log(shift(lagged_dpc,type='lead')/lagged_dpc),by=c('state','county')]
US_counties[,woy:=as.numeric(strftime(date,format='%V'))]
US_counties[,ci_width:=p97.5_growth_rate-p2.5_growth_rate]


ili <- cdcfluview::ilinet(region = 'state') %>% as.data.table()
ili_avg <- ili[,list(hist_ili=mean(ilitotal/total_patients,na.rm=T)),by=week]
ili_avg[,woy:=week]
ili_avg[which.max(hist_ili)] ### week 52 had the maximum historical average of patients visiting doc with ILI
wk=12  ### we'll predict the ERLs for the a priori estimated time of approximate maximum respiratory viral transmission across US counties
setkey(ili_avg,woy)
setkey(US_counties,woy)
US_counties <- US_counties[ili_avg[,c('woy','hist_ili')]]

gr <- gam(dr~hist_ili+te(growth_rate,log(lagged_dpc)),
          data=US_counties[lagged_dpc>0])
gD <- gam(dD~hist_ili+te(growth_rate,log(lagged_dpc)),
          data=US_counties[!is.infinite(dD)])


n=10   ### number of vectors to draw on our plot
grd <- expand.grid('woy'=wk,
                   'hist_ili'=max(ili_avg$hist_ili),
                   'lagged_dpc'=10^seq(-5,log10(4e-3),length.out = n),
                   'growth_rate'=seq(-0.1,0.6,length.out=n)) %>% as.data.table


grd$dr <- predict(gr,newdata=grd)
grd$dD <- predict(gD,newdata=grd)
grd[dD<0,dD:=NA]

nyc <- US_counties[county=='New York City' & date<as.Date('2020-05-01')]
nyc[,logD:=log(lagged_dpc)]
grd[,logD:=log(lagged_dpc)]   ## by log transform, making vector fields with vectors of fixed length is easier.
grd <- grd[growth_rate-approx(nyc$lagged_dpc,nyc$p97.5_growth_rate,xout = grd$lagged_dpc)$y<0.2]

### Below, we impute the points on the NYC line for a finer estimation of where the ERL is rejected
ny_impute <- data.table('logD'=seq(min(nyc[!is.infinite(logD)]$logD,na.rm=T),max(nyc$logD,na.rm=T),length.out=1e3),
                        'hist_ili'=max(ili_avg$hist_ili))
ny_impute[,lagged_dpc:=exp(logD)]
ny_impute[,woy:=wk]
ny_impute[,growth_rate:=approx(nyc$lagged_dpc,nyc$growth_rate,xout = lagged_dpc)$y]
ny_impute[,slp:=(shift(growth_rate,type='lead')-growth_rate)/(shift(logD,type='lead')-logD)]
ny_impute$pred_slope <- predict(gr,newdata=ny_impute)/predict(gD,newdata=ny_impute)
ny_impute[,slope_diff:=pred_slope-slp] ### If this is positive, the vector field crosses the NYC line from below at this point.
ny_impute[,ERL:=c('ERL','Rejected')[as.numeric(slope_diff>0)+1]]


### unit vectors: (dd,rr) in the (log(D),r) phase space.
grd[,dd:=dD/sqrt((dD^2+dr^2))]
grd[,rr:=dr/sqrt((dD^2+dr^2))]

# g_vec <- ggplot(grd,aes(logD,growth_rate))+
#   geom_segment(aes(xend=logD+dd/n,yend=growth_rate+rr/n),arrow=arrow(length = unit(0.04,'inches')))+
#   geom_line(data=nyc,lwd=2,col='red')+
#   geom_ribbon(data=nyc,aes(ymin=p2.5_growth_rate,ymax=p97.5_growth_rate),fill='red',alpha=0.2)+
#   scale_x_continuous('D(t+11)',breaks=log(10^(-5:-2)),limits=log(c(1e-5,4e-3)),labels = paste('1e-',5:2))+
#   scale_y_continuous('r(t)',limits=c(-0.1,0.6))+
#   theme_bw(base_size=12)+
#   geom_vline(xintercept = log(3.6e-3),lty=2,lwd=2)+
#   geom_point(data=ny_impute[ERL=='Rejected'],aes(color=ERL),cex=2)+
#   geom_hline(yintercept = 0)+
#   theme(legend.position=c(0.75,0.8))+
#   scale_color_manual(values='yellow')+
#   ggtitle('Expected r(D) Trajectories')

### add lines of expected trajectories
integrate_rd <- function(r0=0.5,D0=1e-5,Dmax=2e-3,woy=52,hist_ili=0.03782747,rmin=-0.1,gr.=gr,gD.=gD){
  dd <- data.table('woy'=woy,'hist_ili'=hist_ili,'growth_rate'=r0,'lagged_dpc'=D0)
  D <- D0
  r <- r0
  while(D<Dmax & r>rmin){
    dr <- predict(gr,newdata=dd[.N])
    dD <- predict(gD,newdata=dd[.N])
    r <- r+dr
    D <- D*exp(dD)
    dd <- rbind(dd,data.table('woy'=woy,'hist_ili'=hist_ili,'growth_rate'=r,'lagged_dpc'=D))
  }
  dd[,logD:=log(lagged_dpc)] %>% return
}

### compute trajectories
cl <- makeCluster(ncores)
clusterEvalQ(cl,expr={library(data.table)
                      library(magrittr)
                      library(mgcv)})
clusterExport(cl,varlist=c('gr','gD'))

r_sampled <- seq(0,1,length.out = 60)
traj_US <- parLapply(cl,r_sampled,integrate_rd,woy=wk)
names(traj_US) <- 1:length(r_sampled)
trajcs <- rbindlist(traj_US,use.names=T,idcol='trajectory')

stopCluster(cl)
rm('cl')
gc()

### find maximum fall growth rates across counties
mx_fall_growth <- US_counties[date>as.Date('2020-08-01'),list(r=max(growth_rate,na.rm=T),
                                                              D=lagged_dpc[which.max(growth_rate)],
                                                              dt=date[which.max(growth_rate)]),by=c('state','county')]
pt=mx_fall_growth[,list(growth_rate=median(r[!is.infinite(r)],na.rm=T),
                     logD=median(log(D[D>0]),na.rm=T),
                     dt=median(dt,na.rm=T))]
### find median

### October 28 was the median date of maximum growth across US counties.
ix <- trajcs[,min((pt$growth_rate-growth_rate)^2+(pt$logD-logD)^2),by=trajectory][,trajectory[which.min(V1)]]
nyc[,trajectory:='NYC ERL']

dum <- rbind(nyc[,c('growth_rate','logD','trajectory')],
             trajcs[trajectory %in% c(1,ix),c('growth_rate','logD','trajectory')])
dum[trajectory==1,trajectory:='predicted']
dum[trajectory==ix,trajectory:='US fall max']

g_vec <- ggplot(grd,aes(logD,growth_rate))+
  geom_line(data=trajcs,aes(group=trajectory),alpha=0.5)+
  geom_segment(aes(xend=logD+dd/n,yend=growth_rate+rr/n),arrow=arrow(length = unit(0.04,'inches')))+
  geom_line(data=nyc,lwd=2,col='red')+
  theme_bw(base_size=12)+
  geom_vline(xintercept = log(3.6e-3),lty=2,lwd=2)+
  geom_line(data=dum,aes(size=trajectory,color=trajectory))+
  geom_ribbon(data=nyc,aes(ymin=p2.5_growth_rate,ymax=p97.5_growth_rate),fill='red',alpha=0.4)+
  scale_color_manual(values=c('red','black','black'))+
  scale_size_manual(values=c(1.5,1,1.5))+
  # geom_point(data=ny_impute[ERL=='Rejected'],aes(color=ERL),cex=2)+
  geom_hline(yintercept = 0)+
  theme(legend.position=c(0.82,0.7))+
  ggtitle('Expected r(D) Trajectories')+
  scale_x_continuous('D(t+11)',breaks=log(10^(-5:-2)),limits=log(c(1e-5,4e-3)),labels = paste('1e-',5:2))+
  scale_y_continuous('r(t)',limits=c(-0.1,0.6))+
  geom_line(data=trajcs[trajectory==ix,],lwd=2)+
  geom_point(data=pt,cex=4)

g_vec
ggsave('figures/r_D_inferred_vector_field.png',height=10,width=10,units='in')

# combining figures -------------------------------------------------------

g_dynamics <- ggarrange(g_season,g_reg,nrow=2,align='v',labels = c("A","B"))
g_endpoints <- ggarrange(g_rD_stop,g_vec,g_final_hist,nrow=3,ncol=1,labels=c('C',"D",'E'))
  
ggarrange(g_dynamics,g_endpoints,widths=c(2.7,1),ncol=2)
ggsave('figures/deaths_at_peaks_reg_season_largest_counties.png',height=9,width=17,units='in')


# mask-wearing ------------------------------------------------------------

final_stop[,deaths:=lagged_dpc*population]
final_stop[,survivors:=population-deaths]
final_stop[,datenum:=as.numeric(date)]

fit <- gam(cbind('Successes'=deaths,'Failures'=survivors)~s(log10(population))+mask_wearing,
           data=final_stop[state!='New York'],family=binomial)

### We exclude New York State because nearly 75% of the deaths in New York occured prior to the implementation of the mask-wearing mandate.

fit$coefficients['mask_wearingTRUE']
# mask_wearingTRUE 
# -0.1320571
### mask-wearing reduced deaths at final stop by roughly 13%... we predict this more fully below

counter_factual <- final_stop[,c('state','county','population','mask_wearing','deaths','lagged_dpc')]
counter_factual[,mask_wearing:=!mask_wearing]

counter_factual$counterfactual_preds <- predict(fit,newdata = counter_factual,type='response')
counter_factual[,lives_saved_by_masks:=(counterfactual_preds-lagged_dpc)*population]
counter_factual[,mask_wearing:=!mask_wearing]

#### summarize counter-factuals:
x <- counter_factual[state!='New York',list(total_lives_saved=sum(lives_saved_by_masks,na.rm=T),
                      total_deaths=sum(deaths,na.rm=T)),by=mask_wearing]
x
#    mask_wearing total_lives_saved total_deaths
# 1:         TRUE         29002.240       205654
# 2:        FALSE         -7699.463        62295

### Estimated 29,000 lives saved by the peak from NPIs colinear with mask-wearing mandates.
### Estimated 7,700 lives lost in states without NPIs colinear with mask-wearing mandates.

x[,total_lives_saved/(total_deaths)]
# 0.1410244 -0.1235968
# This model estimates an 14% in deaths had states outside New York not implemented mask-wearing mandates and colinear NPIs.
# Model estimates 12.4% of deaths in states without mask-wearing mandates would have been saved by mask-wearing mandates.