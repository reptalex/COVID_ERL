rm(list=ls())
gc()
library(data.table)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
source('scripts/utils.R')

rDlag <- 11
max_date <- as.Date('2020-09-01')
swy <- rgb(1,205/255,0)   ## color of swedish line

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

# import & prepare data --------------------------------------------------------
USA <- read.csv('data/nbss_us_states.csv') %>% as.data.table
USA[,date:=as.Date(date)]
USA <- USA[date<max_date]
USA[,peak:=find_peaks(growth_rate),by=state]
USA[,prominence:=prominence(peak,exp(mean_position),log_scale = FALSE),by=state]
USA[,biggest_peak:= as.numeric(prominence==max(prominence,na.rm=T)),by=state]

SWE[,peak:=find_peaks(growth_rate)]
SWE[,prominence:=prominence(peak,exp(mean_position),log_scale=FALSE)]
SWE[prominence==max(prominence,na.rm=T)]   ### June 17th

bgs <- USA[biggest_peak==TRUE]
### looking at the peaks - this works well for the six states we're focusing on. It fails for others, but since the purpose is to show these dates occur after Sweden's
### June 17th peak, I'm going to leave this.
# ggplot(USA,aes(date,new_confirmed))+
#   geom_line()+
#   facet_wrap(.~state,scales = 'free_y')+
#   geom_vline(data=bgs,aes(xintercept=date),col='red')+
#   theme_bw()

setkey(USA,state,date)
X <- read.csv('data/nbss_countries.csv') %>% as.data.table

SWE <- X[country=='Sweden']
SWE[,date:=as.Date(date)]
SWE <- SWE[date<max_date]

nyc <- read.csv('data/nbss_nyc.csv') %>% as.data.table()
nyc[,date:=as.Date(date)]
nyc <- nyc[date<as.Date('2020-05-01')]

USA[,deaths_pc:=shift(deaths_pc,rDlag,type='lead'),by=state]
SWE[,deaths_pc:=shift(deaths_pc,rDlag,type='lead')]
nyc[,deaths_pc:=shift(deaths_pc,rDlag,type='lead')]

AZ <- USA[state=='Arizona',]
FL <- USA[state=='Florida']
TX <- USA[state=='Texas']
GA <- USA[state=='Georgia']
CA <- USA[state=='California']
LA <- USA[state=='Louisiana']
NY <- USA[state=='New York']

# State interventions -----------------------------------------------------

intervention_cols <- data.frame('intervention'=c('none','stay-at-home',
                                                 'reopening1','reopening2','masks','reversal'),
                                'color'=brewer.pal(6,'Paired'))

Tots <- rbind(AZ,FL,TX,GA,CA,LA)
Tots[,intervention:=NA]

# Arizona -----------------------------------------------------------------
AZ[date<as.Date('2020-03-30'),intervention:='none']
AZ[date>=as.Date('2020-03-30') & date<as.Date('2020-05-15'),intervention:='stay-at-home']
AZ[date>=as.Date('2020-05-15') & date<as.Date('2020-06-17'),intervention:='reopening']
AZ[date>=as.Date('2020-06-17'),intervention:='masks']
AZ[date>as.Date('2020-06-29'),intervention:='reversal']
AZ[,intervention:=factor(intervention,levels=unique(intervention))]

lbl <- data.frame('deaths_pc'=c(1e-5,3e-6,1e-4),
                  'growth_rate'=c(0.37,0.1,.42),
                  'label'=c('NY State','Sweden','NYC'))

az <- ggplot(AZ,aes(deaths_pc,growth_rate,color=intervention))+
#   geom_line(data=Tots,alpha=0.2,lwd=2,aes(group=state))+
  geom_line(data=nyc,aes(color='black'),color='red',lwd=2)+
  geom_line(data=NY,aes(color='black'),color='black',lwd=2)+
  geom_line(data=SWE,aes(color='black'),color=swy,lwd=2)+
  geom_line(lwd=3)+
  # geom_line(data=AZ[intervention=='reversal'],color='red',lwd=2)+
  scale_x_continuous(trans='log',limits=c(1e-6,3e-3),breaks=10^(-6:0),labels = NULL,name=NULL)+
  scale_y_continuous(limits=c(-0.1,0.5))+
  theme_bw(base_size=15)+
  scale_color_manual(values=intervention_cols$color[c(1,2,4,5,6)])+
  geom_text(data=lbl,aes(label=label,color=NULL))+
  geom_label(data=lbl,aes(label=label,color=NULL),col=c('black',swy,'red'),size=10)+
  geom_hline(yintercept = 0)+
  ggtitle(paste('Arizona: peaked',bgs[state=='Arizona',date]))+
  theme(legend.position = c(.5,.04),
        legend.direction = 'horizontal',
        legend.background = element_blank(),
        legend.box.background = element_rect(color = 'black'))

az.rt <- ggplot(AZ[date>as.Date('2020-03-01')],
                aes(date,growth_rate,color=intervention))+
  geom_line(lwd=1.5)+
  theme_bw()+
  scale_color_manual(values=intervention_cols$color[c(1:3,5,6)])+
  theme(legend.position='none')+
  geom_vline(xintercept = AZ[grepl('local',intervention),min(date)])+
  scale_y_continuous('r(t)')+  
  scale_x_continuous(limits = as.Date(c('2020-03-01','2020-09-01')),
                     breaks=as.Date(c('2020-03-01','2020-05-01','2020-07-01','2020-09-01')),
                     labels = c('Mar','May','July','Sep'))

az.cs <- ggplot(AZ[date>as.Date('2020-03-01')],
                aes(date,new_confirmed,fill=intervention))+
  geom_bar(stat='identity')+
  theme_bw()+
  scale_fill_manual(values=intervention_cols$color[c(1:3,5,6)])+
  theme(legend.position='none')+
  scale_y_continuous('N(t)')+
  scale_x_continuous(name=NULL,breaks=NULL)

az2=az+annotation_custom(ggplotGrob(az.rt),
                         xmin=log(3e-4),xmax=log(3e-3),
                         ymin=0.17,ymax=0.35)+
  annotation_custom(ggplotGrob(az.cs),
                    xmin=log(2.8e-4),xmax=log(3e-3),
                    ymin=0.35,ymax=0.5)


# Texas -------------------------------------------------------------------

TX[date<as.Date('2020-03-31'),intervention:='none']
TX[date>=as.Date('2020-03-31') & date<as.Date('2020-05-01'),intervention:='stay-at-home']
TX[date>=as.Date('2020-05-01') & date<as.Date('2020-05-18'),intervention:='Phase 1']
# TX[date>=as.Date('2020-05-18') & date<as.Date('2020-06-03'),intervention:='Phase 2']
# TX[date>=as.Date('2020-06-03') & date<as.Date('2020-06-26'),intervention:='Phase 3'] ### optional, but these increase figure complexity
TX[date>=as.Date('2020-05-18') & date<as.Date('2020-06-26'),intervention:='Phase 2-3']
TX[date>=as.Date('2020-06-26'),intervention:='reversal']
TX[,intervention:=factor(intervention,levels=unique(intervention))]

tx <- ggplot(TX,aes(deaths_pc,growth_rate,color=intervention))+
#   geom_line(data=Tots,alpha=0.2,lwd=2,aes(group=state))+
  geom_line(data=nyc,aes(color='black'),color='red',lwd=2)+
  geom_line(data=NY,aes(color='black'),color='black',lwd=2)+
  geom_line(data=SWE,aes(color='black'),color=swy,lwd=2)+
  geom_line(lwd=3)+
  # geom_line(data=TX[intervention=='reversal'],color='red',lwd=2)+
  scale_x_continuous(trans='log',limits=c(1e-6,3e-3),breaks=10^(-6:0),labels = NULL,name=NULL)+
  scale_y_continuous(limits=c(-0.1,0.5),labels = NULL,name=NULL)+
  theme_bw(base_size=15)+
  scale_color_manual(values=intervention_cols$color[c(1:4,6)])+
  geom_hline(yintercept = 0)+
  ggtitle(paste('Texsa: peaked',bgs[state=='Texas',date]))+
  theme(legend.position = c(.5,0.04),
        legend.direction = 'horizontal',
        legend.background = element_blank(),
        legend.box.background = element_rect(color = 'black'))

tx.rt <- ggplot(TX[date>as.Date('2020-03-01')],
                aes(date,growth_rate,color=intervention))+
  geom_line(lwd=1.5)+
  theme_bw()+
  scale_color_manual(values=intervention_cols$color[c(1:4,6)])+
  theme(legend.position='none')+
  geom_vline(xintercept = TX[grepl('reversal',intervention),min(date)])+
  scale_y_continuous('r(t)')+  
  scale_x_continuous(limits = as.Date(c('2020-03-01','2020-09-01')),
                     breaks=as.Date(c('2020-03-01','2020-05-01','2020-07-01','2020-09-01')),
                     labels = c('Mar','May','July','Sep'))

tx.cs <- ggplot(TX[date>as.Date('2020-03-01')],
                aes(date,new_confirmed,fill=intervention))+
  geom_bar(stat='identity')+
  theme_bw()+
  scale_fill_manual(values=intervention_cols$color[c(1:4,6)])+
  theme(legend.position='none')+
  scale_y_continuous('N(t)')+
  scale_x_continuous(name=NULL,breaks=NULL)

tx2=tx+annotation_custom(ggplotGrob(tx.rt),
                         xmin=log(3e-4),xmax=log(3e-3),
                         ymin=0.17,ymax=0.35)+
  annotation_custom(ggplotGrob(tx.cs),
                    xmin=log(2.65e-4),xmax=log(3e-3),
                    ymin=0.35,ymax=0.5)

# Louisiana ---------------------------------------------------------------
LA[date<as.Date('2020-03-23'),intervention:='none']
LA[date>=as.Date('2020-03-23') & date<as.Date('2020-05-15'),intervention:='stay-at-home']
LA[date>=as.Date('2020-05-15') & date<as.Date('2020-06-05'),intervention:='Phase 1']
LA[date>=as.Date('2020-06-05'),intervention:='Phase 2']

LA[,intervention:=factor(intervention,levels=unique(intervention))]

la <- ggplot(LA,aes(deaths_pc,growth_rate,color=intervention))+
#   geom_line(data=Tots,alpha=0.2,lwd=2,aes(group=state))+
  geom_line(data=nyc,aes(color='black'),color='red',lwd=2)+
  geom_line(data=NY,aes(color='black'),color='black',lwd=2)+
  geom_line(data=SWE,aes(color='black'),color=swy,lwd=2)+
  geom_line(lwd=3)+
  scale_x_continuous(trans='log',limits=c(1e-6,3e-3),breaks=10^(-6:0))+
  scale_y_continuous(limits=c(-0.1,0.5),labels = NULL,name=NULL)+
  theme_bw(base_size=15)+
  scale_color_manual(values=intervention_cols$color[1:4])+
  geom_hline(yintercept = 0)+
  ggtitle(paste('Louisiana: peaked',bgs[state=='Louisiana',date]))+
  theme(legend.position = c(.5,0.04),
        legend.direction = 'horizontal',
        legend.background = element_blank(),
        legend.box.background = element_rect(color = 'black'))

la.rt <- ggplot(LA[date>as.Date('2020-03-01')],
                aes(date,growth_rate,color=intervention))+
  geom_line(lwd=1.5)+
  theme_bw()+
  scale_color_manual(values=intervention_cols$color[1:4])+
  theme(legend.position='none')+
  scale_y_continuous('r(t)')+  
  scale_x_continuous(limits = as.Date(c('2020-03-01','2020-09-01')),
                     breaks=as.Date(c('2020-03-01','2020-05-01','2020-07-01','2020-09-01')),
                     labels = c('Mar','May','July','Sep'))

la.cs <- ggplot(LA[date>as.Date('2020-03-01')],
                aes(date,new_confirmed,fill=intervention))+
  geom_bar(stat='identity')+
  theme_bw()+
  scale_fill_manual(values=intervention_cols$color[1:4])+
  theme(legend.position='none')+
  scale_y_continuous('N(t)')+
  scale_x_continuous(name=NULL,breaks=NULL)

la2 <- la+annotation_custom(ggplotGrob(la.rt),
                            xmin=log(3e-4),xmax=log(3e-3),
                            ymin=0.17,ymax=0.35)+
  annotation_custom(ggplotGrob(la.cs),
                    xmin=log(2.8e-4),xmax=log(3e-3),
                    ymin=0.35,ymax=0.5)


# Florida -----------------------------------------------------------------
FL[date<as.Date('2020-03-30'),intervention:='none']
FL[date>=as.Date('2020-03-30') & date<as.Date('2020-06-01'),intervention:='stay-at-home']
FL[date>=as.Date('2020-06-01'),intervention:='reopening']
FL[,intervention:=factor(intervention,levels=unique(intervention))]
fl <- ggplot(FL,aes(deaths_pc,growth_rate,color=intervention))+
#   geom_line(data=Tots,alpha=0.2,lwd=2,aes(group=state))+
  geom_line(data=nyc,aes(color='black'),color='red',lwd=2)+
  geom_line(data=NY,aes(color='black'),color='black',lwd=2)+
  geom_line(data=SWE,aes(color='black'),color=swy,lwd=2)+
  geom_line(lwd=3)+
  scale_x_continuous(trans='log',limits=c(1e-6,3e-3),breaks=10^(-6:0),labels = NULL,name=NULL)+
  scale_y_continuous(limits=c(-0.1,0.5))+
  theme_bw(base_size=15)+
  scale_color_manual(values=intervention_cols$color[c(1,2,4)])+
  geom_hline(yintercept = 0)+
  ggtitle(paste('Florida: peaked',bgs[state=='Florida',date]))+
  theme(legend.position = c(.5,0.04),
        legend.direction = 'horizontal',
        legend.background = element_blank(),
        legend.box.background = element_rect(color = 'black'))
  

fl.rt <- ggplot(FL[date>as.Date('2020-03-01')],
                aes(date,growth_rate,color=intervention))+
  geom_line(lwd=1.5)+
  theme_bw()+
  scale_color_manual(values=intervention_cols$color[c(1,2,4)])+
  theme(legend.position='none')+
  scale_y_continuous('r(t)')+  
  scale_x_continuous(limits = as.Date(c('2020-03-01','2020-09-01')),
                     breaks=as.Date(c('2020-03-01','2020-05-01','2020-07-01','2020-09-01')),
                     labels = c('Mar','May','July','Sep'))

fl.cs <- ggplot(FL[date>as.Date('2020-03-01')],
                aes(date,new_confirmed,fill=intervention))+
  geom_bar(stat='identity')+
  theme_bw()+
  scale_fill_manual(values=intervention_cols$color[c(1,2,4)])+
  theme(legend.position='none')+
  scale_y_continuous('N(t)')+
  scale_x_continuous(name=NULL,breaks=NULL)

fl2=fl+annotation_custom(ggplotGrob(fl.rt),
                         xmin=log(3e-4),xmax=log(3e-3),
                         ymin=0.17,ymax=0.35)+
  annotation_custom(ggplotGrob(fl.cs),
                    xmin=log(2.6e-4),xmax=log(3e-3),
                    ymin=0.35,ymax=0.5)


# Georgia ---------------------------------------------------------------
GA[date<as.Date('2020-04-08'),intervention:='none']
GA[date>=as.Date('2020-04-08') & date<as.Date('2020-04-24'),intervention:='stay-at-home']
GA[date>=as.Date('2020-04-24'),intervention:='reopening']
GA[,intervention:=factor(intervention,levels=unique(intervention))]

ga <-  ggplot(GA,aes(deaths_pc,growth_rate,color=intervention))+
#   geom_line(data=Tots,alpha=0.2,lwd=2,aes(group=state))+
  geom_line(data=nyc,aes(color='black'),color='red',lwd=2)+
  geom_line(data=NY,aes(color='black'),color='black',lwd=2)+
  geom_line(data=SWE,aes(color='black'),color=swy,lwd=2)+
  geom_line(lwd=3)+
  scale_x_continuous(trans='log',limits=c(1e-6,3e-3),breaks=10^(-6:0))+
  scale_y_continuous(limits=c(-0.1,0.5))+
  theme_bw(base_size=15)+
  scale_color_manual(values=intervention_cols$color[c(1,2,4)])+
  geom_hline(yintercept = 0)+
  ggtitle(paste('Georgia: peaked',bgs[state=='Georgia',date]))+
  theme(legend.position = c(.5,0.04),
        legend.direction = 'horizontal',
        legend.background = element_blank(),
        legend.box.background = element_rect(color = 'black'))

ga.rt <- ggplot(GA[date>as.Date('2020-03-01')],
                aes(date,growth_rate,color=intervention))+
  geom_line(lwd=1.5)+
  theme_bw()+
  scale_color_manual(values=intervention_cols$color[c(1,2,4)])+
  theme(legend.position='none')+
  scale_y_continuous('r(t)')+  
  scale_x_continuous(limits = as.Date(c('2020-03-01','2020-09-01')),
                     breaks=as.Date(c('2020-03-01','2020-05-01','2020-07-01','2020-09-01')),
                     labels = c('Mar','May','July','Sep'))

ga.cs <- ggplot(GA[date>as.Date('2020-03-01')],
                aes(date,new_confirmed,fill=intervention))+
  geom_bar(stat='identity')+
  theme_bw()+
  scale_fill_manual(values=intervention_cols$color[c(1,2,4)])+
  theme(legend.position='none')+
  scale_y_continuous('N(t)')+
  scale_x_continuous(name=NULL,breaks=NULL)

ga2=ga+annotation_custom(ggplotGrob(ga.rt),
                         xmin=log(3e-4),xmax=log(3e-3),
                         ymin=0.17,ymax=0.35)+
  annotation_custom(ggplotGrob(ga.cs),
                    xmin=log(2.8e-4),xmax=log(3e-3),
                    ymin=0.35,ymax=0.5)


# California --------------------------------------------------------------
CA[date<as.Date('2020-03-19'),intervention:='none']
CA[date>=as.Date('2020-03-19') & date<as.Date('2020-05-07'),intervention:='stay-at-home']
CA[date>=as.Date('2020-05-07') & date<as.Date('2020-06-07'),intervention:='Stage 2']
CA[date>=as.Date('2020-06-07') & date<as.Date('2020-06-28'),intervention:='Stage 3']
CA[date>=as.Date('2020-06-28'),intervention:='reversal']
CA[,intervention:=factor(intervention,levels=unique(intervention))]

ca <-  ggplot(CA,aes(deaths_pc,growth_rate,color=intervention))+
#   geom_line(data=Tots,alpha=0.2,lwd=2,aes(group=state))+
  geom_line(data=nyc,aes(color='black'),color='red',lwd=2)+
  geom_line(data=NY,aes(color='black'),color='black',lwd=2)+
  geom_line(data=SWE,aes(color='black'),color=swy,lwd=2)+
  geom_line(lwd=3)+
  # geom_line(data=CA[intervention=='reversal'],color='red',lwd=2)+
  scale_x_continuous(trans='log',limits=c(1e-6,3e-3),breaks=10^(-6:0),labels = NULL,name=NULL)+
  scale_y_continuous(limits=c(-0.1,0.5),name=NULL,labels=NULL)+
  theme_bw(base_size=15)+
  scale_color_manual(values=intervention_cols$color[c(1:4,6)])+
  geom_hline(yintercept = 0)+
  ggtitle(paste('California: peaked',bgs[state=='California',date]))+
  theme(legend.position = c(.5,0.04),
        legend.direction = 'horizontal',
        legend.background = element_blank(),
        legend.box.background = element_rect(color = 'black'))

ca.rt <- ggplot(CA[date>as.Date('2020-03-01')],
                aes(date,growth_rate,color=intervention))+
  geom_line(lwd=1.5)+
  theme_bw()+
  scale_color_manual(values=intervention_cols$color[c(1:4,6)])+
  theme(legend.position='none')+
  geom_vline(xintercept = CA[grepl('reversal',intervention),min(date)])+
  scale_y_continuous('r(t)')+  
  scale_x_continuous(limits = as.Date(c('2020-03-01','2020-09-01')),
                     breaks=as.Date(c('2020-03-01','2020-05-01','2020-07-01','2020-09-01')),
                     labels = c('Mar','May','July','Sep'))

ca.cs <- ggplot(CA[date>as.Date('2020-03-01')],
                aes(date,new_confirmed,fill=intervention))+
  geom_bar(stat='identity')+
  theme_bw()+
  scale_fill_manual(values=intervention_cols$color[c(1:4,6)])+
  theme(legend.position='none')+
  scale_y_continuous('N(t)')+
  scale_x_continuous(name=NULL,breaks=NULL)

ca2=ca+annotation_custom(ggplotGrob(ca.rt),
                         xmin=log(3e-4),xmax=log(3e-3),
                         ymin=0.17,ymax=0.35)+
  annotation_custom(ggplotGrob(ca.cs),
                    xmin=log(2.65e-4),xmax=log(3e-3),
                    ymin=0.35,ymax=0.5)

# all ---------------------------------------------------------------------



ggarrange(az2,ca2,fl2,tx2,ga2,la2,nrow=3,ncol=2,align='v')
ggsave('figures/state_intervention_trajectories.png',height=18,width=16)



# plotting all states -----------------------------------------------------

states <- setdiff(unique(USA$state),c('Guam','District of Columbia'))
ny <- NULL
sw <- NULL
dum.nyc <- NULL
for (st in states){
  dum <- NY[,c('deaths_pc','growth_rate')]
  dum[,state:=st]
  ny <- rbind(ny,dum)
  
  dum <- SWE[,c('deaths_pc','growth_rate')]
  dum[,state:=st]
  sw <- rbind(sw,dum)
  
  dum <- nyc[,c('deaths_pc','growth_rate')]
  dum[,state:=st]
  dum.nyc <- rbind(dum.nyc,dum)
}

USA[state %in% states] %>%
  ggplot(aes(deaths_pc,growth_rate))+
  geom_line(data=dum.nyc,aes(color='black'),color='red',lwd=1.2)+
  geom_line(data=ny,aes(color='black'),color='black',lwd=1.2)+
  geom_line(data=sw,aes(color='black'),color=swy,lwd=1.2)+
  geom_line(col='steelblue',lwd=1.2)+
  scale_x_continuous(trans='log',limits=c(1e-5,3e-3),breaks=10^(-5:0))+
  scale_y_continuous(limits=c(-0.1,0.5))+
  theme_bw(base_size=15)+
  geom_hline(yintercept = 0)+
  facet_wrap(.~state)

ggsave('figures/all_states_second_wave_rD.png',height=10,width=15,units='in')
