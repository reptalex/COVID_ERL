library(data.table)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(cdcfluview)

ILI <- cdcfluview::ilinet(region = 'state') %>% as.data.table

ILI <- ILI[!grepl('Mariana',region)][!region=='Florida']
mn <- ILI[year<2020,list(unweighted_ili=mean(unweighted_ili),
                         year='baseline'),by=c('region','week')]
COVID <- rbind(ILI[year %in% c(2019,2020,2021),c('unweighted_ili','region','week','year')],mn)
COVID[,year:=factor(year)]
sts <- c('New York','Illinois','')
ggplot(ILI,aes(week,unweighted_ili,group=year))+
  geom_line(alpha=0.3)+
  geom_line(data=mn,lwd=2)+
  geom_line(data=COVID,aes(color=year),lwd=1.5)+
  theme_bw()+
  facet_wrap(.~region,scales = 'free_y')+
  # scale_y_continuous(trans='log')+
  ggtitle('Influenza-like Illness')+
  scale_color_manual(values=c('red','deepskyblue2','orange','black'))