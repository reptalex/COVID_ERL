library(magrittr)
library(data.table)
library(COVID19)
library(ggpubr)
plot_erl <- function(admin_level_1='India',admin_level_2=NULL,admin_level_3=NULL,X=NULL,
                      risk_rel_nyc=1/3,case_death_lag=11,line_col=rgb(0.1,0.8,0.2),
                      min_date=as.Date('2020-01-01'),max_date=NULL,
                      highlight_start=NULL,highlight_end=Inf,return_subplots=FALSE){
  if (is.null(max_date)){
    max_date <- Sys.Date()
  }
  nyc <- read.csv('data/nbss_nyc.csv') %>% as.data.table()
  nyc[,deaths_pc:=shift(deaths_pc,11,type='lead')]
  nyc <- nyc[date<as.Date('2020-05-01')]
  swe <- read.csv('data/nbss_sweden.csv') %>% as.data.table()
  swe[,deaths_pc:=shift(deaths_pc,11,type='lead')]
  source('scripts/utils.R')
  disps <- read.csv('data/precomputed_dispersion.csv')
  
  if (is.null(X)){
    if (!is.null(admin_level_3)){
      X <- COVID19::covid19(admin_level_1,level=3) %>% as.data.table
      if (nrow(X)==0){
        stop('COVID19 package does not have admin level 3 for input admin_level_1')
      }
      if (length(unique(X$administrative_area_level_2))>1 & is.null(admin_level_2)){
        stop('Must specify admin_level_2 due to non-uniqueness of admin_level_3')
      } else if (!is.null(admin_level_2)){
        X <- X[administrative_area_level_2==admin_level_2 & administrative_area_level_3==admin_level_3]
      } else {
        X <- X[administrative_area_level_3==admin_level_3]
      }
    } else if (!is.null(admin_level_2)){
      X <- COVID19::covid19(admin_level_1,level=2) %>% as.data.table
      if (nrow(X)==0){
        stop('COVID19 package does not have admin level 3 for input admin_level_1')
      }
      X <- X[administrative_area_level_2==admin_level_2]
    } else {
      X <- COVID19::covid19(admin_level_1) %>% as.data.table
    }
  }
  
  X$new_confirmed <- dfs(X$confirmed)
  X$new_deaths <- dfs(X$deaths)
  X <- covid19_nbss(as.data.frame(X),precomputed_dispersions = disps) %>% as.data.table
  X[,deaths_pc:=shift(deaths/population,case_death_lag,type='lead')]
  X <- X[date>=min_date & date<=max_date]
  g_cases <- ggplot(X,aes(date,new_confirmed))+
    geom_bar(stat='identity',col=line_col,fill=line_col)+
    geom_bar(stat='identity',aes(y=new_deaths),col='black',fill='black')+
    scale_y_continuous(trans='log',breaks=10^(0:6))+
    theme_bw(base_size=15)+
    ggtitle("Cases (colored) and Deaths (black)")
  g_erl <- ggplot(X,aes(deaths_pc,growth_rate))+
    geom_hline(yintercept = 0)+
    geom_line(data=nyc,color='red',lwd=2)+
    geom_line(data=swe,color=rgb(1,205/255,0),lwd=2)+
    geom_line(lwd=2,col=line_col)+
    geom_line(col=line_col,lwd=2)+
    geom_line()+
    theme_bw(base_size=15)+
    scale_x_continuous(trans='log',limits=c(1e-7,4e-3),breaks=10^(-7:-3))
  
  ttl <- admin_level_1
  if (!is.null(admin_level_2)){
    ttl <- paste(ttl,admin_level_2,sep=',')
  }
  if (!is.null(admin_level_3)){
    ttl <- paste(ttl,admin_level_3,sep=',')
  }
    g_erl <- g_erl+ggtitle(paste(ttl,'ERL'))
    
  if (risk_rel_nyc!=1){
    g_erl <- g_erl+geom_line(data=nyc,aes(x=deaths_pc*risk_rel_nyc),col='red')
  }
  
  g_combined <- ggarrange(g_cases,g_erl,nrow=2,align='v')
  
  if(return_subplots){
    return(list('cases'=g_cases,'erl'=g_erl,'combined'=g_combined))
  } else{
    return(g_combined)
  }
}
