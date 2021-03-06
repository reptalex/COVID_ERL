library(data.table)
library(deSolve)
library(magrittr)
library(zoo)
library(KFAS)
library(parallel)
library(tidyverse)
library(tsoutliers)
library(progress)
library(lubridate)
library(EpiEstim)
library(mgcv)

dfs <- function(x){
  x <- c(x[1],diff(x))
  x[x<0] <- NA
  return(x)
}
# 
# get_cori <- function(df.in, 
#                      icol_name, 
#                      out_name = 'Re_Cori',
#                      window = 1, 
#                      GI_mean=4, 
#                      GI_var=4.75^2,
#                      wend = TRUE){
#   
#   max.obs.time <- df.in %>% filter(!is.na(!!sym(icol_name))) %>% pull(time) %>% tail(1)
#   
#   
#   idat <- df.in %>%
#     #filter(get(icol_name) > 0 & !is.na(get(icol_name))) %>%
#     complete(time = 2:max.obs.time) %>%
#     arrange(time) %>%
#     filter(time <= max.obs.time)
#   
#   nas <- which(is.na(getElement(df.in,icol_name)))
#   df.in[nas,(icol_name):=0]
#   # idat[icol_name] <- na_to_0(idat[icol_name])
#   #mutate(cleaned = ifelse(is.na(!!sym(icol_name)) & time <= max.obs.time, 0, !!sym(icol_name)))
#   
#   
#   ts <- idat$time
#   ts <- ts[ts > 1 & ts <= (max(ts)-window+1)]
#   te <- ts+(window-1)
#   
#   estimate_R(
#     incid = pull(idat, !!icol_name),
#     method = "uncertain_si",
#     config = make_config(
#       list(
#         mean_si = GI_mean,
#         min_mean_si = GI_mean -1,
#         max_mean_si = GI_mean + 1,
#         std_mean_si = 1.5,
#         std_std_si = 1.5,
#         std_si = sqrt(GI_var),
#         min_std_si = sqrt(GI_var)*.8,
#         max_std_si = sqrt(GI_var)*1.2,
#         n1 = 50,
#         n2 = 100, 
#         t_start=ts,
#         t_end=te
#       )
#     )
#   ) -> outs
#   
#   R <- outs$R %>%
#     mutate(time = if(wend == TRUE) t_end else ceiling((t_end+t_start)/2) ) %>%
#     select(time, `Mean(R)`, `Quantile.0.025(R)`, `Quantile.0.975(R)`) %>%
#     setNames(c('time', paste0(out_name, '.mean'), paste0(out_name, '.025'), paste0(out_name, '.975'))) %>%
#     as.data.table
#   
#   R[,Cori.smooth:=frollapply(Re_Cori.mean,7,mean,align = 'right')]
#   return(R)
# }

# Requires to get sockets parallization to work on os x. 
library(rstudioapi)
if (Sys.getenv("RSTUDIO") == "1" && !nzchar(Sys.getenv("RSTUDIO_TERM")) && 
    Sys.info()["sysname"] == "Darwin" && getRversion() >= "4.0.0") {
  if(versionInfo()$version < "1.3.1056"){
    parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
  }  
}

nbs <- function(x,name='growth_rate',remove_initial_zeros=TRUE,...){
  y <- rep(NA,length(x))
  first_nz <- min(which(x>0))
  if (name=='all'){
    mdl <- nbss(x[first_nz:length(x)],...)
    
    dd <- matrix(NA,nrow=first_nz-1,ncol=ncol(mdl))
    colnames(dd) <- colnames(mdl)
    y <- rbind(dd,mdl)
  } else {
    y[first_nz:length(x)] <- nbss(x[first_nz:length(x)],...)  %>% getElement(name)
  }
  return(y)
}

nbss <- function(x,remove_outliers=TRUE,filtering=FALSE,dispersion=NULL){
  nb_model <- function(x, pars){
    model_nb <- SSModel(x ~ SSMtrend(2, Q=list(0, NA),
                                     P1=diag(c(10, 1)),
                                     a1=c(0, 0),
                                     state_names=c("level", "trend"))+
                          SSMseasonal(7),
                        u=rep(exp(pars[1]), length(x)), distribution="negative binomial")
    fit <- fitSSM(model_nb, c(0), method="L-BFGS-B", control=list(maxit=200))
    return(fit)
  }
  
  if (remove_outliers==TRUE){
    x <- outlier_detection(x)
  }
  if (is.null(dispersion)){
    logLik_nb <- function(x, pars){
      fit <- nb_model(x, pars)
      ll <- logLik(fit$model, marginal = TRUE)
      return(-ll)
    }
    
    res <- tryCatch(optim(c(-1), function(y) logLik_nb(x, y), method="Brent", lower=-2, upper=2) , error=function(e) NULL)
    if (is.null(res)) return(NULL)
    fit <- nb_model(x, res$par)
  } else {
    fit <- tryCatch(nb_model(x,dispersion),error=function(e) NULL)
    res <- NULL
    res$par[1] <- dispersion
  }
  if (filtering==TRUE){
    sm_signal <- KFS(fit$model, filtering="signal",smoothing='none')
    sm_state <- KFS(fit$model, filtering="state",smoothing='none')
    
    out <- data.frame(p2.5_position = c(qnorm(0.025, sm_state$a[-1,'level'], (sqrt(sm_state$P[1,1,-1])))), 
                      p97.5_position = c(qnorm(0.975, sm_state$a[-1,'level'],(sqrt(sm_state$P[1,1,-1])))), 
                      mean_position = (c(sm_state$a[-1,'level'])),
                      p2.5_growth_rate = c(qnorm(0.025, sm_state$a[-1,'trend'], (sqrt(sm_state$P[2,2,-1])))), 
                      p97.5_growth_rate = c(qnorm(0.975, sm_state$a[-1,'trend'], (sqrt(sm_state$P[2,2,-1])))),
                      p25_growth_rate = c(qnorm(0.25, sm_state$a[-1,'trend'], (sqrt(sm_state$P[2,2,-1])))), 
                      p75_growth_rate = c(qnorm(0.75, sm_state$a[-1,'trend'], (sqrt(sm_state$P[2,2,-1])))), 
                      growth_rate = (c(sm_state$a[-1,'trend'])),
                      percentile_0_growth_rate =c(pnorm(0, sm_state$a[-1,'trend'], (sqrt(sm_state$P[2,2,-1])))),
                      growth_rate = (c(sm_state$a[-1,'trend'])),
                      dispersion=res$par[1], 
                      z_score_growth_rate = c(sm_state$a[-1,2]/sqrt(sm_state$P[2,2,-1])))
  } else {
    sm_signal <- KFS(fit$model, smoothing="signal")
    sm_state <- KFS(fit$model, smoothing="state")
    out <- data.frame(p2.5_signal = exp(qnorm(0.025, sm_signal$thetahat, sqrt(c(sm_signal$V_theta)))), 
                      p97.5_signal = exp(qnorm(0.975, sm_signal$thetahat, sqrt(c(sm_signal$V_theta)))), 
                      mean_signal = exp(sm_signal$thetahat), 
                      p2.5_position = c(qnorm(0.025, sm_state$alphahat[,1], (sqrt(sm_state$V[1,1,])))), 
                      p97.5_position = c(qnorm(0.975, sm_state$alphahat[,1],(sqrt(sm_state$V[1,1,])))), 
                      mean_position = (c(sm_state$alphahat[,1])),
                      p2.5_growth_rate = c(qnorm(0.025, sm_state$alphahat[,2], (sqrt(sm_state$V[2,2,])))), 
                      p97.5_growth_rate = c(qnorm(0.975, sm_state$alphahat[,2], (sqrt(sm_state$V[2,2,])))),
                      p25_growth_rate = c(qnorm(0.25, sm_state$alphahat[,2], (sqrt(sm_state$V[2,2,])))), 
                      p75_growth_rate = c(qnorm(0.75, sm_state$alphahat[,2], (sqrt(sm_state$V[2,2,])))), 
                      growth_rate = (c(sm_state$alphahat[,2])),
                      percentile_0_growth_rate =c(pnorm(0, sm_state$alphahat[,2], (sqrt(sm_state$V[2,2,])))),
                      growth_rate = (c(sm_state$alphahat[,2])),
                      dispersion=res$par[1], 
                      z_score_growth_rate = c(sm_state$alphahat[,2]/sqrt(sm_state$V[2,2,])))
  }
  return(out)
}

# Main Model function -- CALL THIS 
# dat should be a data.frame for a single region to be modeled 
#    with column "new_confirmed" which is count data
#    arranged by date (increasing)
#
# return_fit=TRUE is for debugging mostly. 
#
# any other columns in dat are included in output of this function. 
# 
# 
# in addition the following key variables are included as model output
#   growth_rate == exponential growth rate
#   pX_Y is quantile X for the quantitly Y
#   mean_ is the mean of the quantity Y
#   z_score_growth_rate = growth_rate/sd_growth_rate
fit_covid_ssm <- function(d, series="new_confirmed", precomputed_dispersions=NULL,
                          return_fit=FALSE,maxiter=200,filtering=FALSE){
  dat <- d
  
  # Helper functions
  nb_model <- function(dat, pars){
    model_nb <- SSModel(dat[,series] ~ SSMtrend(2, Q=list(0, NA),
                                                     P1=diag(c(10, 1)),
                                                     a1=c(0, 0),
                                                     state_names=c("level", "trend"))+
                          SSMseasonal(7),
                        u=rep(exp(pars[1]), nrow(dat)), distribution="negative binomial")
    fit <- fitSSM(model_nb, c(0), method="L-BFGS-B", control=list(maxit=200))
    return(fit)
  }
  logLik_nb <- function(dat, pars){
    fit <- nb_model(dat, pars)
    ll <- logLik(fit$model, marginal = TRUE)
    return(-ll)
  }
  
  # Remove preceding zeros
  dat <- dat %>% 
    arrange(date)
  dat$cs = cumsum(ifelse(is.na(dat[,series]), 0, dat[,series]))
  pass <- dat$cs > 0 # rows to keep in analysis 
  # dat <- dat %>% filter(cs > 0)
  
  
  if (nrow(dat[pass,]) < 10) return(cbind(d, error="not enough non-zero"))
  
  
  # Remove weekend zeros
  dat[[series]] <- ifelse((dat[[series]] == 0) & (format(dat$date, "%u") %in% c(6,7)), 
                          NA, dat[[series]])
  # Remove holiday zeros
  dat[[series]] <- ifelse((dat[[series]] == 0) & (dat$date %in% c(ymd("2020-11-26"), # Thanksgiving
                                                                  ymd("2020-12-25"), 
                                                                  ymd("2020-12-31"), 
                                                                  ymd("2021-01-01"))),
                          NA, dat[[series]])
  
  
  # outlier detection for early outbreak
  pass <- custom_processors(dat, pass)
  if (sum(dat[pass,series]!=0, na.rm=TRUE) < 10) return(cbind(d, error="not enough non-zero"))
  tryCatch({
    tmp <-getElement(dat,series)[pass]
    #tmp <- ifelse(is.na(tmp), 0, tmp)
    tmp <- na.approx(tmp)
    outlier_filtered_ts <- tmp %>% outlier_detection
    outlier_filtered_ts[is.na(getElement(dat, series)[pass])] <- NA
    tmp <- rep(NA, nrow(dat))
    tmp[pass] <- outlier_filtered_ts
    dat <- mutate(dat, series=tmp)    
  },  error = function(err){
    return(cbind(d, error="outlier detection errored"))
  })
  
  if (length(dat[pass,series][!is.na(dat[pass,series])]) < 10) return(cbind(d, error="too many NA"))
  
  
  if (!is.null(precomputed_dispersions)){
    res <- list()
    res$par <- dplyr::filter(precomputed_dispersions, id==unique(dat$id))$dispersion
    if (length(res$par)==0) return(cbind(d, error="could not find precomputed dispersions")) # no precomputed dispersion (previously was not able to fit likely)
  } else {
    res <- tryCatch(optim(c(-1), function(x) logLik_nb(dat[pass,], x), method="Brent", lower=-2, upper=3), 
                    error = function(e) NULL)
    if (is.null(res)) return(cbind(d, error="dispersion optimization failed"))
    if (res$convergence != 0) return(cbind(d, error="dispersion optimization failed"))
  }
  
  # now fit the model with the optimized dispersion parameters
  fit <- nb_model(dat[pass,], res$par)
  if(return_fit) return(fit)
  if (fit$optim.out$convergence != 0) return(cbind(d, error="model optimiztion (not dispersion) failed"))
  if (filtering==FALSE){
    sm_signal <- KFS(fit$model, smoothing="signal")
    sm_state <- KFS(fit$model, smoothing="state")
    
    if (series=="new_confirmed"){
      out <- data.frame(p2.5_signal = exp(qnorm(0.025, sm_signal$thetahat, sqrt(c(sm_signal$V_theta)))), 
                        p97.5_signal = exp(qnorm(0.975, sm_signal$thetahat, sqrt(c(sm_signal$V_theta)))), 
                        mean_signal = exp(sm_signal$thetahat), 
                        p2.5_position = c(qnorm(0.025, sm_state$alphahat[,1], (sqrt(sm_state$V[1,1,])))), 
                        p97.5_position = c(qnorm(0.975, sm_state$alphahat[,1],(sqrt(sm_state$V[1,1,])))), 
                        mean_position = (c(sm_state$alphahat[,1])),
                        p2.5_growth_rate = c(qnorm(0.025, sm_state$alphahat[,2], (sqrt(sm_state$V[2,2,])))), 
                        p97.5_growth_rate = c(qnorm(0.975, sm_state$alphahat[,2], (sqrt(sm_state$V[2,2,])))),
                        p25_growth_rate = c(qnorm(0.25, sm_state$alphahat[,2], (sqrt(sm_state$V[2,2,])))), 
                        p75_growth_rate = c(qnorm(0.75, sm_state$alphahat[,2], (sqrt(sm_state$V[2,2,])))), 
                        growth_rate = (c(sm_state$alphahat[,2])),
                        percentile_0_growth_rate =c(pnorm(0, sm_state$alphahat[,2], (sqrt(sm_state$V[2,2,])))),
                        dispersion=res$par[1], 
                        z_score_growth_rate = c(sm_state$alphahat[,2]/sqrt(sm_state$V[2,2,])))
    } else if (series == "new_deaths"){
      out <- data.frame(p2.5_signal_deaths = exp(qnorm(0.025, sm_signal$thetahat, sqrt(c(sm_signal$V_theta)))), 
                        p97.5_signal_deaths = exp(qnorm(0.975, sm_signal$thetahat, sqrt(c(sm_signal$V_theta)))), 
                        mean_signal_deaths = exp(sm_signal$thetahat), 
                        p2.5_position_deaths = c(qnorm(0.025, sm_state$alphahat[,1], (sqrt(sm_state$V[1,1,])))), 
                        p97.5_position_deaths = c(qnorm(0.975, sm_state$alphahat[,1],(sqrt(sm_state$V[1,1,])))), 
                        mean_position_deaths = (c(sm_state$alphahat[,1])),
                        p2.5_growth_rate_deaths = c(qnorm(0.025, sm_state$alphahat[,2], (sqrt(sm_state$V[2,2,])))), 
                        p97.5_growth_rate_deaths = c(qnorm(0.975, sm_state$alphahat[,2], (sqrt(sm_state$V[2,2,])))),
                        p25_growth_rate_deaths = c(qnorm(0.25, sm_state$alphahat[,2], (sqrt(sm_state$V[2,2,])))), 
                        p75_growth_rate_deaths = c(qnorm(0.75, sm_state$alphahat[,2], (sqrt(sm_state$V[2,2,])))), 
                        growth_rate_deaths = (c(sm_state$alphahat[,2])),
                        percentile_0_growth_rate_deaths =c(pnorm(0, sm_state$alphahat[,2], (sqrt(sm_state$V[2,2,])))),
                        dispersion_deaths=res$par[1], 
                        z_score_growth_rate_deaths = c(sm_state$alphahat[,2]/sqrt(sm_state$V[2,2,])))
    }
  } else {
    if (series=="new_confirmed"){
      sm_signal <- KFS(fit$model, filtering="signal",smoothing='none')
      sm_state <- KFS(fit$model, filtering="state",smoothing='none')
      
      out <- data.frame(p2.5_position = c(qnorm(0.025, sm_state$a[-1,'level'], (sqrt(sm_state$P[1,1,-1])))), 
                        p97.5_position = c(qnorm(0.975, sm_state$a[-1,'level'],(sqrt(sm_state$P[1,1,-1])))), 
                        mean_position = (c(sm_state$a[-1,'level'])),
                        p2.5_growth_rate = c(qnorm(0.025, sm_state$a[-1,'trend'], (sqrt(sm_state$P[2,2,-1])))), 
                        p97.5_growth_rate = c(qnorm(0.975, sm_state$a[-1,'trend'], (sqrt(sm_state$P[2,2,-1])))),
                        p25_growth_rate = c(qnorm(0.25, sm_state$a[-1,'trend'], (sqrt(sm_state$P[2,2,-1])))), 
                        p75_growth_rate = c(qnorm(0.75, sm_state$a[-1,'trend'], (sqrt(sm_state$P[2,2,-1])))), 
                        growth_rate = (c(sm_state$a[-1,'trend'])),
                        percentile_0_growth_rate =c(pnorm(0, sm_state$a[-1,'trend'], (sqrt(sm_state$P[2,2,-1])))),
                        growth_rate = (c(sm_state$a[-1,'trend'])),
                        dispersion=res$par[1], 
                        z_score_growth_rate = c(sm_state$a[-1,2]/sqrt(sm_state$P[2,2,-1])))
    } else if (series == "new_deaths"){
      
      out <- data.frame(p2.5_position_deaths = c(qnorm(0.025, sm_state$a[-1,'level'], (sqrt(sm_state$P[1,1,-1])))), 
                        p97.5_position_deaths = c(qnorm(0.975, sm_state$a[-1,'level'],(sqrt(sm_state$P[1,1,-1])))), 
                        mean_position_deaths = (c(sm_state$a[-1,'level'])),
                        p2.5_growth_rate_deaths = c(qnorm(0.025, sm_state$a[-1,'trend'], (sqrt(sm_state$P[2,2,-1])))), 
                        p97.5_growth_rate_deaths = c(qnorm(0.975, sm_state$a[-1,'trend'], (sqrt(sm_state$P[2,2,-1])))),
                        p25_growth_rate_deaths = c(qnorm(0.25, sm_state$a[-1,'trend'], (sqrt(sm_state$P[2,2,-1])))), 
                        p75_growth_rate_deaths = c(qnorm(0.75, sm_state$a[-1,'trend'], (sqrt(sm_state$P[2,2,-1])))), 
                        growth_rate_deaths = (c(sm_state$a[-1,'trend'])),
                        percentile_0_growth_rate_deaths =c(pnorm(0, sm_state$a[-1,'trend'], (sqrt(sm_state$P[2,2,-1])))),
                        growth_rate_deaths = (c(sm_state$a[-1,'trend'])),
                        dispersion_deaths=res$par[1], 
                        z_score_growth_rate_deaths = c(sm_state$a[-1,2]/sqrt(sm_state$P[2,2,-1])))
    }
  }
  
  if (any(grepl('administrative_area',colnames(dat)))){
    if (any(out$p97.5_signal > 1e6)) return(cbind(d, error="97.5_signal > 1e6"))
  }
  if (quantile(abs(out$z_score_growth_rate), probs=0.75) < 0.4) return(cbind(d, error="quantile(abs(out$z_score_growth_rate), probs=0.75) < 0.4"))
  
  tack_on <- matrix(NA, nrow=nrow(dat), ncol=ncol(out)) %>% as.data.frame()
  colnames(tack_on) <- colnames(out)
  tack_on[pass,] <- out
  return(cbind(dat, tack_on))
}


covid19_nbss <- function(dat,series="new_confirmed", level='all',
                         mc.cores=1, precomputed_dispersions=NULL,
                         filtering=FALSE){
  if (level=="country"){
    tmp <- filter(dat, administrative_area_level==1)
  } else if (level=="state"){
    tmp <- filter(dat, administrative_area_level==2)
  } else if (level=="all") {
    tmp <- dat
  } else {
    stop("only level variables that are supported are all, country, and state")
  }
  tmp <- dat %>% 
    mutate(new_deaths = ifelse(new_deaths < 0, 0, new_deaths), 
           new_confirmed = ifelse(new_confirmed < 0, 0, new_confirmed)) %>% 
    as.data.frame() %>% 
    split(.$id)
  if (mc.cores == 1){
    fits <- list()
    pb <- progress_bar$new(total = length(tmp), format=" [:bar] :percent eta: :eta")
    for (i in 1:length(tmp)){
      pb$tick()
      fits[[i]] <- fit_covid_ssm(tmp[[i]], series, precomputed_dispersions,filtering=filtering)
      #if (is.null(fits[[i]])) stop("foo")
    }  
  } else {
    cl <- parallel::makeCluster(mc.cores)
    parallel::clusterEvalQ(cl, {
      library(tidyverse)
      library(lubridate)
      library(KFAS)
      library(tsoutliers)
      library(data.table)
    })
    parallel::clusterExport(cl,c("custom_processors", "outlier_detection", "fit_covid_ssm"))
    fits <- parLapply(cl, tmp,  function(x,series,precomputed_dispersions,filtering) fit_covid_ssm(x, series, precomputed_dispersions,filtering=filtering),
                      series=series,precomputed_dispersions=precomputed_dispersions,filtering=filtering)
    stopCluster(cl)
    rm('cl')
  }
  fits <- bind_rows(fits)
  return(fits)
}


custom_processors <- function(dat,pass){
  # Iowa and Indiana
  # if (unique(dat$administrative_area_level_2)%in%c("Iowa", "Indiana", "Kentucky")){
  
  if (unique(dat$administrative_area_level_1)=="United States"){
    if (!is.na(unique(dat$administrative_area_level_2))){
      if (unique(dat$administrative_area_level_2 != "Washington")){
        pass <- pass & (dat$date > ymd("2020-02-28") )
      }
    }
    # Insert new US custom processors here
  }
  return(pass)
}


outlier_detection <- function(x){
  if (sum(x, na.rm=TRUE) < 100 &  sum(x==0, na.rm=TRUE) > length(x)*0.5){
    return(x)
  }
  tryCatch({
    res <- tso(ts(x),
               type="TC", delta=0.1, maxit.iloop = 100, maxit.oloop = 10,  
               #tsmethod = "auto.arima", args.tsmethod = list(allowdrift = FALSE, ic = "bic", stationary=TRUE),
               tsmethod="arima", args.tsmethod=list(order=c(1,1,2), method="ML", transform.pars=TRUE, optim.method="BFGS"),
               cval=4)
  }, error = function(err){
    res <- tso(ts(x),
               type="TC", delta=0.1, maxit.iloop = 100, maxit.oloop = 10, 
               #tsmethod = "auto.arima", args.tsmethod = list(allowdrift = FALSE, ic = "bic", stationary=TRUE),
               tsmethod="arima", args.tsmethod=list(order=c(1,1,2), method="ML", transform.pars=FALSE),
               cval=4)
  })
  res <- res$outliers %>% 
    filter(tstat > 10)
  x[res$ind] <- NA
  return(x)
}


seird_solve_intervention_relaxation <- function(t,state,parameters){
  with(as.list(c(t,state,parameters)),{
    if (t>lag_onset_to_death){
      Dlag=lagvalue(t-lag_onset_to_death)[5]
    } else {
      Dlag=0
    }
    dS=lambda-(as.numeric(Dlag<=D_intervention)*beta_pre+
                 as.numeric(Dlag>D_intervention & Dlag<=D_relaxation)*(beta_intervention)+
                 as.numeric(Dlag>D_relaxation)*beta_relaxation)*S*I-m*S
    dE=(as.numeric(Dlag<=D_intervention)*beta_pre+
          as.numeric(Dlag>D_intervention & Dlag<=D_relaxation)*(beta_intervention)+
          as.numeric(Dlag>D_relaxation)*beta_relaxation)*S*I-a*E-m*E
    dI=a*E-m_inf*I-gamma*I
    dR=gamma*I-m*R
    dD=m_inf*I
    
    list(c(dS, dE, dI, dR, dD))
  })
}

seird <- function(r,cfr=0.006,S0=3.27e8,start_date=as.Date('2020-01-15'),days=200,
                  lag_onset_to_case=7,
                  lag_onset_to_death=18,window_size=21,half_life=NULL,
                  intervention_efficacy=0,relaxation=1,intervention_deaths=Inf,relaxation_deaths=Inf,
                  gamma=1/9,a=1/3,gr_estimation='nbss',growth_rate_deaths=FALSE,
                  day_of_week_effects=data.table('day_of_week'=c('Saturday','Sunday',
                                                                 'Monday','Tuesday',
                                                                 'Wednesday','Thursday','Friday'),
                                                 'const'=c(0.6,0.4,1,1,1,1,.9)),
                  case_detection=1e-4,nb_size=5){
  
  state <- c('S'=S0,'E'=0,'I'=1,'R'=0,'D'=0)
  times <- seq(0, days, by = 0.01)
  m=8.685/100000/365
  
  m_inf <- cfr*gamma/(1-cfr)
  alpha <- m_inf/m
  
  # alpha=1.1
  c=a/(r+gamma+m_inf)   # the ratio of E/I during exponential phase
  beta_pre=(r+a+m)/(c*S0)
  beta_intervention=beta_pre*(1-intervention_efficacy)
  beta_relaxation=beta_pre*relaxation
  lambda=11.8/1000/365
  parameters <- c('m'=m,
                  'lambda'=lambda,
                  'gamma'=gamma,
                  'a'=a,
                  'beta_pre'=beta_pre,
                  'beta_intervention'=beta_intervention,
                  'beta_relaxation'=beta_relaxation,
                  'm_inf'=m_inf,
                  'D_intervention'=intervention_deaths,
                  'D_relaxation'=relaxation_deaths,
                  'r'=r,
                  'lag_onset_to_death'=lag_onset_to_death)
  
  out <- dede(y = state, times = times, func = seird_solve_intervention_relaxation, parms = parameters) %>% as.data.table
  out[,day:=ceiling(time)]
  
  out <- out[,list(S=S[.N],
                   E=E[.N],
                   I=I[.N],
                   R=R[.N],
                   D=D[.N]),by=day]
  out[,date:=seq(start_date,start_date+days,by='day')]
  out[,day_of_week:=weekdays(date)]
  setkey(out,day_of_week)
  setkey(day_of_week_effects,day_of_week)
  
  out <- day_of_week_effects[out]
  setkey(out,date)
  
  
  out[,rt:=a*E/I-m_inf-gamma]
  out[,D:=shift(D,lag_onset_to_death)]
  out[,beta:=(as.numeric(D<=intervention_deaths)*beta_pre+
                as.numeric(D>intervention_deaths & D<relaxation_deaths)*(beta_intervention)+
                as.numeric(D>=relaxation_deaths)*beta_relaxation)]
  out[,new_infections:=beta*S*I]
  
  out[is.na(D),D:=0]
  out[,new_deaths:=c(0,diff(D))]
  # out[,sample_deaths:=rpois(.N,lambda = new_deaths*const)]
  out[,new_confirmed:=shift(rnbinom(.N,mu=case_detection*I*const,
                              size=nb_size),lag_onset_to_case)]
  # out[cumsum(sample_deaths)==0,sample_deaths:=NA]
  out[cumsum(new_confirmed)==0,new_confirmed:=NA]
  out$n <- NULL
  out[,r:=r]
  out[,cfr:=cfr]
  out[,lag_onset_to_case:=lag_onset_to_case]
  out[,lag_onset_to_death:=lag_onset_to_death]
  out[,intervention_efficacy:=intervention_efficacy]
  out[,intervention_deaths:=intervention_deaths]
  out[,relaxation:=relaxation]
  out[,relaxation_deaths:=relaxation_deaths]
  out[,deaths_pc:=D/S0]
  
  return(out)
}

plot_seird <- function(x,sep=FALSE){
  dd <- data.table('N'=c(x$S,x$E,x$I,x$R,x$D))
  dd$day <- rep(x$day,times=5)
  dd$compartment <- rep(c('S','E','I','R','D'),each=nrow(x))
  
  g1 <- ggplot(dd[compartment != "D"],aes(day,N,color=compartment))+
    geom_line(lwd=2)+
    theme_bw(base_size=12)+
    theme(legend.position=c(0.8,0.5))
  g2 <- ggplot(dd[compartment=='D'],aes(day,N))+
    geom_line(lwd=2)+
    theme_bw(base_size=12)
  g3 <- ggplot(x,aes(day,rt))+
    geom_line(col='green',lwd=2)+
    scale_y_continuous('r(t)')+
    theme_bw(base_size=12)
  if (sep){
    return(list('SEIR'=g1,'D'=g2,'r'=g3))
  } else {
    return(ggarrange(g1,g2,g3,nrow=3,align='v'))
  }
}
