library(viridis)
source('scripts/utils.R')

# seird parameters -------------------------------------------------------------------
# S0=8.4e6       #NYC population size
r <- log(2)/3  #exponential growth rate of cases - 3 day doubling time.
gamma <- 1/9   #infectious period
a <- 1/3       #incubation period
ifr <- 0.004   #infection fatality rate
mu <- ifr*gamma/(1-ifr)   #mortality rate of infected individuals
m <- 0         #mortality rate of uninfected individuals - approx 0 for 400-day sims
S0=19e6        #NY state poopulation size
c=a/(r+gamma+mu)   # the ratio of I/E during exponential phase
beta=(r+a+m)/(c*S0)   # transmission rate for SEIR model

# Unmitigated epidemic simulation -----------------------------------------

x <- seird(r,ifr,days=400,S0 = S0)
D <- x$D
Z <- x[,S0-S]
rt <- x$rt
lag_onset_to_death <- unique(x$lag_onset_to_death) ### this was a reporting lag in the simulation above


### below: the nonlinear relationship between deaths and cumulative incidence
plot(shift(D,lag_onset_to_death,type='lead'),log(S0/(S0-Z))*mu/beta,pch=16,lwd=2)
abline(0,1)


# Analytical estimates of ERLs --------------------------------------------
lD <- (1/(gamma+mu)+1/a)  ### lag for per-capita deaths to estimate cumulative incidence
tau <- round(lD)+lag_onset_to_death

### below is the quadratic analytical estimate for the ERL obtained in the form D(r) instead of r(D).
# r_vals <- seq(-gamma-mu,beta*S0/(gamma+mu),length.out=100)
# Quadratic_erl_estimate <- mu*S0/(gamma+mu)-mu/beta*(1+lD*r_vals+r_vals^2/a/(gamma+mu))

# Epidemics with intervenetions & relaxations -----------------------------

efficacy <- 0.9
x_stop <- seird(r,ifr,intervention_deaths=10,
                       intervention_efficacy=efficacy,S0=S0,
                       relaxation_deaths=Inf,days=400)

# plot_seird(x_stop)
x_relax <- seird(r,ifr,intervention_deaths=20,
                 intervention_efficacy=efficacy,S0=S0,
                 relaxation_deaths=x_stop[day==200,D],days=400)
# plot_seird(x_relax)
x_break_erl <-  seird(r,ifr*1.3,intervention_deaths=30,
                       intervention_efficacy=efficacy,S0=S0,
                       relaxation_deaths=1.52e4,days=400)

# plot_seird(x_break_erl)

# Epidemics over Time ----------------------------------------------------------------

cols <- c('red',viridis::viridis(3))

xx <- rbind(x[,c('day','date','new_confirmed','D','rt')],
            x_stop[,c('day','date','new_confirmed','D','rt')],
            x_relax[,c('day','date','new_confirmed','D','rt')],
            x_break_erl[,c('day','date','new_confirmed','D','rt')])

xx$scenario <- c(rep('Unmitigated',nrow(x)),
                 rep('Contained',nrow(x_stop)),
                 rep('Relaxation to ERL',nrow(x_relax)),
                 rep('ERL Rejected',nrow(x_break_erl)))
xx$scenario <- factor(xx$scenario,levels=c('Unmitigated','Contained','Relaxation to ERL','ERL Rejected'))

g_cases <- ggplot(xx[day>50 & day<270 & new_confirmed!=0],aes(date,new_confirmed,fill=scenario))+
  geom_bar(stat='identity',position='identity')+
  scale_fill_manual(values=cols)+
  facet_wrap(.~scenario,nrow=4,scales = 'free_y')+
  theme_bw(base_size=15)+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.position='none'
  )+
  scale_y_continuous('Cases, N(t)')


# ERLs --------------------------------------------------------------------

xx[,lagged_dpc:=shift(D/S0,tau,type='lead'),by=scenario]

um <- xx[scenario=='Unmitigated']
um_pk <- um[day>10 & rt>0,max(lagged_dpc,na.rm=T)]  ##unmitigated peak
ep <- xx[scenario=='Contained' & !is.na(lagged_dpc)][day==max(day)] #current state of contained epidemic
r_max <- approx(um$lagged_dpc,um$rt,xout=ep$lagged_dpc)$y
g_stop <- ggplot(xx[scenario %in% c('Unmitigated','Contained')],aes(lagged_dpc,rt))+
  geom_line(aes(color=scenario),lwd=2)+
  scale_x_continuous(name = NULL,trans='log',limits=c(1e-4,4e-3),breaks=10^(-4:-2),labels=NULL)+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.position='none'
  )+
  scale_color_manual(values=cols[1:2])+
  geom_point(data=xx[scenario=='Contained' & !is.na(lagged_dpc)][day==max(day)],cex=6,col=cols[2])+
  scale_y_continuous('r(t)')+
  geom_segment(aes(x=ep$lagged_dpc,xend=ep$lagged_dpc,
                   y=ep$rt,yend=r_max),
               arrow=arrow(),lty=2)+
  geom_vline(xintercept = um_pk,lty=2,col=cols[1],lwd=2)+
  geom_segment(aes(x=ep$lagged_dpc,xend=um_pk,y=ep$rt,yend=ep$rt),
               arrow=arrow(),lty=2)+
  geom_text(aes(x=ep$lagged_dpc*1.2,y=0.1,label='rho'),parse=TRUE,cex=12)+
  geom_text(aes(x=1e-3,y=-0.01,label='delta'),parse=TRUE,cex=12)+
  theme_bw(base_size=15)+
  theme(legend.position='none',
        strip.background = element_blank(),
        strip.text.x = element_blank())+
  ggtitle('ERLs as Upper Bounds')

g_erl <- ggplot(xx[scenario %in% c('Unmitigated','Relaxation to ERL')],aes(lagged_dpc,rt))+
  geom_line(aes(color=scenario),lwd=2)+
  scale_x_continuous(name = NULL,trans='log',limits=c(1e-4,4e-3),breaks=10^(-4:-2),labels=NULL)+
  scale_y_continuous('r(t)')+
  scale_color_manual(values=cols[c(1,3)])+
  theme_bw(base_size=15)+
  theme(legend.position='none',
        strip.background = element_blank(),
        strip.text.x = element_blank())+
  ggtitle('Corroboration of ERL')

g_reject <- ggplot(xx[scenario %in% c('Unmitigated','ERL Rejected')],aes(lagged_dpc,rt))+
  geom_line(aes(color=scenario),lwd=2)+
  scale_y_continuous('r(t)')+
  scale_x_continuous(name = NULL,trans='log',limits=c(1e-4,6e-3),breaks=10^(-4:-2))+
  scale_color_manual(values=cols[c(1,4)])+
  theme_bw(base_size=15)+
  theme(legend.position='none')+
  ggtitle('Rejection of ERL')

ggarrange(g_cases,ggarrange(g_stop,g_erl,g_reject,labels = c('B','C','D'),nrow=3),
          labels=c("A",NA),ncol=2)
ggsave('figures/ERLs_for_comparative_epidemiology.png',height=8,width=7,units='in')
