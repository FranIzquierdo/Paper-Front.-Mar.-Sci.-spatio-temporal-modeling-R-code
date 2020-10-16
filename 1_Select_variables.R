### ENVIRONMENTAL VARIABLE SELECTION FOR STRUCTURE PROGRESSIVE WITH NOT SHARED EFFECTS ###

# In this step we will fit the model progressive with the covariates SBT, SBS and Bathymetry

# All covariates will be included as smoothed effects (Random walk of second order effect)

# First, we will include each covariate through SHARED effects (between binomial and gamma) and the possible combinations between them

# Second, we will include each covariate through NOT SHARED effects (between binomial and gamma) and the possible double combinations between them

# We will use vague PRIORS from other studies for Binomial and Gamma

# Our dataset contains all the scaled covariates

# We need to apply inla.group (different number tried) into the Stack step because the sampling points are too close for the RW2 

# To avoid a long script (as all the covariates are fitted in the same way), just the example with bathymetry is included

set.seed(1234) 
library(INLA)
inla.setOption(scale.model.default = TRUE) ### set scale.model=TRUE, see scale tutorial

# Data obs  ----------------------------------------------------

data<-read.table(file="./data/datosinlaFin12.txt", dec=".", header=TRUE)

# Mesh--------------------------------------------------------------------

coords<-as.matrix(data[,2:1])
bound=inla.nonconvex.hull(as.matrix(data[,2:1]), convex=-0.075, eps=0.05, resolution=40)
mesh <- inla.mesh.2d(loc=coords, boundary=bound, max.edge=c(0.05, 0.45), offset=c(0.1,0.6), cutoff=0.08, min.angle = 0.05)

# SPDE --------------------------------------------------------------------

spde <- inla.spde2.pcmatern(mesh, prior.range=c(0.5, 0.05), prior.sigma=c(0.6, 0.05))

# Stack ---------------------------------------------------------------

est.temp <- inla.spde.make.A(mesh, loc=matrix(c(data$Long,data$Lat),ncol=2),
                             #index=index.est
                             group=data$Group,
                             n.group=length(unique(data$Group)))

# INDEX bin/con
mesh.index.bin<- inla.spde.make.index("i.bin", n.spde=spde$n.spde, n.group=max(data$Group))
mesh.index.con<- inla.spde.make.index("i.con", n.spde=spde$n.spde, n.group=max(data$Group))

est.bin<-inla.stack(data=list(y=cbind(data$NRec_pres,NA)),
                    A=list(est.temp, 1),
                    effects=list(mesh.index.bin,
                                 list(bin.b0=1,
                                      bath.bin=inla.group(data$Bath, n=15))),
                                      tag='est.bin')
#
                                 #sbt.bin=inla.group(data$SST, n=15),
                                 #sbs.bin=inla.group(data$SSS, n=15))),
                                 
est.con<-inla.stack(data=list(y=cbind(NA,ifelse(data$NRec_pres>0,data$NRec,NA))),
                    A=list(est.temp, 1),
                    effects=list(mesh.index.con,
                                 list(con.b0=1,
                                      bath.con=inla.group(data$Bath, n=15))),
                                                         tag='est.con')
                                 #sbt.con=inla.group(data$SST, n=15),
                                 #sbs.con=inla.group(data$SSS,n=15))),
est.stack<-inla.stack(est.bin,est.con)

# Link
link=c(rep(1,length(est.stack$data$data[,1])/2), rep(2,length(est.stack$data$data[,1])/2))

# Priors Rw2
prior_bin=list(theta=list(prior="pc.prec", param=c(0.5,0.01))) # Good prior for binomial bathy rw2 effect.
prior_con=list(theta=list(prior="pc.prec", param=c(1,0.01))) # Good prior for gamma bathy rw2 effect.


# Model Covariate RW2 shared
prog_est_B_Sin<-y~-1+f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5)+ con.b0+
  f(bath.bin, model="rw2", hyper=prior_bin)+f(bath.con, copy="bath.bin", fixed=F)+
  f(i.bin,model=spde,group = i.bin.group,control.group = list(model="ar1"))+
  f(i.con,model=spde,group = i.con.group,control.group = list(model="ar1"))
saveRDS(prog_est_B, "./rds/selec_var/prog_est_B.rds")

prog_est_B_Sin<-inla(prog_est_B_Sin,
              family=c('binomial',"gamma"),
              #control.fixed=list(expand.factor.strategy="inla"),
              data=inla.stack.data(est.stack), control.compute=list(dic=TRUE,cpo=TRUE,waic=T),
              control.predictor=list(A=inla.stack.A(est.stack), compute=TRUE,link=link),
              verbose=TRUE,control.inla = list(strategy = "gaussian"), num.threads = 2)

saveRDS(prog_est_B_Sin, "./rds/selec_var/prog_est_B.rds")

# Model Covariate rw2 NOT shared

prog_est_B_NS<-y~-1+f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5)+ con.b0+
  f(bath.bin, model="rw2", hyper=prior_bin)+f(bath.con, model="rw2", hyper=prior_con)+
  f(i.bin,model=spde,group = i.bin.group,control.group = list(model="ar1"))+
  f(i.con,model=spde,group = i.con.group,control.group = list(model="ar1"))

prog_est_B_NS<-inla(prog_est_B_NS,
                    family=c('binomial',"gamma"),
                    #control.fixed=list(expand.factor.strategy="inla"),
                    data=inla.stack.data(est.stack), control.compute=list(dic=TRUE,cpo=TRUE,waic=T),
                    control.predictor=list(A=inla.stack.A(est.stack), compute=TRUE,link=link),
                    verbose=TRUE,control.inla = list(strategy = "gaussian"), num.threads = 2)

saveRDS(prog_est_B_NS, "./rds/selec_var/prog_est_B_NS.rds")

# Model Comparison --------------------------------------------------------

mod1<- readRDS("./rds/selec_var/prog_est_B.rds")
mod2<- readRDS("./rds/selec_var/prog_est_S.rds")
mod3<- readRDS("./rds/selec_var/prog_est_T.rds")
mod4<- readRDS("./rds/selec_var/prog_est_BS.rds")
mod5<- readRDS("./rds/selec_var/prog_est_BT.rds")
mod6<- readRDS("./rds/selec_var/prog_est_ST.rds")
mod7<- readRDS("./rds/selec_var/prog_est_B_NS.rds")
mod8<- readRDS("./rds/selec_var/prog_est_S_NS.rds")
mod9<- readRDS("./rds/selec_var/prog_est_T_NS.rds")
mod10<- readRDS("./rds/selec_var/prog_est_BS_NS.rds")
mod11<- readRDS("./rds/selec_var/prog_est_BT_NS.rds")
mod12<- readRDS("./rds/selec_var/prog_est_ST_NS.rds")


SEL_ENVAR<-data.frame(
  model=c("prog_B_s","prog_S_s","prog_T_s","prog_BS_s","prog_BT_s","prog_ST_s","prog_B_ns","prog_S_ns","prog_T_ns","prog_BS_ns","prog_BT_ns","prog_ST_ns"),
  dic=c(mod1$dic$dic,mod2$dic$dic,mod3$dic$dic,mod4$dic$dic,mod5$dic$dic,mod6$dic$dic,mod7$dic$dic,mod8$dic$dic,mod9$dic$dic,mod10$dic$dic,mod11$dic$dic,mod12$dic$dic),
  waic=c(mod1$waic$waic,mod2$waic$waic,mod3$waic$waic,mod4$waic$waic,mod5$waic$waic,mod6$waic$waic,mod7$waic$waic,mod8$waic$waic,mod9$waic$waic,mod10$waic$waic,mod11$waic$waic,mod12$waic$waic),
  lcpo=c(-mean(log(mod1$cpo$cpo),na.rm=T),-mean(log(mod2$cpo$cpo),na.rm=T),-mean(log(mod3$cpo$cpo),na.rm=T),-mean(log(mod4$cpo$cpo),na.rm=T),-mean(log(mod5$cpo$cpo),na.rm=T),-mean(log(mod6$cpo$cpo),na.rm=T),-mean(log(mod7$cpo$cpo),na.rm=T),-mean(log(mod8$cpo$cpo),na.rm=T),-mean(log(mod9$cpo$cpo),na.rm=T),-mean(log(mod10$cpo$cpo),na.rm=T),-mean(log(mod11$cpo$cpo),na.rm=T),-mean(log(mod12$cpo$cpo),na.rm=T)),
  failure=c(sum((mod1$cpo$failure>0)*1,na.rm=T),sum((mod2$cpo$failure>0)*1,na.rm=T),sum((mod3$cpo$failure>0)*1,na.rm=T),sum((mod4$cpo$failure>0)*1,na.rm=T),sum((mod5$cpo$failure>0)*1,na.rm=T),sum((mod6$cpo$failure>0)*1,na.rm=T),sum((mod7$cpo$failure>0)*1,na.rm=T),sum((mod8$cpo$failure>0)*1,na.rm=T),sum((mod9$cpo$failure>0)*1,na.rm=T),sum((mod10$cpo$failure>0)*1,na.rm=T),sum((mod11$cpo$failure>0)*1,na.rm=T),sum((mod12$cpo$failure>0)*1,na.rm=T)),
  time=c(mod1$cpu.used[4],mod2$cpu.used[4],mod3$cpu.used[4],mod4$cpu.used[4],mod5$cpu.used[4],mod6$cpu.used[4],mod7$cpu.used[4],mod8$cpu.used[4],mod9$cpu.used[4],mod10$cpu.used[4],mod11$cpu.used[4],mod12$cpu.used[4])
)

SEL_ENVAR

# SHARED covariate models have shown a better godness of fit compared to NOT SHARED covariate models

# For our particular case, the 3 best models are Bath + SBS, Bath + SBT and Bath

# We select as the BEST MODEL the one which includes just Bathymetry, because SBT and SBS improvement is rather negligible

