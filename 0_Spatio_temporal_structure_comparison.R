### SPATIO-TEMPORAL STRUCTURES COMPARISON ###

# For our response variables (presence/absence and num. hake recruits) we compare the structures:

# Opportunistic, persistent and progressive from Paradinas et al,. 2017.

# We try the 3 structures 1) SHARED and 2) not SHARED between binomial (pre/abs) and gamma (abundance) sub-models.

# After fitting the 6 models, we will make a comparison through WAIC, DIC and LCPO scores

# Finally we select the best spatio-temporal structure

library(INLA)
set.seed(1234)
inla.setOption(scale.model.default = TRUE) ### set scale.model=TRUE, see scale tutorial

# Data obs. ----------------------------------------------------

data<-read.table(file="./data/datosinlaFin12.txt", dec=".", header=TRUE)

# Mesh --------------------------------------------------------------------


coords<-as.matrix(data[,1:2]) # coordinates in lat,long
bound=inla.nonconvex.hull(as.matrix(data[,1:2]), convex=-0.075, eps=0.05, resolution=40)
mesh <- inla.mesh.2d(loc=coords, boundary=bound, max.edge=c(0.05, 0.45), offset=c(0.1,0.6), cutoff=0.08, min.angle = 0.05)

# SPDE --------------------------------------------------------------------

# A sensitivity analysis was performed in terms of spatial range and sigma
spde <- inla.spde2.pcmatern(mesh, prior.range=c(0.5, 0.05), prior.sigma=c(0.6, 0.05))

# S.PERSISTENT --------------------------------------------------------------

# This model includes the TEMPORAL component as identical independent random effects (iid)

# A. Matrix
A.bin <- inla.spde.make.A(mesh, loc=cbind(data$Long, data$Lat))
A.con <- inla.spde.make.A(mesh, loc=cbind(data$Long, data$Lat))

# Stack
est.bin<-inla.stack(data=list(y=cbind(data$NRec_pres,NA)),
                    A=list(A.bin, 1),
                    effects=list(i.bin=1:spde$n.spde,
                                 list(bin.b0=1,
                                      year.bin=data$Year)),
                    tag='est.bin')

est.con<-inla.stack(data=list(y=cbind(NA,ifelse(data$NRec>0,data$NRec,NA))),
                    A=list(A.con, 1),
                    effects=list(i.con=1:spde$n.spde,
                                 list(con.b0=1,
                                      year.con=ifelse(data$NRec>0,data$Year,NA))),
                    tag='est.con')

pers.stack=inla.stack(est.bin,est.con)

# Link
link=c(rep(1,length(pers.stack$data$data[,1])/2),rep(2,length(pers.stack$data$data[,1])/2))

# Model persistent shared
per_S<- inla(y ~ -1 + f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5) + con.b0
            + f(i.bin, model=spde) + f(i.con, copy="i.bin",fixed=F)
            + f(year.bin, model="iid")+ f(year.con, copy="year.bin",fixed=F),
            family=c('binomial',"gamma"),
            data=inla.stack.data(pers.stack),
            control.inla = list(strategy = "gaussian"),
            control.compute=list(dic=TRUE, cpo=TRUE, waic=T, return.marginals=TRUE),
            control.predictor=list(A=inla.stack.A(pers.stack), compute=TRUE,link=link),
            verbose=T , num.threads = 2)
saveRDS(per_S, "./rds/comp_estr/per_S.rds")

# Model persistent not shared
per_NS<- inla(y ~ -1 + f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5) + con.b0
             + f(i.bin, model=spde) + f(i.con, model=spde)
             + f(year.bin, model="iid")+ f(year.con,model="iid"),
             family=c('binomial',"gamma"),
             data=inla.stack.data(pers.stack),
             control.inla = list(strategy = "gaussian"),
             control.compute=list(dic=TRUE, cpo=TRUE, waic=T, return.marginals=TRUE),
             control.predictor=list(A=inla.stack.A(pers.stack), compute=TRUE,link=link),
             verbose=T , num.threads = 2)
saveRDS(per_NS, "./rds/comp_estr/per_NS.rds") 


# S.OPORTUNISTIC ------------------------------------------------------------

# This model includes the TEMPORAL component associated to the SPATIAL component as a replica

# A. Matrix
bin.temp <- inla.spde.make.A(mesh, loc=cbind(data$Long, data$Lat),
                             repl=data$Year,
                             n.repl=length(unique(data$Year)))

con.temp <- inla.spde.make.A(mesh, loc=cbind(data$Long, data$Lat),
                             repl=data$Year,
                             n.repl=length(unique(data$Year)))

# Index
mesh.index.bin<- inla.spde.make.index("i.bin", n.spde=spde$n.spde, n.repl=length(unique(data$Year)))
mesh.index.con<- inla.spde.make.index("i.con", n.spde=spde$n.spde, n.repl=length(unique(data$Year)))

# Stack
bin.stack.repl<-inla.stack(data=list(y=cbind(data$NRec_pres,NA)),
                           A=list(bin.temp, 1),
                           effects=list(mesh.index.bin,
                                        list(bin.b0=1,
                                             prof.bin=inla.group(data$Bath, n=15))),
                           tag='est_bin')

con.stack.repl<-inla.stack(data=list(y=cbind(NA,ifelse(data$NRec>0,data$NRec,NA))),
                           A=list(con.temp, 1),
                           effects=list(mesh.index.con,
                                        list(con.b0=1,
                                             prof.con=inla.group(data$Bath,n=15))),
                           tag='est_con')

repl.stack=inla.stack(bin.stack.repl,con.stack.repl)

# Link
link=c(rep(1,length(repl.stack$data$data[,1])/2),rep(2,length(repl.stack$data$data[,1])/2))

# Model oportunistic shared
opo_S<-inla(y ~ -1 + bin.b0 + con.b0 
            + f(i.bin, model=spde, replicate=i.bin.group) + f(i.con, copy="i.bin", fixed=F),
            family=c('binomial',"gamma"),
            data=inla.stack.data(repl.stack),
            control.inla = list(strategy = "gaussian"),
            control.compute=list(dic=TRUE, cpo=TRUE, waic=T, return.marginals=TRUE),
            control.predictor=list(A=inla.stack.A(repl.stack), compute=TRUE,link=link),
            verbose=T , num.threads = 2)
saveRDS(opo_S, "./rds/comp_estr/opo_S.rds") # OPO Shared

# Model oportunistic not shared
opo_NS<-inla(y ~ -1 + bin.b0 + con.b0 
            + f(i.bin, model=spde, replicate=i.bin.group) + f(i.con, model=spde, replicate=i.con.group),           
            family=c('binomial',"gamma"),
            data=inla.stack.data(repl.stack),
            control.inla = list(strategy = "gaussian"),
            control.compute=list(dic=TRUE, cpo=TRUE, waic=T, return.marginals=TRUE),
            control.predictor=list(A=inla.stack.A(repl.stack), compute=TRUE,link=link),
            verbose=T , num.threads = 2)
saveRDS(opo_NS, "./rds/comp_estr/opo_NS.rds") # OPO Not Shared

# S.PROGRESSIVE --------------------------------------------------------------- 

# The TEMPORAL component is associated to the SPATIAL component and vary in a correlated way between time units (years)

# A. Matrix
bin.temp <- inla.spde.make.A(mesh, loc=matrix(c(data$Long,data$Lat),ncol=2),
                             group=data$Year,
                             n.group=length(unique(data$Year)))
con.temp <- inla.spde.make.A(mesh, loc=matrix(c(data$Long,data$Lat),ncol=2),
                             group=data$Year,
                             n.group=length(unique(data$Year)))

# Index
mesh.index.bin<- inla.spde.make.index("i.bin", n.spde=spde$n.spde, n.group=length(unique(data$Year)))
mesh.index.con<- inla.spde.make.index("i.con", n.spde=spde$n.spde, n.group=length(unique(data$Year)))

# Stack
est.comp.bin<-inla.stack(data=list(y=cbind(data$NRec_pres,NA)),
                         A=list(bin.temp, 1),
                         effects=list(mesh.index.bin,
                                      list(bin.b0=1,data)),
                         tag='est_bin')

est.comp.con<-inla.stack(data=list(y=cbind(NA,ifelse(data$NRec_pres>0,data$NRec,NA))),
                         A=list(con.temp, 1),
                         effects=list(mesh.index.con,
                                      list(con.b0=1,data)),
                         tag='est_con')

est.stack<-inla.stack(est.comp.bin,est.comp.con)

# Link
link=c(rep(1,length(est.stack$data$data[,1])/2), rep(2,length(est.stack$data$data[,1])/2))

# Model progressive shared
prog_S<-y~-1+f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5)+ con.b0+
      f(i.bin,model=spde, replicate= i.bin.group,control.group = list(model="ar1"))+
      f(i.con,copy="i.bin", fixed=F)

prog_S<-inla(prog_S,
      family=c('binomial',"gamma"),
      #control.fixed=list(expand.factor.strategy="inla"),
      data=inla.stack.data(est.stack), control.compute=list(dic=TRUE,cpo=TRUE,waic=T),
      control.predictor=list(A=inla.stack.A(est.stack), compute=TRUE,link=link),
      control.mode=list(theta=prog_S$mode$theta, restart=TRUE),
      verbose=TRUE,control.inla = list(strategy = "gaussian"), num.threads = 2)
saveRDS(prog_S, "./rds/comp_estr/prog_S.rds") # prog Shared

# Model progressive not shared
prog_NS<-y~-1+f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5)+con.b0+
      f(i.bin,model=spde,group = i.bin.group,control.group = list(model="ar1"))+
      f(i.con,model=spde,group = i.con.group,control.group = list(model="ar1"))
      
prog_NS<-inla(prog_NS,     
      family=c('binomial',"gamma"),
      #control.fixed=list(expand.factor.strategy="inla"),
      data=inla.stack.data(est.stack), control.compute=list(dic=TRUE,cpo=TRUE,waic=T),
      control.predictor=list(A=inla.stack.A(est.stack), compute=TRUE,link=link),
      control.mode=list(theta=prog_NS$mode$theta, restart=TRUE),
      verbose=TRUE,control.inla = list(strategy = "gaussian"), num.threads = 2)
saveRDS(prog_NS, "./rds/comp_estr/prog_NS.rds") # prog Not Shared


# STRUCTURE COMPARISON --------------------------------------------------------

mod1<- readRDS("./rds/comp_estr/per_S.rds")# PER S =1
mod2<- readRDS("./rds/comp_estr/per_NS.rds")# PER NS =2
mod3<- readRDS("./rds/comp_estr/opo_S.rds")# opo S =3
mod4<- readRDS("./rds/comp_estr/opo_NS.rds")# opo NS =4
mod5<- readRDS("./rds/comp_estr/prog_S.rds")# prog S =5
mod6<- readRDS("./rds/comp_estr/prog_NS.rds")# prog NS =6


COMP_STR<-data.frame(
  model=c("per_S","per_NS","opo_S","opo_NS","prog_S","prog_NS"),
  dic=c(mod1$dic$dic,mod2$dic$dic,mod3$dic$dic,mod4$dic$dic,mod5$dic$dic,mod6$dic$dic),
  waic=c(mod1$waic$waic,mod2$waic$waic,mod3$waic$waic,mod4$waic$waic,mod5$waic$waic,mod6$waic$waic),
  lcpo=c(-mean(log(mod1$cpo$cpo),na.rm=T),-mean(log(mod2$cpo$cpo),na.rm=T),-mean(log(mod3$cpo$cpo),na.rm=T),-mean(log(mod4$cpo$cpo),na.rm=T),-mean(log(mod5$cpo$cpo),na.rm=T),-mean(log(mod6$cpo$cpo),na.rm=T)),
  failure=c(sum((mod1$cpo$failure>0)*1,na.rm=T),sum((mod2$cpo$failure>0)*1,na.rm=T),sum((mod3$cpo$failure>0)*1,na.rm=T),sum((mod4$cpo$failure>0)*1,na.rm=T),sum((mod5$cpo$failure>0)*1,na.rm=T),sum((mod6$cpo$failure>0)*1,na.rm=T)),
  time=c(mod1$cpu.used[4],mod2$cpu.used[4],mod3$cpu.used[4],mod4$cpu.used[4],mod5$cpu.used[4],mod6$cpu.used[4])
  )
COMP_STR

# The best spatio-temporal structure for our dataset is PROGRESSIVE with NOT SHARED effects between occurrence (binomial) and abundance (gamma) sub-processes
