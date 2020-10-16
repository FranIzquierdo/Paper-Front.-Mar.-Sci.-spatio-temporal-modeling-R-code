### FINAL MODEL ###

# Spatio-temporal structure PROGRESSIVE SHARED (binomial and gamma sub-processess)

# Smoothed NOT SHARED covariate Bathymetry RW2 with inla.group = 13

# INLA does the stimation and prediction of the model at the same time

# However, due to the characteristics of our model and the fine mesh resolution

# We will do a prediction step by taking samples of the posterior distribution and constructing the linear predictor


set.seed(1234)
library(INLA)
inla.setOption(scale.model.default=TRUE)

# datasets ----------------------------------------------------------------

data<-read.table(file="./data/data12G.txt", dec=".", header=TRUE) # obs.data
cov<-read.table(file="./data/cov12G.txt", dec=".", header=TRUE) # pred.data

# MESH --------------------------------------------------------------------

coords<-as.matrix(data[,2:1])
bound=inla.nonconvex.hull(as.matrix(data[,2:1]), convex=-0.075, eps=0.05, resolution=40)
mesh <- inla.mesh.2d(loc=coords, boundary=bound, max.edge=c(0.05, 0.45), offset=c(0.1,0.6), cutoff=0.08, min.angle = 0.05)

# SPDE --------------------------------------------------------------------
spde <- inla.spde2.pcmatern(mesh, prior.range=c(0.5, 0.05), prior.sigma=c(0.6, 0.05))

# make.A.matrix ---------------------------------------------------------------
est.temp <- inla.spde.make.A(mesh, loc=matrix(c(data$Long,data$Lat),ncol=2),
                             #index=index.est
                             group=data$Group,
                             n.group=length(unique(data$Group)))


# INDEX bin/con
mesh.index.bin<- inla.spde.make.index("i.bin", n.spde=spde$n.spde, n.group=max(data$Group))
mesh.index.con<- inla.spde.make.index("i.con", n.spde=spde$n.spde, n.group=max(data$Group))

# Stack
est.bin<-inla.stack(data=list(y=cbind(data$NRec_pres,NA)),
                    A=list(est.temp, 1),
                    effects=list(mesh.index.bin,
                                 list(bin.b0=1,
                                      bath.bin=data$IgBath)),
                    tag='est.bin')

est.con<-inla.stack(data=list(y=cbind(NA,ifelse(data$NRec_pres>0,data$NRec,NA))),
                    A=list(est.temp, 1),
                    effects=list(mesh.index.con,
                                 list(con.b0=1,
                                      bath.con=data$IgBath)),
                    tag='est.con')

est.stack<-inla.stack(est.bin,est.con)

# Link
link=c(rep(1,length(est.stack$data$data[,1])/2), rep(2,length(est.stack$data$data[,1])/2))

# Priors
prior_bin=list(theta=list(prior="pc.prec", param=c(0.5,0.01)))

# FINAL MODEL
modelo12_B_IG<-y~-1+f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5)+ con.b0+
  f(bath.bin, model="rw2", hyper=prior_bin)+f(bath.con, copy="bath.bin", fixed=F)+
  f(i.bin,model=spde,group = i.bin.group,control.group = list(model="ar1"))+
  f(i.con,model=spde,group = i.con.group,control.group = list(model="ar1"))

modelo12_B_IG<-inla(modelo12_B_IG,
                   family=c('binomial',"gamma"),
                   #control.fixed=list(expand.factor.strategy="inla"),
                   data=inla.stack.data(est.stack), control.compute=list(dic=TRUE,cpo=TRUE,waic=T,config=FALSE),
                   control.predictor=list(A=inla.stack.A(est.stack), compute=TRUE,link=link),
                   verbose=TRUE, num.threads = 2)

saveRDS(modelo12_B_IG, "./modelo12_B_IG.rds")

