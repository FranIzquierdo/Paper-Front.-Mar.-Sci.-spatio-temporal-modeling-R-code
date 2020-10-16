### PREDICTION BY COMPONENTS ###

# In this script we will make a prediction by components

# It is worth to mention that INLA does the proces of estimation and prediction at the same time

# This could be easily implemented through a prediction.stack

# However, due to the complexity of the spatio-temporal structure and the finner mesh, we performed the presen approach

# In this approach, we will construct the lineal predictor by suming the correspondent components and applying the des-transform functions to each likelihood (binomial, gamma)


# PREDICTION BY COMPONENTS ------------------------------------------------

shared<- readRDS("./rds/modelo_final/modelo12_B_IG.rds") # charge model rds

dim(shared$summary.random$i.bin)[1]/mesh$n# equivalent to the number of years 

### load pred covariates
covs=read.table("./data/cov12G.txt")

### detransform functions
a_fun=function(x){exp(x)} ## exponent
o_fun=function(x){exp(x)/(1+exp(x))} ### inverse logit


## we will take samples from the marginal posterior distribution of each component into the linear predictor for each node of the mesh

n=100000 ## number of samples to get from each component (the bigger the more accurate predictions)

### spatial
l=length(shared$marginals.random$i.bin)

####depth; put prediciton bathymetry in terms of estimation bathymetry
pred_bathy=covs$IgBath
pred_bathy[which(pred_bathy>max(shared$summary.random$bath.bin[,1]))]=max(shared$summary.random$bath.bin[,1])

### randomise beta0
beta0_bin=inla.rmarginal(n, shared$marginals.fixed$bin.b0) # hist(beta0_bin)
beta0_abun=inla.rmarginal(n, shared$marginals.fixed$con.b0) # hist(beta0_abun)

### create bathy list with random values drawed from marginals of the effect at each bathymetry
b.eff_bin=b.eff_abun=list()
for(i in 1:length(shared$summary.random$bath.bin[,1])){
  b.eff_bin[[i]]=inla.rmarginal(n, shared$marginals.random$bath.bin[[i]])
}
for(i in 1:length(shared$summary.random$bath.con[,1])){
  b.eff_abun[[i]]=inla.rmarginal(n, shared$marginals.random$bath.con[[i]])
}

### bathy index; get the indexation of the bathyemtry to do prediction at each node of the mesh
bathy.idx_bin=bathy.idx_abun=c()
for(i in 1:mesh$n){
  if(is.na(pred_bathy[i])){
    bathy.idx_bin[i]=NA
    bathy.idx_abun[i]=NA
  }else{
    bathy.idx_bin[i]=which.min(abs(shared$summary.random$bath.bin[,1]-pred_bathy[i]))
    bathy.idx_abun[i]=which.min(abs(shared$summary.random$bath.con[,1]-pred_bathy[i]))}  
}

# MEDIAN prediction -----------------------------------------------------------

# We selected the median of the posterior distribution for occurrence and abundance for being more stable than the mean

abun_mean=occur_mean=c()
for (i in 1:(l)){
  year=ceiling(i/mesh$n)
  idx=i-(year-1)*mesh$n
  
  if(is.na(pred_bathy[idx])){
    abun_mean[i]=NA
    occur_mean[i]=NA
  }else{
    ### sum each component on the linear predictor and des-transform with the correspondent function
    occur_mean[i]=median(o_fun(beta0_bin + b.eff_bin[[bathy.idx_bin[idx]]] + inla.rmarginal(n,shared$marginals.random$i.bin[[i]])))
    abun_mean[i]=median(a_fun(beta0_abun + b.eff_abun[[bathy.idx_abun[idx]]] + inla.rmarginal(n,shared$marginals.random$i.con[[i]])))
  }
  if(i%in%seq(0,l,mesh$n)){print(paste("Finished year",year))}
}

occur_mean ## median for occurrences
abun_mean ## median for abundances

write.table(occur_mean, file="./data/meanpredBin.txt", append = FALSE, sep = "\t", dec = ".",row.names = TRUE, col.names = TRUE)
write.table(abun_mean, file="./data/meanpredGam.txt", append = FALSE, sep = "\t", dec = ".",row.names = TRUE, col.names = TRUE)


# QUANTILE prediction ---------------------------------------------------------------------

# A recommended representation to display uncertainty are the posterior predictive quantiles as well as the difference between them

# In this case, we calculated the first and third QUARTILES of the posterior predictive distribution for occurrence and abundance

abun_q=occur_q=c()
for (i in 1:(l)){
  year=ceiling(i/mesh$n)
  idx=i-(year-1)*mesh$n
  
  if(is.na(pred_bathy[idx])){
    abun_q[i]=NA
    occur_q[i]=NA
  }else{
    occur_q[i]=quantile(o_fun(beta0_bin + b.eff_bin[[bathy.idx_bin[idx]]] + inla.rmarginal(n,shared$marginals.random$i.bin[[i]])), prob=0.75) ### change probability
    abun_q[i]=quantile(a_fun(beta0_abun + b.eff_abun[[bathy.idx_abun[idx]]] + inla.rmarginal(n,shared$marginals.random$i.con[[i]])), prob=0.75)
  }
  if(i%in%seq(0,l,mesh$n)){print(paste("Finished year",year))}
}

occur_q## occurrence
abun_q ## abundance


write.table(occur_sd, file="./data/sdpredBinq75.txt", append = FALSE, sep = "\t", dec = ".",row.names = TRUE, col.names = TRUE)
write.table(abun_sd, file="./data/sdpredGamq75.txt", append = FALSE, sep = "\t", dec = ".",row.names = TRUE, col.names = TRUE)
