### INLA GROUP FUNCTION ###

# We have selected the FINAL MODEL

# As we just have one covariate, we will use the data without scaling for an easy interpretation

# However, it is worth to mention that when more covariables are in the model, scaling is important to avoid numerical confounding

# Due to the high computational cost of the model, we won't predict by using a stack.pred

# In order to ensure we have exactly the same RW2 knot or groups for the covariate (Bathy), we will set the same minimum and maximum bathymetry for observed and prediction data

# DATA OBS was obtained by extracting the sampling points from the yearly raster covariate maps

# DATA PRED was obtained by extracting the MESH points from the yearly covariate maps

library(INLA)


# Datasets ----------------------------------------------------------------

data<-read.table(file="./data/dataobs12.txt", dec=".", header=TRUE) # data Obs
cov<-read.table(file="./data/datapred12.txt", dec=".", header=TRUE) #data Pred



# INLA.GROUP --------------------------------------------------------------

which.min(data$Bath)
data[1410,]
which.max(data$Bath)
data[741,]

which.min(cov$Bath)
cov[1085,]
which.max(cov$Bath)
cov[1313,]

# We want the same bathy lower and upper limit for observed (data) and prediction (cov) dataset

cov$Bath[which(cov$Bath>max(data$Bath) )]<-NA
cov$Bath[which(cov$Bath<min(data$Bath))]<-NA

# We join both datasets (data and cov) keeping into account the dimension

all<-rbind(data[,c(9,10,11,12)],cov[,c(4,5,6,7)])
igroupbath<-inla.group(all[,1], n=13, idx.only = TRUE) # Number of groups = 13
groupbath<-inla.group(all[,1], n=13)

allin<-cbind(all, igroupbath, groupbath) # join

length(data$Bath) #1556
length(cov$Bath) #22764
length(allin$Bath)#24320

data<-cbind(data, allin$groupbath[1:1556],allin$igroupbath[1:1556])
cov<-cbind(cov, allin$groupbath[1557:24320],allin$igroupbath[1557:24320])

names(data)<-c("Lat","Long","Year","Trawl","Group","Index","NRec","NRec_pres","Bath","Rugo","SBS","SBT","Depth","IgBath","Igroup")
data<-data[,-c(10,11,12)]

names(cov)<-c("Long","Lat","Group","Bath","Rugo","SBS","SBT","IgBath","Igroup")
cov<-cov[,-c(5,6,7)]

length(unique(data$IgBath))
length(unique(cov$IgBath))

summary(data$Bath)
summary(cov$Bath)

write.table(data, file="./data/data12G.txt", append = FALSE, sep = "\t", dec = ".",row.names = TRUE, col.names = TRUE)
write.table(cov, file="./data/cov12G.txt", append = FALSE, sep = "\t", dec = ".",row.names = TRUE, col.names = TRUE)
