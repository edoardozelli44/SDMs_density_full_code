## In the following section example R code is provided for RF ##

## ===== Modelling the spatial distribution of density of VMEs indicator taxa in New Zealand waters  ===== ##

##------------------------------  Edoardo Zelli
##------------------------------  Start date : September 2022
##------------------------------  End date : August 2023

# Load libraries
library(Metrics) 
library(randomForest) 
library(raster)
library(rgdal)
library(ROCR)
library(pROC) 
library(hydroGOF) 
library(pdp) 
library(ggplot2) 
library(gridExtra)
library(AER)
library(MASS)
library(corrplot)

# Load biological data
Bio.D <- read.csv("C:/Users/ez14/OneDrive - The University of Waikato/EDO/Chapter 1/Training_Fabrice/DTIS_Clipped_Proj.csv", sep=";")

# # Load env data (1)
# # Extracting env data from CSV file and converting into RASTERS
# # e.g., Present_Bathy
# setwd("C:/Users/ez14/OneDrive - The University of Waikato/EDO/Chapter 1/Training_Fabrice/Present_Env_predictors")
# Present <- read.csv("C:/Users/ez14/OneDrive - The University of Waikato/EDO/Chapter 1/Training_Fabrice/ENV/PRESENT_ENV.csv")
# library(tidyverse)
# bathy = Present %>% select(x,y,depth)
# head(bathy)
# lon_ex = range(bathy$x)
# lat_ex = range(bathy$y)
# raster_BATHY = raster::raster()
# class(raster_BATHY)
# raster::crs(raster_BATHY) <- "EPSG:3994" # set the CRS
# raster::extent(raster_BATHY) = raster::extent(lon_ex[1], lon_ex[2],lat_ex[1], lat_ex[2])
# raster::res(raster_BATHY) = 1
# raster_BATHY
# dim(raster_BATHY) = c(length(unique(bathy$y)),length(unique(bathy$x)))
# raster::values(raster_BATHY) = bathy$depth
# raster::plot(raster_BATHY)
# writeRaster(raster_BATHY, 'Present_Bathy.tif', overwrite=TRUE)

# Load env data (2) - (leading all rasters in one go) - NB: not other files in the folder, only the rasters!

# Present
setwd("C:/Users/ez14/OneDrive - The University of Waikato/EDO/Chapter 1/Training_Fabrice/Present_Env_predictors_EEZ")
f1km_Present <- list.files(getwd())
ras1km_Present <- lapply(f1km_Present,raster) # load as raster
PredStack1km_Present <- stack(ras1km_Present) # creating a preidictor stack (a raster stack of env. variables used)
# convert list of rasters to stack: turn rasters (stack) into a dataframe (needed for spatial prediction)
Pred_1km_Present <-  na.omit(as.data.frame(PredStack1km_Present, xy = T))
colnames(Pred_1km_Present) <- c("X","Y","Aragonite","Bathy","Bot_sal","Bpi_broad","Bpi_fine",
                                "Calcite","Detr_flux","Nitro","OXY_C","Prof","Slope")

# # plot environemntal variables (rasters) to check 
plot(ras1km_Present[[2]])
# plot(ras1km_Present[[11]])

# SSP2
setwd("C:/Users/ez14/OneDrive - The University of Waikato/EDO/Chapter 1/Training_Fabrice/Future_Env_predictors_SSP2_EEZ")
f1km_SSP2 <- list.files(getwd())
ras1km_SSP2 <- lapply(f1km_SSP2,raster) 
PredStack1km_SSP2 <- stack(ras1km_SSP2)
Pred_1km_SSP2<- na.omit(as.data.frame(PredStack1km_SSP2, xy = T))
colnames(Pred_1km_SSP2) <- c("X","Y","Aragonite","Bathy","Bot_sal","Bpi_broad","Bpi_fine",
                             "Calcite","Detr_flux","Nitro","OXY_C","Prof","Slope")

# SSP3
setwd("C:/Users/ez14/OneDrive - The University of Waikato/EDO/Chapter 1/Training_Fabrice/Future_Env_predictors_SSP3_EEZ")

f1km_SSP3 <- list.files(getwd())
ras1km_SSP3 <- lapply(f1km_SSP3,raster) 
PredStack1km_SSP3 <- stack(ras1km_SSP3)
Pred_1km_SSP3<- na.omit(as.data.frame(PredStack1km_SSP3, xy = T))
colnames(Pred_1km_SSP3) <- c("X","Y","Aragonite","Bathy","Bot_sal","Bpi_broad","Bpi_fine",
                             "Calcite","Detr_flux","Nitro","OXY_C","Prof","Slope")

# to project Lat Long to X AND y 
MPI.proj <- crs(PredStack1km_Present[[3]])
myproj <- CRS("+proj=merc +lat_ts=-41 +lon_0=100 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

# transform the x y 
xy <- as.data.frame(cbind(Y=Bio.D$Lat_mean, X= Bio.D$Lon_mean))
coordinates(xy) <- c("X","Y")
proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  # for example
xy <- spTransform(xy, MPI.proj)
xy <-  as.data.frame(cbind(X = xy@coords[,1], Y = xy@coords[,2]))

#extract env info (from our predictorStacks - e.g., Bathy or OXY)
Present_env <- as.data.frame(raster::extract(x=PredStack1km_Present, y=xy, method='simple', buffer=NULL, small=FALSE, cellnumbers=FALSE))
colnames(Present_env) <- c("Aragonite","Bathy","Bot_sal","Bpi_broad","Bpi_fine",
                           "Calcite","Detr_flux","Nitro","OXY_C","Prof","Slope")

imp.var <- c("Aragonite","Bathy","Bpi_broad","Bpi_fine","Detr_flux","Nitro","OXY_C","Prof","Slope")

imp.var <- imp.var[order(imp.var)]
order(imp.var) # Order (range alphabetically) the env.varibles

################################################################################

# extract info for species of interest from biological data
taxa.name <- Bio.D[,c("Lon_mean","Lat_mean", "Goniocorella.dumosa")]

# make a PA column (to create a PA column in your data frame all you do is adding a $)
taxa.name$taxa.name.PA <- 0
colnames(taxa.name) <- c("Lon_mean","Lat_mean", "taxa.name.D","taxa.name.PA")
taxa.name[taxa.name$taxa.name.D>0,]$taxa.name.PA <- 1

# for RF presence-absence needs to be saved as factor!!
# for other methods it's fine to leave as 0,1
taxa.name$taxa.name.PA <- as.factor(taxa.name$taxa.name.PA)

# format the density
# Conforming dimension (1000 m2) - taxa.name.D*1000
taxa.name$taxa.name.D = taxa.name$taxa.name.D*1000
# Replacing 0 with NA 
taxa.name[taxa.name$taxa.name.D==0,]$taxa.name.D <- NA
# LOG transform the density
taxa.name$taxa.name.D = log(taxa.name$taxa.name.D)
hist(taxa.name$taxa.name.D, main = "Density", ylab="Frequency", xlab="n? of individuals")

#Final object for modelling (AB)
taxa.name.F <- cbind(xy,taxa.name,Present_env)

# removing NAs values
taxa.name.F <- subset(taxa.name.F, !is.na(taxa.name.F$Bathy&taxa.name.F$Prof)) # by colunms (with NAs)
# rows_to_delete <- c(115,445,446,453,457,465,466,467,471,472,683,684,685:688,692:697,702,703,710,711) # by rows
# taxa.name.F <- taxa.name.F[-rows_to_delete, ] 

# SAVE OBJECTS
setwd("C:/Users/ez14/OneDrive - The University of Waikato/EDO/Chapter 1/Training_Fabrice/R_files/taxa.name_final")
save(taxa.name.F, file = "taxa.name.F.RData")
save(Pred_1km_Present, file = "Pred_1km_Present.RData")
save(Pred_1km_SSP2, file = "Pred_1km_SSP2.RData")
save(Pred_1km_SSP3, file = "Pred_1km_SSP3.RData")
load("Pred_1km_Present.RData")
load("Pred_1km_SSP2.RData")
load("Pred_1km_SSP3.RData")
load("taxa.name.F.RData")

################################################################################

# MULTICOLLINEARITY

# Assess collinearity with Pearson correlation coefficient
taxa.name.F$taxa.name.P <- 0
taxa.name.F[taxa.name.F$taxa.name.PA==1,]$taxa.name.P <- 1
taxa.name.F.1 <- taxa.name.F[taxa.name.F$taxa.name.P==1, ]
taxa.name.F.2  <- taxa.name.F.1[, -c(1,2,3,4,5,6,9,12,18)] # remove the columns and the variables that are not important
coeff<-cor(taxa.name.F.2, method="pearson",use="pairwise.complete.obs")
coeff[which(abs(coeff[])>0.8)] # no correlation > 0.8
# Correlation plot
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(coeff,method="number", col=col(200),tl.cex=0.8,number.cex = 0.7,
         tl.col="black", tl.srt=45)

png("Collinearity_taxa.name.png")
corrplot(coeff,method="number", col=col(200),tl.cex=0.8,number.cex = 0.7,
         tl.col="black", tl.srt=45)
dev.off()

###=============================================================================
###=============================================================================

# BOOTSTRAPPING 
set.seed(2022)
n.boot <- 150

# setup objects to save information from bootstraps

# for AUC and P
deviance_mat <- array(0, c(n.boot,3)) # for saving model fit metrics

# for variable relative influence 
var_mat_PA <- array(0, c(length(imp.var),n.boot)) 
var_mat_AB <- array(0, c(length(imp.var),n.boot))
rownames(var_mat_PA) <- c("Aragonite","Bathy","Bpi_broad","Bpi_fine","Detr_flux","Nitro","OXY_C","Prof","Slope")
rownames(var_mat_AB) <- c("Aragonite","Bathy","Bpi_broad","Bpi_fine","Detr_flux","Nitro","OXY_C","Prof","Slope")

# for spatial prediction

# Present
length.map_present <- nrow(Pred_1km_Present)
boot_mat.PA_present <- array(0, c(length.map_present,n.boot)) # matrix where rows will store predictions for each bootstrap (columns)
boot_mat.AB_present <- array(0, c(length.map_present,n.boot)) # matrix where rows will store predictions for each bootstrap (columns)

# SSP2
length.map_SSP2 <- nrow(Pred_1km_SSP2)
boot_mat.PA_SSP2 <- array(0, c(length.map_SSP2,n.boot)) # matrix where rows will store predictions for each bootstrap (columns)
boot_mat.AB_SSP2 <- array(0, c(length.map_SSP2,n.boot)) # matrix where rows will store predictions for each bootstrap (columns)

# SSP3
length.map_SSP3 <- nrow(Pred_1km_SSP3)
boot_mat.PA_SSP3 <- array(0, c(length.map_SSP3,n.boot)) # matrix where rows will store predictions for each bootstrap (columns)
boot_mat.AB_SSP3 <- array(0, c(length.map_SSP3,n.boot)) # matrix where rows will store predictions for each bootstrap (columns)

# LOOP

for(i in 1:n.boot){
  
  # Data partitioning - creating training and evaluation (test) presence-absence and abundance dataframes
  train_taxa.name.F.raw <- sample(seq_len(nrow(taxa.name.F)), size = nrow(taxa.name.F), replace = T) # index of rows for of presence data
  train_taxa.name.F <- taxa.name.F[train_taxa.name.F.raw, ] #train for PA
  train_taxa.name.F_1 <- na.omit(train_taxa.name.F)  #train for Abundance
  test_taxa.name.F <- taxa.name.F[-train_taxa.name.F.raw, ] # test for PA
  test_taxa.name.F_1 <- na.omit(test_taxa.name.F)    # test for Abundance 
  
  # FIT MODEL
  
  model_PA <- tuneRF(x = train_taxa.name.F [,imp.var],
                     y = train_taxa.name.F$taxa.name.PA,
                     mtryStart = 2,
                     ntreeTry = 1000,
                     stepFactor = 2,
                     improve = 0.005,
                     trace = F, plot = F, doBest = T) # tuning for PA
  
  model_AB <- tuneRF(x = train_taxa.name.F_1[,imp.var],
                     y = train_taxa.name.F_1$taxa.name.D,
                     mtryStart = 2,
                     ntreeTry = 1000,
                     stepFactor = 2,
                     improve = 0.005,
                     trace = F, plot = F, doBest = T) # tuning for Abundance
  
  ################################################################################ 
  
  # MODEL FIT METRICS 
  
  # Presence-absence (AUC) 
  
  pred.PA <- predict(model_PA, test_taxa.name.F[,imp.var], type = "prob")
  pred.1 <- prediction(pred.PA[,2], test_taxa.name.F$taxa.name.PA) 
  #evaluation <- performance(pred.1, "acc") 
  #roc <- performance(pred.1, "tpr", "fpr")
  
  # results AUC
  auc <- performance(pred.1, "auc")
  AUC <- unlist(slot(auc, "y.values")) 
  deviance_mat[i,1] <- AUC

  # TSS
  actuals <- test_taxa.F$taxa.PA
  predicted <- predict(model_PA, test_taxa.F[,imp.var], type = "prob")[,2]
  myROC <- pROC::roc(actuals, predicted, quiet = T)
  myROC$auc # this is the AUC
  Sens_spec <- pROC::coords(myROC, x="best", input="threshold", best.method="youden",transpose = FALSE)
  TSS <- mean(Sens_spec[,2]) + mean(Sens_spec[,2])-1 
  deviance_mat[i,2] <- TSS
    
  # Abundance (Pearson's correlation)
  
  # observed abundance = test_taxa.name.F_1$taxa.name.D
  # predicted abundance = Abundance model used to predict to train_taxa.name.F_1
  pred.AB <- predict(model_AB, test_taxa.name.F_1[,imp.var], type = "response")
  actual <- as.vector(test_taxa.name.F_1$taxa.name.D)
  predicted <- as.vector(pred.AB)
  P <- cor(actual, predicted)
  deviance_mat[i,3] <- P
  
  ################################################################################
  
  # VARIABLES RELATIVE INFLUENCE
  
  imp.pred.PA <- importance(model_PA)
  imp.pred.AB <- importance(model_AB)
  for (m in 1:length(imp.pred.PA)){imp.pred.PA[m] <- (imp.pred.PA[m]/sum(imp.pred.PA))*100}
  for (m in 1:length(imp.pred.AB)){imp.pred.AB[m] <- (imp.pred.AB[m]/sum(imp.pred.AB))*100}
  var_mat_PA[,i] <- imp.pred.PA
  var_mat_AB[,i] <- imp.pred.AB
  
  ################################################################################  
  
  # SPATIAL PREDICTIONS
  
  # PRESENT  
  
  # PA
  pa_present <- predict(model_PA, Pred_1km_Present, type = "prob")[,2]
  boot_mat.PA_present[,i] <- pa_present
  # AB
  abu_present <- predict(model_AB, Pred_1km_Present, type = "response")
  boot_mat.AB_present[,i] <- pa_present * exp(abu_present)
  
  #SSP2
  
  # PA
  pa_SSP2 <- predict(model_PA, Pred_1km_SSP2, type = "prob")[,2]
  boot_mat.PA_SSP2[,i] <- pa_SSP2
  # AB
  abu_SSP2 <- predict(model_AB, Pred_1km_SSP2, type = "response")
  boot_mat.AB_SSP2[,i] <- pa_SSP2 * exp(abu_SSP2)
  
  # SSP3
  
  # PA
  pa_SSP3 <- predict(model_PA, Pred_1km_SSP3, type = "prob")[,2]
  boot_mat.PA_SSP3[,i] <- pa_SSP3
  # AB
  abu_SSP3 <- predict(model_AB, Pred_1km_SSP3, type = "response")
  boot_mat.AB_SSP3[,i] <- pa_SSP3 * exp(abu_SSP3)
  
  # keep a track of iterations
  print(paste("Iteration ", i, " out of ", n.boot, sep =""))
}

# calculate mean and SD of model fit metrics
mean_AUC <- mean(deviance_mat[,1])
mean_TSS<- mean(deviance_mat[,2])
mean_p <- mean(deviance_mat[,3])
SD_AUC <- sd(deviance_mat[,1])
SD_TSS <- sd(deviance_mat[,2])
SD_p <- sd(deviance_mat[,3])

# calculate mean and SD of variables relative influence 
mean_imp.var_PA <- apply(var_mat_PA, 1, mean)
mean_imp.var_AB <- apply(var_mat_AB, 1, mean)
SD_imp.var_PA <- apply(var_mat_PA, 1, sd)
SD_imp.var_AB <- apply(var_mat_AB, 1, sd)

# calculate mean and SD of spatial prediction  

# Present
boot_mat.PA_mean_present <- apply(boot_mat.PA_present, 1, mean)
boot_mat.AB_mean_present <- apply(boot_mat.AB_present, 1, mean)
boot_mat.PA_SD_present <- apply(boot_mat.PA_present, 1, sd)
boot_mat.AB_SD_present <- apply(boot_mat.AB_present, 1, sd)

# SSP2
boot_mat.PA_mean_SSP2 <- apply(boot_mat.PA_SSP2, 1, mean)
boot_mat.AB_mean_SSP2 <- apply(boot_mat.AB_SSP2, 1, mean)
boot_mat.PA_SD_SSP2 <- apply(boot_mat.PA_SSP2, 1, sd)
boot_mat.AB_SD_SSP2 <- apply(boot_mat.AB_SSP2, 1, sd)

# SSP3
boot_mat.PA_mean_SSP3 <- apply(boot_mat.PA_SSP3, 1, mean)
boot_mat.AB_mean_SSP3 <- apply(boot_mat.AB_SSP3, 1, mean)
boot_mat.PA_SD_SSP3 <- apply(boot_mat.PA_SSP3, 1, sd)
boot_mat.AB_SD_SSP3 <- apply(boot_mat.AB_SSP3, 1, sd)

# SAVE OBJECTS
setwd("C:/Users/ez14/OneDrive - The University of Waikato/EDO/Chapter 1/Training_Fabrice/R_files/taxa.name_final")
save(boot_mat.PA_present, file = "boot_mat.PA_present.taxa.name.RData")
save(boot_mat.AB_present, file = "boot_mat.AB_present.taxa.name.RData")
save(boot_mat.PA_mean_present, file = "boot_mat.PA_mean_present.taxa.name.RData")
save(boot_mat.AB_mean_present, file = "boot_mat.AB_mean_present.taxa.name.RData")

save(boot_mat.PA_SSP2, file = "boot_mat.PA_SSP2.taxa.name.RData")
save(boot_mat.AB_SSP2, file = "boot_mat.AB_SSP2.taxa.name.RData")
save(boot_mat.PA_mean_SSP2, file = "boot_mat.PA_mean_SSP2.taxa.name.RData")
save(boot_mat.AB_mean_SSP2, file = "boot_mat.AB_mean_SSP2.taxa.name.RData")

save(boot_mat.PA_SSP3, file = "boot_mat.PA_SSP3.taxa.name.RData")
save(boot_mat.AB_SSP3, file = "boot_mat.AB_SSP3.taxa.name.RData")
save(boot_mat.PA_mean_SSP3, file = "boot_mat.PA_mean_SSP3.taxa.name.RData")
save(boot_mat.AB_mean_SSP3, file = "boot_mat.AB_mean_SSP3.taxa.name.RData")

###=============================================================================
###=============================================================================
###=============================================================================

# RESULTS - calculations (model fit & predictions)

# turn mean and uncertainty (SD) predictions into raster - PA & ABU 
# NB: no need hurdle because already done in the loop

# PRESENT

# means
pred.PA_mean_present <- rasterFromXYZ(data.frame(x = Pred_1km_Present[,1],
                                                 y = Pred_1km_Present[,2],
                                                 z = boot_mat.PA_mean_present),
                                      crs = myproj)

pred.AB_mean_present <- rasterFromXYZ(data.frame(x = Pred_1km_Present[,1],
                                                 y = Pred_1km_Present[,2],
                                                 z = boot_mat.AB_mean_present),
                                      crs = myproj)
# SD
pred.PA_SD_present <- rasterFromXYZ(data.frame(x = Pred_1km_Present[,1],
                                               y = Pred_1km_Present[,2],
                                               z = boot_mat.PA_SD_present),
                                    crs = myproj)

pred.AB_SD_present <- rasterFromXYZ(data.frame(x = Pred_1km_Present[,1],
                                               y = Pred_1km_Present[,2],
                                               z = boot_mat.AB_SD_present),
                                    crs = myproj)

# SSP2

# means
pred.PA_mean_SSP2 <- rasterFromXYZ(data.frame(x = Pred_1km_SSP2[,1],
                                              y = Pred_1km_SSP2[,2],
                                              z = boot_mat.PA_mean_SSP2),
                                   crs = myproj)

pred.AB_mean_SSP2 <- rasterFromXYZ(data.frame(x = Pred_1km_SSP2[,1],
                                              y = Pred_1km_SSP2[,2],
                                              z = boot_mat.AB_mean_SSP2),
                                   crs = myproj)
# SD
pred.PA_SD_SSP2 <- rasterFromXYZ(data.frame(x = Pred_1km_SSP2[,1],
                                            y = Pred_1km_SSP2[,2],
                                            z = boot_mat.PA_SD_SSP2),
                                 crs = myproj)

pred.AB_SD_SSP2 <- rasterFromXYZ(data.frame(x = Pred_1km_SSP2[,1],
                                            y = Pred_1km_SSP2[,2],
                                            z = boot_mat.AB_SD_SSP2),
                                 crs = myproj)
# SSP3

# means
pred.PA_mean_SSP3 <- rasterFromXYZ(data.frame(x = Pred_1km_SSP3[,1],
                                              y = Pred_1km_SSP3[,2],
                                              z = boot_mat.PA_mean_SSP3),
                                   crs = myproj)

pred.AB_mean_SSP3 <- rasterFromXYZ(data.frame(x = Pred_1km_SSP3[,1],
                                              y = Pred_1km_SSP3[,2],
                                              z = boot_mat.AB_mean_SSP3),
                                   crs = myproj)
# SD
pred.PA_SD_SSP3 <- rasterFromXYZ(data.frame(x = Pred_1km_SSP3[,1],
                                            y = Pred_1km_SSP3[,2],
                                            z = boot_mat.PA_SD_SSP3),
                                 crs = myproj)

pred.AB_SD_SSP3 <- rasterFromXYZ(data.frame(x = Pred_1km_SSP3[,1],
                                            y = Pred_1km_SSP3[,2],
                                            z = boot_mat.AB_SD_SSP3),
                                 crs = myproj)

##==============================================================================

# CLIPPING TO DEPTH EXTENT  (-100/-1500)

# creating mask file from bathymetry layer to clip to depth range

# Present

setwd("C:/Users/ez14/OneDrive - The University of Waikato/EDO/Chapter 1/Training_Fabrice/Present_Env_predictors_EEZ")
bathy.mask_present <- raster("Present_Bathymetry.tif")
bathy.mask_present[bathy.mask_present > -100 | bathy.mask_present <= -1500] <- NA
bathy.mask_present[bathy.mask_present <= -100 | bathy.mask_present > -1500] <- 1
plot(bathy.mask_present)

pred.PA_mean_present <- pred.PA_mean_present*bathy.mask_present
pred.AB_mean_present <- pred.AB_mean_present*bathy.mask_present
pred.PA_SD_present <- pred.PA_SD_present*bathy.mask_present
pred.AB_SD_present <- pred.AB_SD_present*bathy.mask_present

# SPATIAL PREDICTION: map of the mean PA & uncertainty of the PA
plot(pred.PA_mean_present)
plot(pred.PA_SD_present)

# SPATIAL PREDICTION: map of the mean AB (conditional on PA) & uncertainty of the AB
plot(pred.AB_mean_present)
plot(pred.AB_SD_present)

# SSP2

setwd("C:/Users/ez14/OneDrive - The University of Waikato/EDO/Chapter 1/Training_Fabrice/Future_Env_predictors_SSP2_EEZ")
bathy.mask_SSP2 <- raster("SSP2_Bathymetry.tif")
bathy.mask_SSP2[bathy.mask_SSP2 > -100 | bathy.mask_SSP2 <= -1500] <- NA
bathy.mask_SSP2[bathy.mask_SSP2 <= -100 | bathy.mask_SSP2 > -1500] <- 1
plot(bathy.mask_SSP2)

pred.PA_mean_SSP2 <- pred.PA_mean_SSP2*bathy.mask_SSP2
pred.AB_mean_SSP2 <- pred.AB_mean_SSP2*bathy.mask_SSP2
pred.PA_SD_SSP2 <- pred.PA_SD_SSP2*bathy.mask_SSP2
pred.AB_SD_SSP2 <- pred.AB_SD_SSP2*bathy.mask_SSP2

# SPATIAL PREDICTION: map of the mean PA & uncertainty of the PA
plot(pred.PA_mean_SSP2)
plot(pred.PA_SD_SSP2)

# SPATIAL PREDICTION: map of the mean AB (conditional on PA) & uncertainty of the AB
plot(pred.AB_mean_SSP2)
plot(pred.AB_SD_SSP2)

#SSP3

setwd("C:/Users/ez14/OneDrive - The University of Waikato/EDO/Chapter 1/Training_Fabrice/Future_Env_predictors_SSP3_EEZ")
bathy.mask_SSP3 <- raster("SSP3_Bathymetry.tif")
bathy.mask_SSP3[bathy.mask_SSP3 > -100 | bathy.mask_SSP3 <= -1500] <- NA
bathy.mask_SSP3[bathy.mask_SSP3 <= -100 | bathy.mask_SSP3 > -1500] <- 1
plot(bathy.mask_SSP3)

pred.PA_mean_SSP3 <- pred.PA_mean_SSP3*bathy.mask_SSP3
pred.AB_mean_SSP3 <- pred.AB_mean_SSP3*bathy.mask_SSP3
pred.PA_SD_SSP3 <- pred.PA_SD_SSP3*bathy.mask_SSP3
pred.AB_SD_SSP3 <- pred.AB_SD_SSP3*bathy.mask_SSP3

# SPATIAL PREDICTION: map of the mean PA & uncertainty of the PA
plot(pred.PA_mean_SSP3)
plot(pred.PA_SD_SSP3)

# SPATIAL PREDICTION: map of the mean AB (conditional on PA) & uncertainty of the AB
plot(pred.AB_mean_SSP3)
plot(pred.AB_SD_SSP3)

################################################################################

# SAVE OBJECTS
setwd("C:/Users/ez14/OneDrive - The University of Waikato/EDO/Chapter 1/Training_Fabrice/R_files/taxa.name_final")

# PRESENT
writeRaster(pred.PA_mean_present, "pred.PA_mean.taxa.name_present.RData", overwrite=TRUE)
writeRaster(pred.PA_SD_present, "pred.PA_SD.taxa.name_present.RData", overwrite=TRUE)
writeRaster(pred.AB_mean_present, "pred.AB_mean.taxa.name_present.RData", overwrite=TRUE)
writeRaster(pred.AB_SD_present, "pred.AB_SD.taxa.name_present.RData", overwrite=TRUE)
save(pred.PA_mean_present, file = "pred.PA_mean.taxa.name_present.RData")
save(pred.PA_SD_present, file = "pred.PA_SD.taxa.name_present.RData")
save(pred.AB_mean_present, file = "pred.AB_mean.taxa.name_present.RData")
save(pred.AB_SD_present, file = "pred.AB_SD.taxa.name_present.RData")

# SSP2
writeRaster(pred.PA_mean_SSP2, "pred.PA_mean.taxa.name_SSP2.RData", overwrite=TRUE)
writeRaster(pred.PA_SD_SSP2, "pred.PA_SD.taxa.name_SSP2.RData", overwrite=TRUE)
writeRaster(pred.AB_mean_SSP2, "pred.AB_mean.taxa.name_SSP2.RData", overwrite=TRUE)
writeRaster(pred.AB_SD_SSP2, "pred.AB_SD.taxa.name_SSP2.RData", overwrite=TRUE)
save(pred.PA_mean_SSP2, file = "pred.PA_mean.taxa.name_SSP2.RData")
save(pred.PA_SD_SSP2, file = "pred.PA_SD.taxa.name_SSP2.RData")
save(pred.AB_mean_SSP2, file = "pred.AB_mean.taxa.name_SSP2.RData")
save(pred.AB_SD_SSP2, file = "pred.AB_SD.taxa.name_SSP2.RData")

# SSP3
writeRaster(pred.PA_mean_SSP3, "pred.PA_mean.taxa.name_SSP3.RData", overwrite=TRUE)
writeRaster(pred.PA_SD_SSP3, "pred.PA_SD.taxa.name_SSP3.RData", overwrite=TRUE)
writeRaster(pred.AB_mean_SSP3, "pred.AB_mean.taxa.name_SSP3.RData", overwrite=TRUE)
writeRaster(pred.AB_SD_SSP3, "pred.AB_SD.taxa.name_SSP3.RData", overwrite=TRUE)
save(pred.PA_mean_SSP3, file = "pred.PA_mean.taxa.name_SSP3.RData")
save(pred.PA_SD_SSP3, file = "pred.PA_SD.taxa.name_SSP3.RData")
save(pred.AB_mean_SSP3, file = "pred.AB_mean.taxa.name_SSP3.RData")
save(pred.AB_SD_SSP3, file = "pred.AB_SD.taxa.name_SSP3.RData")

# FINAL FILES
writeRaster(pred.PA_mean_present, 'taxa.name.pred.PA_mean_present.tif', overwrite=TRUE)
writeRaster(pred.AB_mean_present, 'taxa.name.pred.AB_mean_present.tif', overwrite=TRUE)
writeRaster(pred.PA_mean_SSP2, 'taxa.name.pred.PA_mean_SSP2.tif', overwrite=TRUE)
writeRaster(pred.AB_mean_SSP2, 'taxa.name.pred.AB_mean_SSP2.tif', overwrite=TRUE)
writeRaster(pred.PA_mean_SSP3, 'taxa.name.pred.PA_mean_SSP3.tif', overwrite=TRUE)
writeRaster(pred.AB_mean_SSP3, 'taxa.name.pred.AB_mean_SSP3.tif', overwrite=TRUE)

writeRaster(pred.PA_SD_present, 'taxa.name.pred.PA_SD_present.tif', overwrite=TRUE)
writeRaster(pred.AB_SD_present, 'taxa.name.pred.AB_SD_present.tif', overwrite=TRUE)
writeRaster(pred.PA_SD_SSP2, 'taxa.name.pred.PA_SD_SSP2.tif', overwrite=TRUE)
writeRaster(pred.AB_SD_SSP2, 'taxa.name.pred.AB_SD_SSP2.tif', overwrite=TRUE)
writeRaster(pred.PA_SD_SSP3, 'taxa.name.pred.PA_SD_SSP3.tif', overwrite=TRUE)
writeRaster(pred.AB_SD_SSP3, 'taxa.name.pred.AB_SD_SSP3.tif', overwrite=TRUE)

###=============================================================================
###=============================================================================
###=============================================================================

# PRIMARY HABITATS

setwd("C:/Users/ez14/OneDrive - The University of Waikato/EDO/Chapter 1/Training_Fabrice/R_files/taxa.name_final")

##### taxa.name_Present
taxa.name.present <- raster("taxa.name.pred.AB_mean_present.tif")
plot(taxa.name.present, zlim = c(0,400000))
Q_pres <- quantile(taxa.name.present, na.rm=TRUE, probs = c(0.98))

# Primary habitats maps (occurrence at 98th percentile)
taxa.name.present.Q <- taxa.name.present
taxa.name.present.Q[values(taxa.name.present.Q) < Q_pres] <- 0
taxa.name.present.Q[values(taxa.name.present.Q) >= Q_pres] <- 1
plot(taxa.name.present.Q)
writeRaster(taxa.name.present.Q, "taxa.name.present.Q.tif", overwrite=TRUE)

# Primary habitats maps (density at 98th percentile)
taxa.name.present.Q.D <- taxa.name.present*taxa.name.present.Q
plot(taxa.name.present.Q.D)
writeRaster(taxa.name.present.Q.D, "taxa.name.present.Q.D.tif", overwrite=TRUE)
values(taxa.name.present.Q.D)[values(taxa.name.present.Q.D) == 0 ] <- NA

# stats
sum(values(taxa.name.present.Q.D), na.rm = T)
mean(values(taxa.name.present.Q.D), na.rm = T)
max(values(taxa.name.present.Q.D), na.rm = T)
min(values(taxa.name.present.Q.D), na.rm = T)

################################################################################

#####  taxa.name_Future (SSP2)
taxa.name.SSP2 <- raster("taxa.name.pred.AB_mean_SSP2.tif")
plot(taxa.name.SSP2, zlim = c(0,400000))

# Primary habitats maps (occurrence at 98th percentile)
taxa.name.SSP2.Q <- taxa.name.SSP2
taxa.name.SSP2.Q[values(taxa.name.SSP2.Q) < Q_pres ] <- 0
taxa.name.SSP2.Q[values(taxa.name.SSP2.Q) >= Q_pres ] <- 1
plot(taxa.name.SSP2.Q)
writeRaster(taxa.name.SSP2.Q, "taxa.name.SSP2.Q.tif", overwrite=TRUE)

# Primary habitats maps (density at 98th percentile)
taxa.name.SSP2.Q.D <- taxa.name.SSP2*taxa.name.SSP2.Q
plot(taxa.name.SSP2.Q.D)
writeRaster(taxa.name.SSP2.Q.D, "taxa.name.SSP2.Q.D.tif", overwrite=TRUE)
values(taxa.name.SSP2.Q.D)[values(taxa.name.SSP2.Q.D) == 0 ] <- NA

# stats
sum(values(taxa.name.SSP2.Q.D), na.rm = T)
mean(values(taxa.name.SSP2.Q.D), na.rm = T)
max(values(taxa.name.SSP2.Q.D), na.rm = T)
min(values(taxa.name.SSP2.Q.D), na.rm = T)

################################################################################

#####  taxa.name_Future (SSP3)
taxa.name.SSP3 <- raster("taxa.name.pred.AB_mean_SSP3.tif")
plot(taxa.name.SSP3, zlim = c(0,400000))

# Primary habitats maps (occurrence at 98th percentile)
taxa.name.SSP3.Q <- taxa.name.SSP3
taxa.name.SSP3.Q[values(taxa.name.SSP3.Q) < Q_pres ] <- 0
taxa.name.SSP3.Q[values(taxa.name.SSP3.Q) >= Q_pres ] <- 1
plot(taxa.name.SSP3.Q)
writeRaster(taxa.name.SSP3.Q, "taxa.name.SSP3.Q.tif", overwrite=TRUE)

# Primary habitats maps (density at 98th percentile)
taxa.name.SSP3.Q.D <- taxa.name.SSP3*taxa.name.SSP3.Q
plot(taxa.name.SSP3.Q.D)
writeRaster(taxa.name.SSP3.Q.D, "taxa.name.SSP3.Q.D.tif", overwrite=TRUE)
values(taxa.name.SSP3.Q.D)[values(taxa.name.SSP3.Q.D) == 0 ] <- NA

# stats
sum(values(taxa.name.SSP3.Q.D), na.rm = T)
mean(values(taxa.name.SSP3.Q.D), na.rm = T)
max(values(taxa.name.SSP3.Q.D), na.rm = T)
min(values(taxa.name.SSP3.Q.D), na.rm = T)

################################################################################

# Summarizing density stats

density_stats_matrix <- array(0, c(3,5))
colnames(density_stats_matrix) <- c("sum","mean","median","max","min")
row.names(density_stats_matrix) <- c("Present", "SSP2", "SSP3")

density_stats_matrix[1,1] <- sum(values(taxa.name.present.Q.D), na.rm = T)
density_stats_matrix[2,1] <- sum(values(taxa.name.SSP2.Q.D), na.rm = T)
density_stats_matrix[3,1] <- sum(values(taxa.name.SSP3.Q.D), na.rm = T)

density_stats_matrix[1,2] <- mean(values(taxa.name.present.Q.D), na.rm = T)
density_stats_matrix[2,2] <- mean(values(taxa.name.SSP2.Q.D), na.rm = T)
density_stats_matrix[3,2] <- mean(values(taxa.name.SSP3.Q.D), na.rm = T)

density_stats_matrix[1,3] <- median(values(taxa.name.present.Q.D), na.rm = T)
density_stats_matrix[2,3] <- median(values(taxa.name.SSP2.Q.D), na.rm = T)
density_stats_matrix[3,3] <- median(values(taxa.name.SSP3.Q.D), na.rm = T)

density_stats_matrix[1,4] <- max(values(taxa.name.present.Q.D), na.rm = T)
density_stats_matrix[2,4] <- max(values(taxa.name.SSP2.Q.D), na.rm = T)
density_stats_matrix[3,4] <- max(values(taxa.name.SSP3.Q.D), na.rm = T)

density_stats_matrix[1,5] <- min(values(taxa.name.present.Q.D), na.rm = T)
density_stats_matrix[2,5] <- min(values(taxa.name.SSP2.Q.D), na.rm = T)
density_stats_matrix[3,5] <- min(values(taxa.name.SSP3.Q.D), na.rm = T)

# FINAL RESULTS
print(density_stats_matrix)

###=============================================================================
###=============================================================================

# QUANTIFYING DENSITY CHANGE IN REGIONS OF PRIMARY HABITATS

# primary habitat overlap raster = 4 categories coded: -1=contraction; 1=no change, 2=expansion; 0=no relevant

# creating a mask to define the categories (NB: no need to create a mask layer for Present)

# SSP2
mask_taxa.name.SSP2 <- taxa.name.SSP2
mask_taxa.name.SSP2.Q <- mask_taxa.name.SSP2
mask_taxa.name.SSP2.Q[values(mask_taxa.name.SSP2.Q) < Q_pres ] <- 0
mask_taxa.name.SSP2.Q[values(mask_taxa.name.SSP2.Q) >= Q_pres ] <- 2
plot(mask_taxa.name.SSP2.Q)
writeRaster(mask_taxa.name.SSP2.Q, "mask_taxa.name.SSP2.Q.tif", overwrite=TRUE)

# SSP3
mask_taxa.name.SSP3 <- taxa.name.SSP3
mask_taxa.name.SSP3.Q <- mask_taxa.name.SSP3
mask_taxa.name.SSP3.Q[values(mask_taxa.name.SSP3.Q) < Q_pres ] <- 0
mask_taxa.name.SSP3.Q[values(mask_taxa.name.SSP3.Q) >= Q_pres ] <- 2

plot(mask_taxa.name.SSP3.Q)
writeRaster(mask_taxa.name.SSP3.Q, "mask_taxa.name.SSP3.Q.tif", overwrite=TRUE)

##==============================================================================
##==============================================================================

# CREATING FINAL LAYERS FOR PRESENT VS FUTURE SCENARIOS (subtracting future to present)

SSP2vsPres <- mask_taxa.name.SSP2.Q-taxa.name.present.Q
SSP3vsPres <- mask_taxa.name.SSP3.Q-taxa.name.present.Q
writeRaster(SSP2vsPres, "SSP2vsPres.tif", overwrite=TRUE)
writeRaster(SSP3vsPres, "SSP3vsPres.tif", overwrite=TRUE)

SSP2vsPres <- raster("SSP2vsPres.tif")
SSP3vsPres <- raster("SSP3vsPres.tif")

SSP2vsPres # check the object - there is no "values" section 
raster::plot(SSP2vsPres) # plot - it plots even without "values" which is weird
table(raster::values(SSP2vsPres)) # check if it can extract values - it can... so why isn't it shown in the object?
is.na(SSP2vsPres) # this now produces a "values" section in the object
SSP2vsPres[is.na(SSP2vsPres)] <- NA  # redefine NAs to make sure that they are compatible with R
SSP2vsPres # now the values section is shown in the object

plot(SSP2vsPres) # plot looks fine
plot(SSP3vsPres) # plot looks fine

################################################################################

# CATEGORY LOOP

# Creating the categories
cat <- c(-1,1,2)

# creating bootstrap matrix to save stats values
overlap_stats_matrix_SSP2 <- array(0, c(6,3))
overlap_stats_matrix_SSP3 <- array(0, c(6,3))
overlap_stats_matrix_SSP2.1 <- array(0, c(6,3))
overlap_stats_matrix_SSP3.1 <- array(0, c(6,3))

colnames(overlap_stats_matrix_SSP2) <- c("contraction","no_change_pres","expansion")
row.names(overlap_stats_matrix_SSP2) <- c("extent","mean","median","max","min","sum")
colnames(overlap_stats_matrix_SSP3) <- c("contraction","no_change_pres","expansion")
row.names(overlap_stats_matrix_SSP3) <- c("extent","mean","median","max","min","sum")
colnames(overlap_stats_matrix_SSP2.1) <- c("contraction","no_change_future","expansion")
row.names(overlap_stats_matrix_SSP2.1) <- c("extent","mean","median","max","min","sum")
colnames(overlap_stats_matrix_SSP3.1) <- c("contraction","no_change_future","expansion")
row.names(overlap_stats_matrix_SSP3.1) <- c("extent","mean","median","max","min","sum")

i=1

for (i in 1:length(cat)){
  
  # SSP2 vs Present
  
  mask_SSP2vsPres <- SSP2vsPres
  mask_SSP2vsPres[mask_SSP2vsPres == cat[i]] <- 10
  mask_SSP2vsPres[mask_SSP2vsPres <10] <- 0
  mask_SSP2vsPres[mask_SSP2vsPres ==10] <- 1
  plot(mask_SSP2vsPres)
  
  overlap_stats_matrix_SSP2[1,i] <- sum(na.omit(values(mask_SSP2vsPres)))
  
  temp <- taxa.name.present.Q.D * mask_SSP2vsPres
  temp[temp == 0] <- NA
  
  overlap_stats_matrix_SSP2[2,i] <- mean(values(temp), na.rm = T)
  overlap_stats_matrix_SSP2[3,i] <- median(values(temp), na.rm = T)
  overlap_stats_matrix_SSP2[4,i] <- max(values(temp), na.rm = T)
  overlap_stats_matrix_SSP2[5,i] <- min(values(temp), na.rm = T)
  overlap_stats_matrix_SSP2[6,i] <- sum(values(temp), na.rm = T)
  
  temp <- taxa.name.SSP2.Q.D * mask_SSP2vsPres
  temp[temp == 0] <- NA
  
  overlap_stats_matrix_SSP2.1[1,i] <- sum(na.omit(values(mask_SSP2vsPres)))
  overlap_stats_matrix_SSP2.1[2,i] <- mean(values(temp), na.rm = T)
  overlap_stats_matrix_SSP2.1[3,i] <- median(values(temp), na.rm = T)
  overlap_stats_matrix_SSP2.1[4,i] <- max(values(temp), na.rm = T)
  overlap_stats_matrix_SSP2.1[5,i] <- min(values(temp), na.rm = T)
  overlap_stats_matrix_SSP2.1[6,i] <- sum(values(temp), na.rm = T)
  
  # SSP3 vs Present
  
  mask_SSP3vsPres <- SSP3vsPres
  mask_SSP3vsPres[mask_SSP3vsPres == cat[i]] <- 10
  mask_SSP3vsPres[mask_SSP3vsPres <10] <- 0
  mask_SSP3vsPres[mask_SSP3vsPres ==10] <- 1
  plot(mask_SSP3vsPres)
  
  overlap_stats_matrix_SSP3[1,i] <- sum(na.omit(values(mask_SSP3vsPres)))
  
  temp.1 <- taxa.name.present.Q.D * mask_SSP3vsPres
  temp.1[temp.1 == 0] <- NA
  
  overlap_stats_matrix_SSP3[2,i] <- mean(values(temp.1), na.rm = T)
  overlap_stats_matrix_SSP3[3,i] <- median(values(temp.1), na.rm = T)
  overlap_stats_matrix_SSP3[4,i] <- max(values(temp.1), na.rm = T)
  overlap_stats_matrix_SSP3[5,i] <- min(values(temp.1), na.rm = T)
  overlap_stats_matrix_SSP3[6,i] <- sum(values(temp.1), na.rm = T)
  
  temp.1 <- taxa.name.SSP3.Q.D * mask_SSP3vsPres
  temp.1[temp.1 == 0] <- NA
  
  overlap_stats_matrix_SSP3.1[1,i] <- sum(na.omit(values(mask_SSP3vsPres)))
  overlap_stats_matrix_SSP3.1[2,i] <- mean(values(temp.1), na.rm = T)
  overlap_stats_matrix_SSP3.1[3,i] <- median(values(temp.1), na.rm = T)
  overlap_stats_matrix_SSP3.1[4,i] <- max(values(temp.1), na.rm = T)
  overlap_stats_matrix_SSP3.1[5,i] <- min(values(temp.1), na.rm = T)
  overlap_stats_matrix_SSP3.1[6,i] <- sum(values(temp.1), na.rm = T)
  
}

# Summarizing density stats
print(overlap_stats_matrix_SSP2)
print(overlap_stats_matrix_SSP2.1)
print(overlap_stats_matrix_SSP3)
print(overlap_stats_matrix_SSP3.1)

# SSP2 

future_ssp2 <- as.data.frame(overlap_stats_matrix_SSP2)
future_ssp2.1 <- as.data.frame(overlap_stats_matrix_SSP2.1)
future_ssp2 <- future_ssp2[, c(1, 2)]
future_ssp2.1 <- future_ssp2.1[, c(2, 3)]

new_overlap_stats_matrix_SSP2 <- array(0, c(6,4))
colnames(new_overlap_stats_matrix_SSP2) <- c("contraction","no_change_pres","no_change_future","expansion")
row.names(new_overlap_stats_matrix_SSP2) <- c("extent","mean","median","max","min","sum")

new_overlap_stats_matrix_SSP2[,1] <- future_ssp2$contraction
new_overlap_stats_matrix_SSP2[,2] <- future_ssp2$no_change_pres
new_overlap_stats_matrix_SSP2[,3] <- future_ssp2.1$no_change_future
new_overlap_stats_matrix_SSP2[,4] <- future_ssp2.1$expansion


# SSP3 

future_SSP3 <- as.data.frame(overlap_stats_matrix_SSP3)
future_SSP3.1 <- as.data.frame(overlap_stats_matrix_SSP3.1)
future_SSP3 <- future_SSP3[, c(1, 2)]
future_SSP3.1 <- future_SSP3.1[, c(2, 3)]

new_overlap_stats_matrix_SSP3 <- array(0, c(6,4))
colnames(new_overlap_stats_matrix_SSP3) <- c("contraction","no_change_pres","no_change_future","expansion")
row.names(new_overlap_stats_matrix_SSP3) <- c("extent","mean","median","max","min","sum")

new_overlap_stats_matrix_SSP3[,1] <- future_SSP3$contraction
new_overlap_stats_matrix_SSP3[,2] <- future_SSP3$no_change_pres
new_overlap_stats_matrix_SSP3[,3] <- future_SSP3.1$no_change_future
new_overlap_stats_matrix_SSP3[,4] <- future_SSP3.1$expansion

# FINAL RESULTS
print(new_overlap_stats_matrix_SSP2)
print(new_overlap_stats_matrix_SSP3)

###=============================================================================
###=============================================================================

# DENSITY WITHIN EXPANSION AREAS

density.exp_stats_matrix <- array(0, c(5,2))
colnames(density.exp_stats_matrix) <- c("SSP2","SSP3")
row.names(density.exp_stats_matrix) <- c("mean","median","max","min","sum")

# SSP2

mask_SSP2vsPres <- SSP2vsPres
mask_SSP2vsPres[mask_SSP2vsPres == cat[3]] <- 10
mask_SSP2vsPres[mask_SSP2vsPres <10] <- 0
mask_SSP2vsPres[mask_SSP2vsPres ==10] <- 1
plot(mask_SSP2vsPres)

temp <- taxa.name.present * mask_SSP2vsPres
temp[temp == 0] <- NA

# stats (mean; max; min; sum)
density.exp_stats_matrix[1,1] <- mean(na.omit(values(temp)))
density.exp_stats_matrix[2,1] <- median(na.omit(values(temp)))
density.exp_stats_matrix[3,1] <- max(na.omit(values(temp)))
density.exp_stats_matrix[4,1] <- min(na.omit(values(temp)))
density.exp_stats_matrix[5,1] <- sum(na.omit(values(temp)))

###=============================================================================

# SSP3

mask_SSP3vsPres <- SSP3vsPres
mask_SSP3vsPres[mask_SSP3vsPres == cat[3]] <- 10
mask_SSP3vsPres[mask_SSP3vsPres <10] <- 0
mask_SSP3vsPres[mask_SSP3vsPres ==10] <- 1
plot(mask_SSP3vsPres)

temp.1 <- taxa.name.present * mask_SSP3vsPres
temp.1[temp.1 == 0] <- NA

# stats (mean; max; min; sum)

density.exp_stats_matrix[1,2] <- mean(na.omit(values(temp.1)))
density.exp_stats_matrix[2,2] <- median(na.omit(values(temp.1)))
density.exp_stats_matrix[3,2] <- max(na.omit(values(temp.1)))
density.exp_stats_matrix[4,2] <- min(na.omit(values(temp.1)))
density.exp_stats_matrix[5,2] <- sum(na.omit(values(temp.1)))

# FINAL RESULTS
print(density.exp_stats_matrix)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# RECAP FINAL RESULTS

# MODEL PERFORMANCE (AUC, p)
round(mean_AUC,2)
round(SD_AUC,2)
round(mean_p,2)
round(SD_p,2)

# QUANTIFYNG OVERALL DENSITY CHANGE OVER TIME
round(sum(values(taxa.name.present), na.rm = T))
round(sum(values(taxa.name.SSP2), na.rm = T))
round(sum(values(taxa.name.SSP3), na.rm = T))

# QUANTIFYNG PRIMARY DENSITY HABIAT CHANGE OVER TIME 
round(density_stats_matrix)

# QUANTIFYING PRIMARY DENSITY HABIAT LOSS AND GAIN OVER TIME
round(new_overlap_stats_matrix_SSP2)
round(new_overlap_stats_matrix_SSP3)

# QUANTIFYING PRESENT DENSITY IN FUTURE EXPANSION AREAS
round(density.exp_stats_matrix)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


save.image("C:/Users/ez14/OneDrive - The University of Waikato/EDO/Chapter 1/Training_Fabrice/R_files/taxa.name_final/1.SCRIPT")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# PARTIAL DEPENDENCE PLOTS (PDPs)

# PROCESS
# create data frame for each important variable and put in a list (dataframe, colnames = the variable names)

Aragonite <- as.data.frame(taxa.name.F$Aragonite)
colnames(Aragonite) <- "Aragonite"
Bathy <- as.data.frame(taxa.name.F$Bathy)
colnames(Bathy) <- "Bathy"
Bot_sal <- as.data.frame(taxa.name.F$Bot_sal)
colnames(Bot_sal) <- "Bot_sal"
Bpi_broad <- as.data.frame(taxa.name.F$Bpi_broad)
colnames(Bpi_broad) <- "Bpi_broad"
Bpi_fine <- as.data.frame(taxa.name.F$Bpi_fine)
colnames(Bpi_fine) <- "Bpi_fine"
Calcite <- as.data.frame(taxa.name.F$Calcite)
colnames(Calcite) <- "Calcite"
Detr_flux <- as.data.frame(taxa.name.F$Detr_flux)
colnames(Detr_flux) <- "Detr_flux"
Nitro <- as.data.frame(taxa.name.F$Nitro)
colnames(Nitro) <- "Nitro"
OXY_C <- as.data.frame(taxa.name.F$OXY_C)
colnames(OXY_C) <- "OXY_C"
Prof <- as.data.frame(taxa.name.F$Prof)
colnames(Prof) <- "Prof"
Slope <- as.data.frame(taxa.name.F$Slope)
colnames(Slope) <- "Slope"

pdp.list <- as.data.frame(list(Aragonite,Bathy,Bpi_broad,Bpi_fine,Detr_flux,Nitro,OXY_C,Prof,Slope))

# run bootstraps
set.seed(2022)
n.boot <- 100
n.var <- length(imp.var)

# create bootstrap matrix for each variable (and for each model - M1 and M2)
Bootsrap.matrix_var <- array(0,c(nrow(taxa.name.F),length(imp.var),n.boot))
dimnames(Bootsrap.matrix_var)[[2]] <- imp.var
dimnames(Bootsrap.matrix_var)[[3]] <- paste('Rep_',seq(1,n.boot),sep="")

Bootsrap.matrix_var2 <- array(0,c(nrow(taxa.name.F),length(imp.var),n.boot))
dimnames(Bootsrap.matrix_var2)[[2]] <- imp.var
dimnames(Bootsrap.matrix_var2)[[3]] <- paste('Rep_',seq(1,n.boot),sep="")

i = 1; j = 1
for (i in 1:n.boot){
  
  train_taxa.name.F.raw <- sample(seq_len(nrow(taxa.name.F)), size = nrow(taxa.name.F), replace = T) # index of rows for presence data
  train_taxa.name.F <- taxa.name.F[train_taxa.name.F.raw, ] #train for PA
  train_taxa.name.F_1 <- na.omit(train_taxa.name.F ) #train for Abundance
  
  # run both models
  M1 <- randomForest(x = train_taxa.name.F [,imp.var],
                     y = train_taxa.name.F$taxa.name.PA,
                     mtryStart = 2,
                     ntreeTry = 2000,
                     stepFactor = 2,
                     improve = 0.0005,
                     trace = T, plot = T, doBest = T) # for PA
  
  M2 <- randomForest(x = train_taxa.name.F_1[,imp.var],
                     y = train_taxa.name.F_1$taxa.name.D,
                     mtryStart = 2,
                     ntreeTry = 1000,
                     stepFactor = 2,
                     improve = 0.005,
                     trace = T, plot = T, doBest = T) # for Abundance
  
  # Nested loop for your variables (one for each model)
  
  # pd.1 for M1 (PA)
  # pd.2 for M2 (AB)
  
  for (j in 1:n.var){ 
    
    pd.1 <- pdp::partial(M1, pred.grid = pdp.list[j], pred.var = c(imp.var[j]), chull = TRUE, prob = TRUE)
    pd.2 <- pdp::partial(M2, pred.grid = pdp.list[j], pred.var = c(imp.var[j]), chull = TRUE, prob = F)
    
    Bootsrap.matrix_var[,j,i] <- pd.1[,2]
    Bootsrap.matrix_var2[,j,i] <- pd.1[,2]*exp(pd.2[,2])
    
  }
}

Bootsrap.matrix_var[,2,]

################################################################################

# plot PDs + 95 PI and save to file
library(devEMF)
emf(file = "taxa.name.imp.var.PA_Present.emf", emfPlus = T)
par(mar=c(4, 2, 1, 1))
par(mfrow=c(4,3))

# PRESENCE-ABSENCE

# change this value to fit all Present range or change to have order of var importance

for (i in c(1:length(imp.var))) {
  plot(sort(pdp.list[,i]),apply(Bootsrap.matrix_var[,i,],1,mean), col = "black",type='l',
       xlab = paste(imp.var[i], " (", round(mean_imp.var_PA[i],1), " % ± ", round(SD_imp.var_PA[i]), ")", sep=""), 
       ylab = '',
       ylim = c(0, 1)) #max(boot_array_EnvTran[,,])))
  # 95% PI
  UC <- na.omit(cbind(sort(pdp.list[,imp.var[i]]),
                      apply(Bootsrap.matrix_var[,i,],1, sd),
                      apply(Bootsrap.matrix_var[,i,],1, sd)))
  UC[,2] <- apply(Bootsrap.matrix_var[,i,],1,mean) + UC[,2]
  UC[,3] <- apply(Bootsrap.matrix_var[,i,],1,mean) - UC[,3]
  polygon(c(UC[,1], rev(UC[,1])),c(UC[,2], rev(UC[,3])), col = rgb(0,0,0, 0.25), border = NA)
  # 95% CI
  # lines(EnvRanges.DF[,i],apply(boot_array_EnvTran.DF[,i,],1,mean) + 2 * (sqrt(apply(boot_array_EnvTran.DF[,i,],1,var))/sqrt(5)), lty = 'dashed')
  # lines(EnvRanges.DF[,i],apply(boot_array_EnvTran.DF[,i,],1,mean) - 2 * (sqrt(apply(boot_array_EnvTran.DF[,i,],1,var))/sqrt(5)), lty = 'dashed')
  # lines(EnvRanges[,i],apply(boot_array_EnvTran[,i,],1,mean),col = 2)
  rug(quantile(taxa.name.F[,imp.var[i]],seq(0,1,0.1), na.rm = T), ticksize = 0.05, side = 1, lwd = 0.75)
  # if (i == 1) title('Sample size - 1000')
} 
dev.off()

emf(file = "taxa.name.imp.var.AB_Present.emf", emfPlus = T)
par(mar=c(4, 2, 1, 1))
par(mfrow=c(4,3))

# ABUNDANCE 

for (i in c(1:length(imp.var))) {
  # 95% PI
  UC <- na.omit(cbind(sort(pdp.list[,imp.var[i]]),
                      apply(Bootsrap.matrix_var2[,i,],1, sd),
                      apply(Bootsrap.matrix_var2[,i,],1, sd)))
  UC[,2] <- apply(Bootsrap.matrix_var2[,i,],1,mean) + UC[,2]
  UC[,3] <- apply(Bootsrap.matrix_var2[,i,],1,mean) - UC[,3]
  
  plot(sort(pdp.list[,i]),apply(Bootsrap.matrix_var2[,i,],1,mean), col = "black",type='l',
       xlab = paste(imp.var[i], " (", (round(mean_imp.var_PA[i],1) + round(mean_imp.var_AB[i],1))/2, " % ± ", (round(SD_imp.var_PA[i],1) + round(SD_imp.var_AB[i],1))/2, ")", sep=""),  
       ylab = '',
       ylim = c(min(UC[,3]), max(UC[,2]))) #max(boot_array_EnvTran[,,])))
  
  polygon(c(UC[,1], rev(UC[,1])),c(UC[,2], rev(UC[,3])), col = rgb(0,0,0, 0.25), border = NA)
  # 95% CI
  # lines(EnvRanges.DF[,i],apply(boot_array_EnvTran.DF[,i,],1,mean) + 2 * (sqrt(apply(boot_array_EnvTran.DF[,i,],1,var))/sqrt(5)), lty = 'dashed')
  # lines(EnvRanges.DF[,i],apply(boot_array_EnvTran.DF[,i,],1,mean) - 2 * (sqrt(apply(boot_array_EnvTran.DF[,i,],1,var))/sqrt(5)), lty = 'dashed')
  # lines(EnvRanges[,i],apply(boot_array_EnvTran[,i,],1,mean),col = 2)
  rug(quantile(taxa.name.F[,imp.var[i]],seq(0,1,0.1), na.rm = T), ticksize = 0.05, side = 1, lwd = 0.75)
  # if (i == 1) title('Sample size - 1000')
} 
dev.off()

###=============================================================================
###=============================================================================

# ENV COVERAGE

####==================    LOAD FILES AND PACKAGES         ==========================####
require(raster); require(pROC); require(dismo)

# load env data
setwd("C:/Users/ez14/OneDrive - The University of Waikato/EDO/Chapter 1/Training_Fabrice/R_files/GD_final")
load("Pred_1km_Present.RData")
load("Pred_1km_SSP2.RData")
load("Pred_1km_SSP3.RData")

#load one of your rasters which you will use a various points
setwd("C:/Users/ez14/OneDrive - The University of Waikato/EDO/Chapter 1/Training_Fabrice/Present_Env_predictors_EEZ")
R <- raster("Present_bot_sal_EEZ.tif")
R[R>0] <- 1 # turn this into a mask
plot(R)

# load biological data - ALL THE DTIS
DF <- Bio.D

####==================    PREPARE THE DATA FOR MODELLING          ==========================####
# assign the unique spatial identifier from your raster mask to all DTIS locations - this is to 
# make sure when you randomly sample the EEZ for your absences you don't pick the same place as your 
# DTIS location
DF$FID <- cellFromXY(R, DF[,c("X","Y")])
# set all your DTIS locations to "present" - i.e., you know have sampled that environmental condition 
# and at some point used it in your models.
DF$sample.sites <- 1

# ABSENCES
Pred_1km_U <-  Pred_1km_Present
Pred_1km_U$FID <-  cellFromXY(R, as.matrix(cbind(Pred_1km_Present$X, Pred_1km_Present$Y)))
Pred_1km_U <- na.omit(Pred_1km_Present[!Pred_1km_U$FID %in% DF$FID,]) # no overlap with presences
set.seed(5)
# sampling randomly from the full environemtnal space (minus the DTIS locations; with 4 times the absences as presences)
train_ind <- sample(seq_len(nrow(Pred_1km_U)), size = nrow(DF)*4)
preddat_EC <- Pred_1km_U[train_ind, ]
preddat_EC <- preddat_EC[,-c(1:2)] # remove the x y coords - don't need them - just keep the environemtnal variable values

# set the sampled absences as 0 (using the same col name as for the presences)
preddat_EC$sample.sites <- 0 # absent

#extract env info (from our predictorStacks - e.g., Bathy or OXY)
# Present_env <- as.data.frame(raster::extract(x=PredStack1km_Present, y=xy, method='simple', buffer=NULL, small=FALSE, cellnumbers=FALSE))
# colnames(Present_env) <- c("Aragonite","Bathy","Bot_sal","Bpi_broad","Bpi_fine",
#                            "Calcite","Detr_flux","Nitro","OXY_C","Prof","Slope")

# bind DTIS biodata with env variables
DF <- cbind(DF,Present_env)
DF_TRY <- na.omit(DF[,colnames(preddat_EC)])

# bind by row the DTIS data and the randomly sampled absences - make sure you have the environemtnal data from your models 
# and that the names for these in the collumns of the dataframes are the same.
DF_ES <- rbind(preddat_EC,DF_TRY) # only keep the 'sample.sites' colum and the environemtnal variables column
DF_ES$sample.sites <- as.factor(DF_ES$sample.sites) # turn the P/A to factor for RF analysis

# all the environmental variables you used in your models
all.var <- c("Aragonite","Bathy","Bot_sal","Bpi_broad","Bpi_fine",
             "Calcite","Detr_flux","Nitro","OXY_C","Prof","Slope") 

####==================    MODELLING AND PREDICTIONS         ==========================####
# RANDOM FOREST MODEL FOR YOUR Presence/absence ~ environmental variables
# make sure this a really big model - you're going to run it once so need to use lots of trees etc. 
M1 <- randomForest(x = DF_ES[,all.var],
                   y = DF_ES$sample.sites,
                   mtryStart = 2,
                   ntreeTry = 2000,
                   stepFactor = 2,
                   improve = 0.0005,
                   trace = T, plot = T, doBest = T) # for PA

# Check AUC to make sure it's ok (> 0.8 should be fine)
preds <- predict(M1, type="prob")[,2]
pROC::roc(DF_ES$sample.sites, preds, quiet = T)$auc

#Present day prediction
Env.Cov.Present <- as.numeric(round(predict(M1, Pred_1km_Present, type = "prob")[,2], digits = 2))

# future prediction
Env.Cov.SSP2 <- as.numeric(round(predict(M1, Pred_1km_SSP2, type = "prob")[,2], digits = 2))
Env.Cov.SSP3 <- as.numeric(round(predict(M1, Pred_1km_SSP3, type = "prob")[,2], digits = 2))

# turn each of these predictions (this is just for the present day env coverage) into a raster,plot them and save them 
t_PRESENT <- rasterFromXYZ(data.frame(x = Pred_1km_Present[,1],
                              y = Pred_1km_Present[,2],
                              z = Env.Cov.Present),
                              crs = crs(R)) # use whatever projection use have for your models
plot(t_PRESENT)
writeRaster(t_PRESENT,filename = "Env_cov_PRESENT.tif", overwrite=T) # write raster file

#SSP2
t_SSP2 <- rasterFromXYZ(data.frame(x = Pred_1km_SSP2[,1],
                              y = Pred_1km_SSP2[,2],
                              z = Env.Cov.SSP2),
                              crs = crs(R)) # use whatever projection use have for your models
plot(t_SSP2)


#SSP3
t_SSP3 <- rasterFromXYZ(data.frame(x = Pred_1km_SSP3[,1],
                              y = Pred_1km_SSP3[,2],
                              z = Env.Cov.SSP3),
                             crs = crs(R)) # use whatever projection use have for your models
plot(t_SSP3)

# SAVE

setwd("C:/Users/ez14/OneDrive - The University of Waikato/EDO/Chapter 1/Training_Fabrice/R_files/2_Env_Coverage")
writeRaster(t_PRESENT,filename = "Env_cov_PRESENT.tif", overwrite=T) # write raster file
writeRaster(t_SSP2,filename = "Env_cov_SSP2.tif", overwrite=T) # write raster file
writeRaster(t_SSP3,filename = "Env_cov_SSP3.tif", overwrite=T) # write raster file
