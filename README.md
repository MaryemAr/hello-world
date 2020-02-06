Ge# hello-world
R code for Digital soil mapping models

### Statistical model

# Packages
install.packages("ithir")
install.packages("raster")
install.packages("rgdal")
install.packages("sp")
install.packages("rpart")
install.packages("gstat")
install.packages("Cubist")
install.packages("randomForest")
install.packages("e1071")

## libraries
library(ithir)
library(raster)
library(rgdal)
library(sp)
library(rpart)
library(gstat)
library(Cubist)
library(randomForest)
library(e1071)

## set working directory
setwd("H:/2019_Manuscripts/2019_Cubist_RK/After_US_Analysis/Cb_RF_SVM/Analysis/Tom_sampling/Analysis/Gamma_EM/SVM_RK")

## input data
cDat <- read.table("Pro_soil_cali.txt", header=T, sep = ",") #calibration dataset
vDat <- read.table("Pro_soil_vali.txt", header=T, sep = ",") #validation dataset

## Selection of digital data using forward stepwise linear regression
## fit the model to calibration datasets using full set of covariates
soil.c <- lm(Ca_30 ~Pcon_1m+Hcon_1m+Pcon_2m+Hcon_2m ,data = cDat)

## Select the best covariates. Remove the variables with less p value one by one
soil.c.c <- stepAIC(soil.c, direction = "forward", 
                    trace = FALSE)
summary(soil.c.c)

## Fit the SVM model with best selcted covariates
tuneResult <- tune(svm, Mg_30 ~TC+Pcon_2m+K, data = cDat, ranges = list(epsilon = seq(0.1,0.2,0.02), cost = c(5,7,15,20)))
soil.c.c <- tuneResult$best.model

# SVM model only validation
SVM.pred.v <- predict(soil.c.c, newdata = vDat)
goof1<- goof(observed = vDat$Mg_30,predicted = SVM.pred.v)
print(goof1)
write.table(goof1, "goof_SVM_Gamma_EM_Mg_30.txt", sep = ",")
write.table(SVM.pred.v, "Measured_pred_SVM_Gamma_EM_Mg_30.txt", sep = ",")

## Derive regression residual (RR)
cDat$residual <- cDat$Mg_30 - predict(soil.c.c, cDat)
mean(cDat$residual)

## Convert data to geodata
coordinates(cDat) <- ~x + y

## Fit variogram
vgm1 <- variogram(residual ~ 1, cDat)
mod <- vgm(psill = var(cDat$residual), "Sph", range = 450,
           nugget = 0.5)
model_1 <- fit.variogram(vgm1, mod)
model_1

# Residual kriging model
gRK <- gstat(NULL, "RKresidual", residual ~ 1, cDat,
             model = model_1)


## SVM model with residual variogram using ordinary kriging and vDat
coordinates(vDat) <- ~x + y
RK.preds.V <- as.data.frame(krige(residual ~ 1, cDat, model = model_1,
                                  newdata = vDat))

## Sum the two components together
RK.preds.fin2 <- SVM.pred.v + RK.preds.V[, 3]


##  SVMRR model validation
goof2<- goof(observed = vDat$Mg_30,
     predicted = RK.preds.fin2)
print(goof2)

write.table(goof2, "goof_RK_SVM_Gamma_EM_Mg_30.txt", sep = ",")
write.table(RK.preds.fin2, "Measured_pred_RK_SVM_Gamma_EM_Mg_30.txt", sep = ",")

## Prediction onto grid (where you only have digital data)
grid <- read.table("Grid.txt", header=T, sep = ",")
lm.pred.g <- predict(soil.c.c, newdata = grid)
coordinates(grid) <- ~x + y
RK.preds.g <- as.data.frame(krige(residual ~ 1, cDat, model = model_1,
                                  newdata = grid))

## Sum the two components together
RK.preds.fin <- lm.pred.g + RK.preds.g[, 3]
error <- RK.preds.g$var1.pred+RK.preds.g$var1.var


write.table(RK.preds.fin, "Grid_Pred_LM_RK_Mg_30.txt", sep = ",")
write.table(error, "Error_Pred_LM_RK_Mg_30.txt", sep = ",")



