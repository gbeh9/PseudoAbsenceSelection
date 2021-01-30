## ================================================================================================================== ##
##                                                                                                                    ##
## Script name:   Applying Barbet-Massin et al.'s pseudo-absences selection in SDMs                                   ##
##                                                                                                                    ##
## Purpose of script:                                                                                                 ##
##                                                                                                                    ##
## Author: MSc. Luíz Fernando Esser                                                                                   ##
##                                                                                                                    ##
## Date Created: 2021-01-26                                                                                           ##
##                                                                                                                    ##
## This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.   ##
## To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/                            ##
## or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.                                   ##
##                                                                                                                    ##
## Email: luizesser@gmail.com                                                                                         ##
##                                                                                                                    ##
## ================================================================================================================== ##
##                                                                                                                    ##
## Notes:                                                                                                             ##
##                                                                                                                    ##
##                                                                                                                    ##
## ================================================================================================================== ##
library(biomod2)
library(dplyr)
library(ecospat)
library(sdm)
library(rgdal)
library(raster)

set.seed(1)

### Read the species data:
sp_dat <- shapefile(system.file("external/po_spatial_points.shp", package="sdm"))
sp_dat2 <- sp_dat@data

### Make a stack of predictors:
current <- stack(system.file("external/predictors.grd", package="sdm"))

### Generate Backgrund Data:
# Generate Background for Group 1: GLM, GAM (random, PA = 1000)
d1 <- sdmData(reformulate(termlabels = names(current), response = names(sp_dat)[1]), 
              train = sp_dat, predictors = current,
              bg=list(n=1000,method='eRandom',remove=TRUE))

# Generate Background for Group 2: MARS (random, PA = 100)
d2 <- sdmData(reformulate(termlabels = names(current), response = names(sp_dat)[1]), 
              train = sp_dat, predictors = current,
              bg=list(n=100,method='eRandom',remove=TRUE))

# 
# Generate Background for Group 3: MDA (SRE, PA = 100)
sp_format <- BIOMOD_FormatingData(resp.var = rep( 1, nrow(sp_dat2) ),
                                  expl.var = current,
                                  resp.xy = sp_dat2[,c("coords_x1","coords_x2")], # colnames from coordinates in sp_dat2
                                  resp.name = "sp4", # species name in sp_dat2
                                  PA.strategy = "sre",
                                  PA.nb.rep = 1,
                                  PA.nb.absences = 100)
bg <- cbind(sp = c(rep(1,nrow(sp_dat2)),rep(0, 100)),
            sp_format@coord)
names(bg)[1] <- names(sp_dat2)[1]
coordinates(bg) <- 2:3
d3 <- sdmData(reformulate(termlabels = names(current), response = names(bg)[1]),
              train = bg, predictors = current)

# Generate Background for Group 4: CTA, BRT, RF (2-degree far, PA = presence)
print("Preparing sdmData for G4 (SRE, PA = presence)") 
sp_format <- BIOMOD_FormatingData(resp.var = rep( 1, nrow(sp_dat2) ),
                                  expl.var = current,
                                  resp.xy = sp_dat2[,c("coords_x1","coords_x2")],
                                  resp.name = "sp4",
                                  PA.strategy = "disk",
                                  PA.dist.min = 220000, # this is approximately what 2 degree is in meters.
                                  PA.nb.rep = 1,
                                  PA.nb.absences = nrow(sp_dat2))
bg <- cbind(sp = c(rep(1,nrow(sp_dat2)),rep(0, nrow(sp_dat2))),
            sp_format@coord)
names(bg)[1] <- names(sp_dat2)[1]
coordinates(bg) <- 2:3
d4 <- sdmData(reformulate(termlabels = names(current), response = names(bg)[1]),
              train = bg, predictors = current)

### Generate models
nruns <- 1
folds <- 4
print("Running Group 1...")
m1 <- sdm(reformulate(termlabels = names(current)), 
          data=d1, methods = c('gam','glm'),
          replicate.method = 'cv', cv.folds=folds, n=nruns)

print("Running Group 2...")
m2 <- sdm(reformulate(termlabels = names(current)), 
          data=d2, methods = c('mars'), 
          replicate.method = 'cv', cv.folds=folds, n=nruns)

print("Running Group 3...")
m3 <- sdm(reformulate(termlabels = names(current)), 
          data=d3, methods = c('mda'), 
          replicate.method = 'cv', cv.folds=folds, n=nruns)

print("Running Group 4...")
m4 <- sdm(reformulate(termlabels = names(current)), 
          data=d4, methods = c('brt', 'rf', 'rpart'),
          replicate.method = 'cv', cv.folds=folds, n=nruns)

### Evaluate Models
print("Evaluating AUC and TSS...")
eval <- lapply(c("m1","m2","m3","m4"), 
               function(x) { getEvaluation(get(x),wtest='training',stat=c('AUC','TSS'),opt=1) } )
n <- nruns*folds # number of models in each algorithm
algo <- c('gam','glm','mars','mda','brt', 'rf', 'rpart')
eval_algo <- data.frame(algo = c(rep(algo, each=n)),bind_rows(eval))
good_models <- eval_algo[eval_algo[,"AUC"] >= 0.8 & eval_algo[,"TSS"] >= 0.5,]
good_models

### Project Models
dir.create(m1@data@species.names)
ensemble_by_algo <- function(m){
  m <- get(m)
  algo_m <- as.character(unique(m@run.info$method))
  for (i in 1:length(algo_m)){
    algo_gm <- good_models[good_models[,"algo"] == algo_m[i],]
    if(nrow(algo_gm) == 0){print(paste0("No good models...")) } 
    if(nrow(algo_gm) == 1){print(paste0("One good model, projecting..."))
      sp <- m@data@species.names
      e <- predict(m,
                   newdata=current,
                   filename=paste0(sp,"/",sp,"_",algo_m[i],".img"),
                   w=algo_gm$modelID)
    }
    if(nrow(algo_gm) > 1){print(paste0("Multiple good models, ensembling..."))
      sp <- m@data@species.names
      e <- ensemble(m,
                    newdata=current,
                    filename=paste0(sp,"/",sp,"_",algo_m[i],".img"),
                    id=algo_gm$modelID,
                    setting=list(method='weighted'))
    }
  }
}
x <- lapply(c("m1","m2","m3","m4"), ensemble_by_algo)

sp <- m1@data@species.names
for(s in sp){
  print(s)
  assign(s, stack(list.files(s, pattern = ".img$", rec = T, full.names = T)))
}

result <- sum(sp4)

plot(result)

### Response Curves:
library(gridExtra)

rc_m1 <- rcurve(m1, gg=T, mean=T, confidence=T, main = 'Response Curve')+ 
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        text=element_text(size=10))
rc_m2 <- rcurve(m2, gg=T, mean=T, confidence=T, main = 'Response Curve')+ 
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        text=element_text(size=10))
rc_m3 <- rcurve(m3, gg=T, mean=T, confidence=T, main = 'Response Curve')+ 
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        text=element_text(size=10))
rc_m4 <- rcurve(m4, gg=T, mean=T, confidence=T, main = 'Response Curve')+ 
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        text=element_text(size=10))
rc <- rbind(rc_m1$data,rc_m2$data,rc_m3$data,rc_m4$data)

rc_b15 <- rc[rc$variable=="b15",]
rc_NDVI <- rc[rc$variable=="NDVI",]

varm1 <- getVarImp(m1)
varm2 <- getVarImp(m2)
varm3 <- getVarImp(m3)
varm4 <- getVarImp(m4)

varimp_b15 <- (varm1@varImportanceMean$corTest$corTest[1] + 
               varm2@varImportanceMean$corTest$corTest[1] + 
               varm3@varImportanceMean$corTest$corTest[1] + 
               varm4@varImportanceMean$corTest$corTest[1])/4
varimp_NDVI <- (varm1@varImportanceMean$corTest$corTest[2] + 
                varm2@varImportanceMean$corTest$corTest[2] + 
                varm3@varImportanceMean$corTest$corTest[2] + 
                varm4@varImportanceMean$corTest$corTest[2])/4

p1 <- rc_NDVI %>% ggplot(aes(Value, Response)) +
  ylim(0,1) +
  geom_text(x = ((max(rc_NDVI$Value)-min(rc_NDVI$Value))*0.2 + min(rc_NDVI$Value)), y = 0.9, 
            label = paste0("VarImp. = ",
                           round(varimp_NDVI, 3)), 
            parse = F) +
  labs(x = "NDVI", y = "Probability") +
  geom_smooth(color="Black", span = 0.50, method = "loess", method.args = list(degree=1))

p2 <- rc_b15 %>% ggplot(aes(Value, Response)) +
  ylim(0,1) +
  geom_text(x = ((max(rc_b15$Value)-min(rc_b15$Value))*0.2 + min(rc_b15$Value)), y = 0.9, 
            label = paste0("VarImp. = ",
                           round(varimp_b15, 3)), 
            parse = F) +
  theme(axis.text.y=element_blank())+
  labs(x = "b15", y = "Probability") +
  geom_smooth(color="Black", span = 0.50, method = "loess", method.args = list(degree=1))

grid.arrange(p1, p2,
             nrow = 1,
             top = "Response Curves from Ensemble Models")

# End