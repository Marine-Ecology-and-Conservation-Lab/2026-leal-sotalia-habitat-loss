if(!require(rgdal)) install.packages("rgdal", dependencies = T)
if(!require(sf)) install.packages("sf", dependencies = T)
if(!require(raster)) install.packages("raster", dependencies = T)
if(!require(spThin)) install.packages("spThin", dependencies = T)
if(!require(usdm)) install.packages("usdm", dependencies = T)
if(!require(biomod2)) install.packages("biomod2", dependencies = T)
if(!require(dplyr)) install.packages("dplyr", dependencies = T)
if(!require(ggplot2)) install.packages("ggplot2", dependencies = T)

#### SCRIPT PATH ####
getwd()
setwd("C:/pres")

#### NAMING THE SPECIES ####
Species <- "Sotalia_guianensis"

#### Loading RData ####
# load(paste(Species,".RData", sep = ""))

#### Getting data ####


#Occurrence data
sotalia=read_sf("./pontos_pres.shp")
View(sotalia)
#### Environmental variables ####

chl=raster("./cloro.tif")
dep=raster("./prof.tif")
dist=raster("./coast_R_final2.tif")
slo=raster("./Slope.tif")
#sst=raster("./sst_sir_sep1R.tif")
#fun=raster("./fundeio.tif")
#canal=raster("./canal.tif")
#porto=raster("./porto.tif")
antr=raster("./antr.tif")

clo = (resample(chl,dep,method="ngb"))
slop = (resample(slo,dep,method="ngb"))
dista= (resample(dist,dep,method="ngb"))
antro= (resample(antr,dep,method="ngb"))

bbox(clo); ncol(clo); nrow(clo) ; res(clo)
bbox(slop); ncol(slop); nrow(slop) ; res(slop)
bbox(dep); ncol(dep); nrow(dep) ; res(dep)
bbox(dista); ncol(dista); nrow(dista) ; res(dista)
bbox(antro); ncol(antro); nrow(antro) ; res(antro)

sotalia_transformado <- st_transform(sotalia, st_crs(dep))

# Merging layers
biostack <- stack(clo, dep, slop, dista, antro)
# biostack <- stack(clo, dep, slop, dista, antro)
pal <- colorRampPalette(c("yellow", "blue"))

names(biostack)[1] <- "Chlorophyll"
names(biostack)[2] <- "Depth"
names(biostack)[3] <- "Slope"
names(biostack)[4] <- "Distance"
#names(biostack)[5] <- "Temperature"
names(biostack)[5] <- "Cumulative"

plot(biostack, col = pal(20))

# Correlation tests -------------------------------------------------------


pairs(biostack) # Correlation
biostack_df <- as.data.frame(biostack)
vif(biostack_df) # Multicollinearity
warnings()
#vifstep(biostack_df)
vifcor(biostack_df, th = 0.7)

# Merging layers without collinearity problem
#biostack <- stack(chl, dep, sal, slo)
png("1. Biostack.png", width = 3200, height = 2400, res = 300); par(oma = c(1,1,1,1)); plot(biostack, col = pal(20)); dev.off()


# Presence data -----------------------------------------------------------
View(sotalia_transformado)
sot<-cbind(sotalia_transformado$Lati, sotalia_transformado$Long) 
colnames(sot)<-c("Latitude","Longitude")
sot<-as.data.frame(sot) 
sot$Especie="Sotalia guianensis"

sot_na <- sot %>% tidyr::drop_na(Longitude, Latitude)
nrow(sot)-nrow(sot_na)

#Checking inconsistencies
sot_na$Longitude<-as.numeric(sot_na$Longitude)
sot_na$Latitude<-as.numeric(sot_na$Latitude)
flags_spatial<-CoordinateCleaner::clean_coordinates(x = sot_na, species = "Especie",lon = "Longitude", lat = "Latitude", tests = c("capitals", # raio ao redor de capitais
                                                                                                                                   "centroids", # raio ao redor de centroides de paises e provincias
                                                                                                                                   "duplicates", # duplicatas
                                                                                                                                   "equal", # coordenadas iguais
                                                                                                                                   "gbif", # raio ao redor da sede da GBIF
                                                                                                                                   "institutions", # raio ao redor de instituicoes de pesquisa em biodiversidade
                                                                                                                                   "urban", # pontos dentro de areas urbanas
                                                                                                                                   "validity", # ponto de fora do sistema de coordenadas
                                                                                                                                   "zeros" # zeros e pontos onde lat = lon 
)
)

#FILTERATING INCONSISTENCIES#
sot_f <- sot_na %>% 
  dplyr::filter(flags_spatial$.summary == TRUE)

#SELECTING COLUMNS#
sot_f = sot_f %>%
  dplyr::select(Especie, Longitude, Latitude)

PresenceData <- sot_f

head(PresenceData)

# Plotting presence data
plot(chl, xlim=c(-44.9,-43.4),ylim=c(-23.4,-22.9), main="Pontos de ocorrencia por Temperatura")
points(PresenceData$Longitude, PresenceData$Latitude, col="red", cex=1)

write.csv2(PresenceData, "0. pontosfiltrados.csv")

# Thinning presence data --------------------------------------------------

dir.create("Thinned")

PresenceDataThinned <- thin( loc.data = PresenceData,
                             lat.col = "Latitude",
                             long.col = "Longitude",
                             spec.col = "Especie",
                             thin.par = 1,
                             reps = 5,
                             locs.thinned.list.return = TRUE,
                             write.files = TRUE,
                             max.files = 5,
                             out.dir = "Thinned",
                             out.base = "Presence_thinned",
                             write.log.file = TRUE,
                             log.file = "Thinned/Presence_thinned_log_file.txt" )

#Loading thinned presences file#
sot_thin_1km=read.csv("Thinned/Presence_thinned_thin1.csv", sep=",")
plot(sot_thin_1km)
getwd()
View(sot_thin_1km)

sot_thin_1km["Occ"] <- c("1")

#### Prepare data & parameters ####

# Select the name of the studied species

myRespName <- "Sotalia_guianensis"

# Get corresponding presence/absence data

myResp <- as.numeric(sot_thin_1km$Occ)
str(myResp)

myRespXY <- sot_thin_1km[which(myResp==1),c("Longitude", "Latitude")]
colnames(myRespXY)


# Pseudo-absences extraction ----------------------------------------------

bm.Species <- BIOMOD_FormatingData(resp.var = myResp,
                                   expl.var = biostack,
                                   resp.xy = myRespXY,
                                   resp.name = myRespName,
                                   PA.nb.rep = 3,
                                   PA.nb.absences = (nrow(sot_thin_1km)*3),
                                   PA.strategy = "disk", 
                                   PA.dist.min = 1000)

# Summary
bm.Species
sum(bm.Species@PA.table[["PA1"]] == TRUE)-nrow(sot_thin_1km) # n of PAs

# Plotting pseudo-absences
plot(bm.Species)
png("2. Pseudo-absences.png", width = 2400, height = 2400, res = 300); plot(bm.Species); dev.off()

#### Parameterize modeling options ####

# Preparing for Maxent

# Saving explanatory variables in .ascii for MaxEnt algorithm

dir.create("maxent_bg")
maxent.background.dat.dir <- "maxent_bg"

for(var_ in names(biostack)){
  cat("/n> saving", paste0(var_, ".asc"))
  writeRaster(subset(biostack, var_), 
              filename = file.path(maxent.background.dat.dir, paste0(var_, ".asc")),
              overwrite = TRUE)
}

# Defining the path for maxent.jar file 

path.to.maxent.jar <- file.path("C:/pres", "maxent.jar")

# Defining Models Options

bm.opt <- BIOMOD_ModelingOptions(GAM = list(k = 4),
                                 MAXENT = list(path_to_maxent.jar = path.to.maxent.jar,
                                               background_data_dir = maxent.background.dat.dir,
                                               maximumbackground = 10000))


#### Running single models ####

{sotaliamodel <- BIOMOD_Modeling(bm.Species,
                                 modeling.id = paste("model_", Species, sep=""),
                                 models = c("GLM", "GAM", "RF", "GBM", "MAXENT"),
                                 bm.options = bm.opt,
                                 CV.nb.rep = 10,
                                 CV.strategy = "kfold",
                                 CV.k = 5,
                                 CV.perc = 0.7,
                                 do.full.models = F,
                                 prevalence = 0.5,
                                 var.import = 3,
                                 metric.eval = c("ROC"))

beepr::beep(8)}
warnings()

# Get evaluation scores & variables importance

ModelEvalsot <- get_evaluations(sotaliamodel)
write.csv2(ModelEvalsot, "3. ModelEval_single.csv")

VarImportsot <- get_variables_importance(sotaliamodel)
write.csv(VarImportsot, "4. VarImport_single.csv")
imagem=bm_PlotVarImpBoxplot(bm.out = sotaliamodel, group.by = c('expl.var', 'algo', 'algo'))
png(" algor.png", width = 2400, height = 2400, res = 300); imagem; dev.off()

# Represent evaluation scores & variables importance

Graph_EvalSinglesot <- bm_PlotEvalMean(bm.out = sotaliamodel, dataset = "validation")
png("5. EvalSingle.png", width = 2400, height = 2400, res = 300); Graph_EvalSinglesot; dev.off()

Graph_BPEvalSinglesot <- bm_PlotEvalBoxplot(bm.out = sotaliamodel, dataset = "validation", group.by = c('algo', 'algo'))
png("6. BoxplotEvalSingle.png", width = 2400, height = 2400, res = 300); Graph_BPEvalSinglesot; dev.off()

ModelEval_ROC <- ModelEvalsot %>% filter(metric.eval == "ROC")
newalgoROC <- factor(ModelEval_ROC$algo, levels = c("GLM", "GAM", "RF", "GBM", "MAXENT"))

NewGraph_EvalSinglesotROC <-  ggplot(ModelEval_ROC, aes(x = newalgoROC, y = validation, fill = algo)) +
  stat_boxplot(geom = "errorbar") + 
  geom_hline(yintercept = 0.8, linetype='longdash', col='red') +
  labs(title = paste(Species), x = "Algorithms", y = "ROC")+
  guides(fill = "none", color = "none") +
  theme_classic() +
  geom_boxplot() +
  ylim(0.15, 1) +
  theme(legend.title = element_blank(),
        title = element_text(size = 16, color = "gray33"),
        axis.text = element_text(size = 12, color = "gray33"),
        axis.title = element_text(size = 12, color = "gray33"))

plot(NewGraph_EvalSinglesotROC)

png("7. NewGraph_EvalSingleROC.png", width = 1800, height = 1800, res = 300); NewGraph_EvalSinglesotROC; dev.off()

#### Running ensemble model ####

{sot_EM <- BIOMOD_EnsembleModeling(sotaliamodel,
                                   models.chosen = "all",
                                   em.by = "all",
                                   metric.select = c("ROC"),
                                   em.algo = c("EMmean", "EMwmean", "EMca", "EMcv"),
                                   metric.select.thresh = c(0.8),
                                   metric.eval = c("ROC"),
                                   var.import = 3)

beepr::beep(8)}

?BIOMOD_EnsembleModeling()
summary(sot_EM)
View(sot_EM)
get_predictions(sot_EM)


{sot_EM1 <- BIOMOD_EnsembleModeling(sotaliamodel,
                                    models.chosen = "all",
                                    em.by = "PA+run",
                                    metric.select = c("ROC"),
                                    em.algo = c("EMmean", "EMwmean", "EMca", "EMcv"),
                                    metric.select.thresh = c(0.8),
                                    metric.eval = c("ROC"),
                                    var.import = 3)
  
  beepr::beep(8)}

?BIOMOD_EnsembleModeling()
summary(sot_EM1)
get_predictions(sot_EM1)

#### Ensemble models Evaluation ####
evalsotensemble = get_evaluations(sot_EM)
evalsotensemble2 = get_evaluations(sot_EM1)
View(evalsotensemble)
write.csv2(evalsotensemble, "9. evalsotensemble.csv")
View(evalsotensemble2)
write.csv2(evalsotensemble2, "11. evalsotensemble2.csv")

#PLOTTING PERFORMANCE#
#MODELS GLM, GAM, GBM, RF, MAXENT
ens_sot_eval =  bm_PlotEvalMean(bm.out = sot_EM)

length(sot_EM)
summary(sot_EM@em.models_kept)

# Get evaluation scores & variables importance

ModelEval_EM <- get_evaluations(sot_EM)
write.csv2(ModelEval_EM,paste0("7. ModelEval_EM.csv"))

VarImport_EM <- get_variables_importance(sot_EM)
write.csv2(VarImport_EM, "8. VarImport_EM.csv")

# Represent evaluation scores & variables importance

Graph_BoxplotVarAlgo <- bm_PlotVarImpBoxplot(bm.out = sot_EM, group.by = c('expl.var', 'algo', 'algo'))
png("9. BoxplotVarAlgo.png", width = 2400, height = 2400, res = 300); Graph_BoxplotVarAlgo; dev.off()
?bm_PlotVarImpBoxplot()

Graph_BoxplotAlgoVar <- bm_PlotVarImpBoxplot(bm.out = sot_EM, group.by = c('algo', 'expl.var', 'algo'))
png("10. BoxplotAlgoVar.png", width = 2400, height = 2400, res = 300); Graph_BoxplotAlgoVar; dev.off()

NewGraph_BoxplotVarImp <-  Graph_BoxplotVarAlgo$tab %>%
  filter(algo == "EMwmean") %>%
  ggplot() +
  geom_boxplot(aes(x = expl.var, y = var.imp, fill = factor(expl.var), color = factor(expl.var))) +
  scale_fill_manual(values = c("white", "white", "white", "white", "white", "white")) +
  scale_color_manual(values = c("red", "blue", "green", "orange", "cyan3", "brown4")) +
  theme_classic() +
  labs(title = paste(Species), x = "Explanatory variables", y = "Variable importance") +
  guides(fill = "none", color = "none") +
  ylim(0, 1) +
  theme(title = element_text(size = 16, color = "gray33"),
        axis.text = element_text(size = 12, color = "gray33"),
        axis.title = element_text(size = 12, color = "gray33"))

plot(NewGraph_BoxplotVarImp)

png("11. NewBoxplot_wmean.png", width = 1800, height = 1800, res = 300); NewGraph_BoxplotVarImp; dev.off()

# Represent response curves

Graph_ResponseCurves_wmean <- bm_PlotResponseCurves(bm.out = sot_EM,
                                                    models.chosen = get_built_models(sot_EM)[2],
                                                    fixed.var = 'median')
png("12. ResponseCurves_wmean.png", width = 2400, height = 2400, res = 300); Graph_ResponseCurves_wmean; dev.off()

ggdat <- Graph_ResponseCurves_wmean$tab
ggdat$algo <- ifelse(grepl("_GLM", ggdat$pred.name, ignore.case = T), "GLM",
                     ifelse(grepl("_GAM", ggdat$pred.name, ignore.case = T ), "GAM",
                            ifelse(grepl("_GBM", ggdat$pred.name, ignore.case = T ), "GBM",
                                   ifelse(grepl("_RF", ggdat$pred.name, ignore.case = T ), "RF",
                                          ifelse(grepl("_MAXENT", ggdat$pred.name, ignore.case = T ), "MAXENT", "NA")))))

ggdat$algo <- as.factor(ggdat$algo)

ggdat$id <- as.factor(ggdat$id)

NewResponseCurve <- ggplot(ggdat, aes(x = expl.val, y = pred.val)) +
  geom_smooth(aes(group = algo, color = algo), fill = "indianred1", se = TRUE, method = "gam", alpha = 0.01) +
  geom_smooth(color = "red", fill = "indianred1", se = TRUE, method = "gam") +
  facet_wrap(~ expl.name, scales = "free_x") +
  labs(title = paste(Species), x = "Explanatory variables value", y = "Suitability value") +
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.position = "",
        title = element_text(size = 16, color = "gray33"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12, color = "gray33"),
        strip.text = element_text(size = 12),
        panel.spacing.x = unit(2, "lines"),
        panel.spacing.y = unit(1, "lines"))

plot(NewResponseCurve)

png("13. NewResponseCurves_wmean.png", width = 2400, height = 1800, res = 300); NewResponseCurve; dev.off()


#### Projecting single models BIOMOD  ####

{Proj_sot <- BIOMOD_Projection(bm.mod = sotaliamodel,
                               proj.name = "Current",
                               new.env = biostack,
                               models.chosen = "all",
                               metric.binary = "all",
                               metric.filter = "all",
                               build.clamping.mask = F,
                               keep.in.memory = F)

beepr::beep(8)}


#### Projecting ensemble model ####

{Proj_sot_EM <- BIOMOD_EnsembleForecasting(bm.em = sot_EM,
                                           bm.proj = Proj_sot,
                                           proj.name = "Current_EM",
                                           models.chosen = "all",
                                           metric.binary = "all",
                                           metric.filter = "all",
                                           output.format = ".img")

beepr::beep(8)}


save.image(file = (paste("../../", Species, ".RData", sep = "")))

