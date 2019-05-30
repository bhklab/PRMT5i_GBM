

library(IncucyteDRC)
library(PharmacoGx)



## IMPORTANT - you need to set the directory to where you downloaded the github code [setwd("PRMT5i_GMB-master")]. 
setwd(".") 
source("incuCyte_helper_functions.R")


##########################################
##########################################
# Processing GXXX lines


test_pm_df <- importPlatemap(input = "Data/plateMap_G_data.txt") # The plate design of the incuCyte experiments
directory <- "Data/G_sensitivity_data/"

files.ls <- dir(directory)

SensitivityData <- lapply(files.ls, function(x){
  
  outputFile <- paste(directory,x,sep = "")
  
  # fitting the growth curves from incuCyte
  test_data <- importIncucyteData(outputFile, metric='pc')
  test_idrcset2 <- splitIncucyteDRCPlateData(test_pm_df, test_data, group_columns='celltype')
  test_idrcset2 <- fitGrowthCurvesGrouped(test_idrcset2)
  test_idrcset2 <- fitGrowthCurvesIndividual(test_idrcset2)

  
  
  th <- 0.95 # 95th percentile as a threshold of the total time to reach 100% confluency of the DMSO controls to allow for maximum drug effect and to avoid saturation which will be used to extract concentration and viabilities at that point.
  
  ibx <- which(test_idrcset2$fitted_data_grouped$sampleid=="DMSO")
  cuttime <- as.numeric(th)*length(ibx)
  test_idrcset2 <- calculateDRCData(test_idrcset2, cut_time = cuttime)

  test_idrcset2 <- fitDoseResponseCurve(test_idrcset2, include_control = F)
  controlGrowth3 <- median(test_idrcset2$drc_data[which(test_idrcset2$drc_data$sampleid=="DMSO"),"cut_val"])
  concViability3 <- get_Conc_Viability(test_idrcset2$drc_models$drc_model[[1]]$origData)
  
  # extracting IC50 and AAC (sensitivity measures) using PharmacoGx package.
  
  # GSK
  IC50_2.1 <- computeIC50(concViability3$conc,concViability3$viab/controlGrowth3, conc_as_log=FALSE, viability_as_pct=F)
  AAC_2.1 <- computeAUC(concViability3$conc,concViability3$viab/controlGrowth3, conc_as_log=FALSE, viability_as_pct=F)
  
  
  # LLY
  concViability4 <- get_Conc_Viability(test_idrcset2$drc_models$drc_model[[2]]$origData)
  
  IC50_2.2 <- computeIC50(concViability4$conc,concViability4$viab/controlGrowth3, conc_as_log=FALSE, viability_as_pct=F)
  AAC_2.2 <- computeAUC(concViability4$conc,concViability4$viab/controlGrowth3, conc_as_log=FALSE, viability_as_pct=F)
  
 
  
  return(c(IC50_2.1,AAC_2.1,IC50_2.2,AAC_2.2,cuttime))
  
})



metrics_for_sensitivity <- do.call(rbind,SensitivityData)

colnames(metrics_for_sensitivity) <- c("IC50_GSK591_0.95_cuttime","AAC_GSK591_0.95_cuttime","IC50_LLY283_0.95_cuttime","AAC_LLY283_0.95_cuttime","cuttime_95")


sampleIDs <- unlist(lapply(files.ls,function(x){
  strsplit(x,split = ".",fixed = T)[[1]][1]
}))

rownames(metrics_for_sensitivity) <- sampleIDs

# change 514 to 361
rownames(metrics_for_sensitivity)[3] <- "361"

rownames(metrics_for_sensitivity) <- paste0("G",rownames(metrics_for_sensitivity))




###########################################
###########################################
# Processing BTXXX lines

source("ComputeSensitivity.R")

viability1 <-"Data/BT_sensitivity_data.txt"

concentration <-  "Data/BT_conc.txt"

# compute sensitivity from almar blue assay using PharmacoGx
sensitivity <- computeSensitivity(viability1, concentration, fold.change=3, dilution.no=9, experiment.dim=c(11, 6))

finalProfile <- cbind(do.call(rbind,strsplit(rownames(sensitivity$profile),"_")),sensitivity$profile[,c("IC50","AUC")])
rownames(finalProfile) <- NULL
colnames(finalProfile) <- c("line","drug","IC50","AAC")
finalProfile$AAC <- finalProfile$AAC/100

sensitivityBT <- matrix(nrow = length(unique(finalProfile$line)),ncol = 4)
pos=1
for (i in seq(1,33,2)) {
  sensitivityBT[pos,] <- c(unlist(finalProfile[i+1,c(3,4),drop=T]),unlist(finalProfile[i,c(3,4),drop=T]))
  pos = pos + 1
}

rownames(sensitivityBT) <- unique(finalProfile$line)

#############################
#Integrating both data

metrics_for_sensitivity_final <- rbind(metrics_for_sensitivity[,1:4],sensitivityBT)

write.table(metrics_for_sensitivity_final,file = "Sensitivity_GBM_PRMT5i.csv",quote = F,sep = ",",col.names = NA,row.names = T)



