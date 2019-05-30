require(PharmacoGx)
require(magicaxis)
options(stringsAsFactors=FALSE)
#' @examples 
#' GDSCauc <- computeSensitivity(viability, concentration, fold.change=2, dilution.no=18)
#'
#' @param viability A text file that is exported from the plate reader 
#' Some information about these files:
#' the values are actually absorbance values.
#' the assay is a colorimetric assay so it doesn’t provide cell counts. 
#' it is needed to calculate a ““background” absorbance based on the outermost wells that have PBS and 
#' no cells that we use to background subtract the readings for wells in the 6 x 10 grid. 
#'  files with two blocks of data (each blocks corresponds to absorbance reads at a different wavelength)
#'  generally we should just have one block of data for one wavelength

#' @param concentration A text file manually created containing the maximum concentration of each drug
#' @param fold.change A value indicates the dilution fold 
#' @param dilution.no A value indicates the number of steps to create the serial dilution for treatment
#' @param replicates How many technical replicates placed on one plate
#' @param experiment.dim The design of the plate in terms of (# of rows, # of columns) including all background, control and treated wells
#' @param experiment.index.start A value indicates that in viability text file, experiments located how many lines after plate identifier
#' @param cutoff A threshold for detecting noisy wells
#' @return [matrix] A list of two matrices, sensitivity.info: annotation of experiments, 
#' sensitivity profiles: experiments summarized in metrics like IC50, AAC,
#' sensitivity.med: the concentartions and their corresponding median vaiability values for each experiment
#' raw.complete: the concentartions and their corresponding vaiability values (all technical replicates are kept) for each experiment
#' @export
computeSensitivity <- function(viability1, concentration, fold.change=3, dilution.no=9, replicates, experiment.dim=c(11, 6), experiment.index.start=3, cutoff=2) {
  viability <- read.csv(viability1, stringsAsFactors=FALSE, sep="\t", header=FALSE, na.strings=c("", " ", "NA"))
  concentration <- read.csv(concentration, stringsAsFactors=FALSE, sep="\t", na.strings=c("", " ", "NA"))
  if(missing(replicates)){
    if(experiment.dim[2] == 6){
      replicates <- experiment.dim[2] / 2
    }else{
      replicates <- experiment.dim[2]
    }
  }
  print(replicates)
  ##this function assumes in the text file produced by the machine always there is field called "Plate:" before plate names  
  indicies <- grep("Plate", viability$V1)
  plates <- viability$V2[indicies]
  
  ##plate names should match between viability and concentration files
#   if(!all(plates == concentration$plate)){
#     message("plates are dicordant between viability and concentration files!")
#     exit()
#   }
  data.blocks <- as.numeric(table(is.na(viability[indicies[1] + experiment.index.start, ]))["FALSE"] / experiment.dim[1])
  wells.with.cells.no <- as.numeric(table(is.na(viability[indicies[1] + experiment.index.start, ]))["FALSE"])
  
  xx <- is.na(viability[indicies[1] + experiment.index.start, ])
  yy <- NULL
  yy[1] <- 0
  for(i in 2:length(xx)){
    if((xx[i - 1] == TRUE) && (xx[i] == FALSE)){
      yy[i] <- 1
     }else{
        yy[i] <- 0
    }
   }
  starts <- which(yy==1)
#     raw.sensitivity.rows <- do.call(c, lapply(1:data.blocks, function(x){paste(plates, x, sep="_")}))
#   }else {
#     raw.sensitivity.rows <- plates
#   }
  
  raw.sensitivity.rows <- NULL
  for(plate in plates) {
    
 #   xx <- unlist(strsplit(gsub(sprintf("_[1-%s]$", dilution.no), "", plate), split=' (?=[^ ]+$)', perl=TRUE))
#    cell <- paste(xx[1:(length(xx)-1)], sep="-")
#    drugs <- unlist(strsplit(xx[length(xx)], split='/'))
    
    cell <- strsplit(plate,split = "-")[[1]][1]
    drugs <- unlist(strsplit(strsplit(plate,split = "-")[[1]][2],split = "/"))
    
    if(data.blocks > 1){
      raw.sensitivity.rows <- c(raw.sensitivity.rows, do.call(c, lapply(1:data.blocks, function(x){paste(paste(cell, drugs, sep="_"), x, sep="_")})))
    }else{
      raw.sensitivity.rows <- c(raw.sensitivity.rows, paste(cell, drugs, sep="_"))
    }
  }
    
   ##Assumption: 6 wells are assigned to control and 6 to controls and treated wells are 18 in triplicares (altogether 66 wells)
  #data.blocks * 2 = PBS + CONTROL(UNTREATED)
  raw.sensitivity <- array(NA, dim=c(length(raw.sensitivity.rows), dilution.no + data.blocks * 2, replicates + 1), dimnames=list(raw.sensitivity.rows, NULL, c("Dose", paste0("Viability", 1:replicates))))
  ##this function assumes all the actual viability values appear from the experiment.index.start lines belower than the Plate name
  ##with 18 concentrations and two controls in triplicates 
  #indicies <- indicies + experiment.index.start
  for(index in indicies) {
    plate <- viability$V2[index]
    
  #  xx <- unlist(strsplit(gsub(sprintf("_[1-%s]$", dilution.no), "", plate), split=' (?=[^ ]+$)', perl=TRUE))
  #  cell <- paste(xx[1:(length(xx)-1)], sep=" ")
  #  drugs <- unlist(strsplit(xx[length(xx)], split='/'))
    
    cell <- strsplit(plate,split = "-")[[1]][1]
    drugs <- unlist(strsplit(strsplit(plate,split = "-")[[1]][2],split = "/"))
    
    for(j in 1:data.blocks){
      start <- starts[j]
      for(dd in 1:length(drugs)) {
        drug <- drugs[dd]
        mm <- as.matrix(viability[(index + experiment.index.start + (dd - 1) * replicates * data.blocks):
                          (index + experiment.index.start + dd * replicates * data.blocks - 1),
                        start:
                          (start + experiment.dim[1] - 1)])
        if(data.blocks > 1){
          exp <- paste(paste(cell, drug, sep="_"), j, sep="_")
        }else{
          exp <- paste(cell, drug, sep="_")
        }
        doses <- concentration$max.concentration.uM.[which(concentration$drug == drug)]
        tt <- as.numeric(doses)
        for(k in 1:(dilution.no - 1)){
          tt <- tt/fold.change
          doses <- c(doses, tt)
        }
        raw.sensitivity[exp, , "Dose"] <- c(rep("pbs", times=data.blocks), rep("control", times=data.blocks), rev(doses) )
        pbs <- grep("pbs", raw.sensitivity[exp, , "Dose"])
        control <- grep("control", raw.sensitivity[exp, , "Dose"])
        treated <- match(rev(doses), raw.sensitivity[exp, , "Dose"])
        
        raw.sensitivity[exp, pbs, paste0("Viability", 1:replicates)] <- matrix(mm[,pbs], ncol=replicates, nrow=data.blocks)
        raw.sensitivity[exp, control, paste0("Viability", 1:replicates)] <- matrix(as.numeric(mm[(nrow(mm) - (replicates - 1)): nrow(mm) ,control]), ncol=replicates, nrow=data.blocks)
        tt <- t(mm[1:replicates, treated])
        if(nrow(mm) > replicates){
          tt <- rbind(t(mm[1:replicates, 2:ncol(mm)]), t(mm[(replicates + 1):nrow(mm), 2:(ncol(mm) - (data.blocks))]))
        }
        raw.sensitivity[exp, treated, paste0("Viability", 1:replicates)] <- tt
      }
    }
  }
  
  viability.col <- grep("Viability", unlist(dimnames(raw.sensitivity)[3]))
  pbs.row <- grep("pbs", raw.sensitivity[1, , "Dose"])
  control.row <- grep("control", raw.sensitivity[1, , "Dose"])
  treated.row <- setdiff(1:nrow(raw.sensitivity[1, , ]), union(pbs.row, control.row))
  pbs <- mean(as.numeric(raw.sensitivity[, pbs.row, viability.col]), na.rm=TRUE)
  
  raw.sensitivity.clean <- array(NA, dim=c(length(raw.sensitivity.rows), length(treated.row), dim(raw.sensitivity)[3]), dimnames=list(raw.sensitivity.rows, NULL, unlist(dimnames(raw.sensitivity)[3])))
  raw.sensitivity.med <- array(NA, dim=c(length(raw.sensitivity.rows), length(treated.row), 2), dimnames=list(raw.sensitivity.rows, NULL, c("Dose", "Viability")))
  sensitivity.info <- data.frame(matrix(NA, ncol=2, nrow=length(raw.sensitivity.rows), dimnames=list(raw.sensitivity.rows, c("cell", "drug"))))
  sensitivity.profile <- data.frame(matrix(NA, ncol=5, nrow=length(raw.sensitivity.rows), dimnames=list(raw.sensitivity.rows, c("IC50", "AUC", "Hill.Slope", "Einf", "EC50"))))
  for(exp in rownames(raw.sensitivity)) {

    ##remove background noise
    tt <- apply(raw.sensitivity[exp, , viability.col], MARGIN = 2, as.numeric)
    tt <- tt - median(tt[pbs.row,])
    
    ##flag and remove outliers
    tt <- cbind(tt, "stdev"= apply(tt, MARGIN=1, sd, na.rm=TRUE))
    tt <- cbind(tt, "coef.var"= tt[, "stdev"]/ apply(tt, MARGIN=1, mean, na.rm=TRUE))


    tt[control.row, "stdev"] <- sd(tt[control.row, 1:length(viability.col)], na.rm=TRUE)
    tt[control.row, "coef.var"] <- tt[control.row, "stdev"]/ mean(tt[control.row, 1:3], na.rm=TRUE)
    
    noisy.points <- which(abs(tt[c(control.row, treated.row), "stdev"]) > cutoff | abs(tt[c(control.row, treated.row), "coef.var"]) > cutoff)
    for(p in noisy.points) {
      tt[p, which.max((tt[p, 1:length(viability.col)] - mean(tt[p, 1:length(viability.col)], na.rm=TRUE)) ^ 2)] <- NA
    }
    
    ##normalize by control
    control <- median(tt[control.row, 1:length(viability.col)], na.rm=TRUE)
    xx <- tt[treated.row, 1:length(viability.col)] / control
    
    raw.sensitivity.clean[exp, , "Dose"] <- raw.sensitivity[exp, treated.row, "Dose"] 
    raw.sensitivity.clean[exp, , viability.col] <- xx
    
    raw.sensitivity.med[exp, , "Dose"] <- as.numeric(rev(raw.sensitivity.clean[exp, , "Dose"]))
    raw.sensitivity.med[exp, , "Viability"] <- as.numeric(rev(apply(xx, MARGIN=1, median, na.rm=TRUE) * 100))
    
    xx <- unlist(strsplit(gsub(sprintf("_[1-%s]$", dilution.no), "", exp), split='_(?=[^_]+$)', perl=TRUE))
    sensitivity.info[exp, "cell"] <- xx[1]
    sensitivity.info[exp, "drug"] <- xx[2]
    
    sensitivity.profile[exp, "IC50"] <- PharmacoGx::computeIC50(concentration=raw.sensitivity.med[exp, , "Dose"], viability=raw.sensitivity.med[exp, , "Viability"])
    sensitivity.profile[exp, "AUC"] <- PharmacoGx::computeAUC(concentration=raw.sensitivity.med[exp, , "Dose"], viability=raw.sensitivity.med[exp, , "Viability"])
    sensitivity.profile[exp, c("Hill.Slope", "Einf", "EC50")] <- unlist(PharmacoGx::logLogisticRegression(conc=raw.sensitivity.med[exp, , "Dose"], viability=raw.sensitivity.med[exp, , "Viability"]))
    
  }
  return(list("info"=sensitivity.info, "profile"=sensitivity.profile, "raw"=raw.sensitivity.med, "raw.complete"=raw.sensitivity.clean))
}

drugDoseResponseCurve <- 
  function(drug, cellline, conc, viability, ylim, xlim, mycol, plot.type=c("Fitted","Actual", "Both"), legend) {
    
    doses <- list(); responses <- list(); legend.values <- list(); j <- 0; pSetIndex <- list()
    doses[[1]] <- conc
    responses[[1]] <- viability

    if (missing(mycol)) {
      require(RColorBrewer) || stop("Library RColorBrewer is not available!")
      mycol <- RColorBrewer::brewer.pal(n=7, name="Set1")
    }
    
    dose.range <- c(10^100 , 0)
    viability.range <- c(0 , 10)
    for(i in 1:length(doses)) {
      dose.range <- c(min(dose.range[1], min(doses[[i]], na.rm=TRUE), na.rm=TRUE), max(dose.range[2], max(doses[[i]], na.rm=TRUE), na.rm=TRUE))
      viability.range <- c(0, max(viability.range[2], max(responses[[i]], na.rm=TRUE), na.rm=TRUE))
    }
    x1 <- 10 ^ 10; x2 <- 0
    
    if(length(doses) > 1) {
      common.ranges <- .getCommonConcentrationRange(doses)
      
      for(i in 1:length(doses)) {
        x1 <- min(x1, min(common.ranges[[i]]))
        x2 <- max(x2, max(common.ranges[[i]]))
      }
    }
    if (!missing(xlim)) {
      dose.range <- xlim
    }
    if (!missing(ylim)) {
      viability.range <- ylim
    }
    
    plot(NA, xlab="Concentration (uM)", ylab="% Viability", axes =FALSE, main=sprintf("%s:%s", drug, cellline), log="x", ylim=viability.range, xlim=dose.range, cex=.7, cex.main=.9)
    magicaxis::magaxis(side=1:2, frame.plot=TRUE, tcl=-.3, majorn=c(5,3), minorn=c(5,2))
    legends <- NULL
    legends.col <- NULL
    if (length(doses) > 1) {
      rect(xleft=x1, xright=x2, ybottom=viability.range[1] , ytop=viability.range[2] , col=rgb(240, 240, 240, maxColorValue = 255), border=FALSE)
    }
    
    for (i in 1:length(doses)) {
      points(doses[[i]],responses[[i]],pch=20,col = mycol[i])
      
      switch(plot.type , "Actual"={
        lines(doses[[i]], responses[[i]], lty=1, lwd=.5, col=mycol[i])
      }, "Fitted"={ 
        log_logistic_params <- logLogisticRegression(conc=doses[[i]], viability=responses[[i]])
        log10_x_vals <- PharmacoGx:::.GetSupportVec(log10(doses[[i]]))
        lines(10 ^ log10_x_vals, PharmacoGx:::.Hill(log10_x_vals, pars=c(log_logistic_params$HS, log_logistic_params$E_inf/100, log10(log_logistic_params$EC50))) * 100 ,lty=1, lwd=.5, col=mycol[i])
      },"Both"={
        lines(doses[[i]],responses[[i]],lty=1,lwd=.5,col = mycol[i])
        log_logistic_params <- PharmacoGx::logLogisticRegression(conc = doses[[i]], viability = responses[[i]])
        log10_x_vals <- PharmacoGx:::.GetSupportVec(log10(doses[[i]]))
        lines(10 ^ log10_x_vals, PharmacoGx:::.Hill(log10_x_vals, pars=c(log_logistic_params$HS, log_logistic_params$E_inf/100, log10(log_logistic_params$EC50))) * 100 ,lty=1, lwd=.5, col=mycol[i])
      })
    }
    
    legend("topright", legend=legend, col=mycol, bty="n", cex=.7, pch=c(15,15))
  }

filterNoisyCurves <- function(sensitivity, epsilon=25 , positive.cutoff.percent=.80, mean.viablity=200, nthread=1) {
  
  acceptable <- mclapply(rownames(sensitivity$info), function(xp) {
    #for(xp in rownames(sensitivityInfo(pSet))){
    drug.responses <- as.data.frame(apply(sensitivity$raw[xp , ,], 2, as.numeric), stringsAsFactors=FALSE)
    drug.responses <- drug.responses[complete.cases(drug.responses), ]
    doses.no <- nrow(drug.responses)
    
    drug.responses[,"delta"] <- .computeDelta(drug.responses$Viability)
    
    delta.sum <- sum(drug.responses$delta, na.rm = TRUE)
    
    max.cum.sum <- .computeCumSumDelta(drug.responses$Viability)
    
    if ((table(drug.responses$delta < epsilon)["TRUE"] >= (doses.no * positive.cutoff.percent)) &
        (delta.sum < epsilon) &
        (max.cum.sum < (2 * epsilon)) &
        (mean(drug.responses$Viability) < mean.viablity)) {
      return (xp)
    }
  }, mc.cores=nthread)
  acceptable <- unlist(acceptable)
  noisy <- setdiff(rownames(sensitivity$info), acceptable)
  return(list("noisy"=noisy, "ok"=acceptable))
}

.computeDelta <- function(xx ,trunc = TRUE) {
  xx <- as.numeric(xx)
  if(trunc)
  {
    return(c(pmin(100, xx[2:length(xx)]) - pmin(100, xx[1:length(xx)-1]), 0))
  }else{
    return(c(xx[2:length(xx)] - xx[1:length(xx)-1]), 0)
  }
}

.computeCumSumDelta <- function(xx, trunc = TRUE) {
  xx <- as.numeric(xx)
  if(trunc) {
    xx <- pmin(xx, 100)
  }
  tt <- t(combn(1:length(xx), 2 , simplify = T))
  tt <- tt[which(((tt[,2] - tt[,1]) >= 2) == T),]
  cum.sum <- unlist(lapply(1:nrow(tt), function(x){xx[tt[x,2]]-xx[tt[x,1]]}))
  return(max(cum.sum))
}

handleNoisyCurves <- function(sensitivity, xp, epsilon=25 , positive.cutoff.percent=.80, mean.viablity=200) {
  
    handler <- NULL
    auc <- NULL
    original.auc <- PharmacoGx::computeAUC(concentration=sensitivity$raw[xp , , "Dose"], viability=sensitivity$raw[xp , , "Viability"])
    for(i in 1:ncol(sensitivity$raw)) {
      tt <- setdiff(1:ncol(sensitivity$raw), i)
      drug.responses <- as.data.frame(apply(sensitivity$raw[xp , tt,], 2, as.numeric), stringsAsFactors=FALSE)
      drug.responses <- drug.responses[complete.cases(drug.responses), ]
      doses.no <- nrow(drug.responses)
      
      drug.responses[,"delta"] <- .computeDelta(drug.responses$Viability)
      
      delta.sum <- sum(drug.responses$delta, na.rm = TRUE)
      
      max.cum.sum <- .computeCumSumDelta(drug.responses$Viability)
      
      if ((table(drug.responses$delta < epsilon)["TRUE"] >= (doses.no * positive.cutoff.percent)) &
          (delta.sum < epsilon) &
          (max.cum.sum < (2 * epsilon)) &
          (mean(drug.responses$Viability) < mean.viablity)) {
        handler <- c(handler, i)
        auc <- c(auc, PharmacoGx::computeAUC(concentration=sensitivity$raw[xp , tt, "Dose"], viability=sensitivity$raw[xp , tt, "Viability"]))
      }
    }
    for(i in 1:ncol(sensitivity$raw)) {
      for(j in setdiff(i:ncol(sensitivity$raw), i))
      {
        tt <- setdiff(1:ncol(sensitivity$raw), c(i,j))
        drug.responses <- as.data.frame(apply(sensitivity$raw[xp , tt,], 2, as.numeric), stringsAsFactors=FALSE)
        drug.responses <- drug.responses[complete.cases(drug.responses), ]
        doses.no <- nrow(drug.responses)
        
        drug.responses[,"delta"] <- .computeDelta(drug.responses$Viability)
        
        delta.sum <- sum(drug.responses$delta, na.rm = TRUE)
        
        max.cum.sum <- .computeCumSumDelta(drug.responses$Viability)
        
        if ((table(drug.responses$delta < epsilon)["TRUE"] >= (doses.no * positive.cutoff.percent)) &
            (delta.sum < epsilon) &
            (max.cum.sum < (2 * epsilon)) &
            (mean(drug.responses$Viability) < mean.viablity)) {
          handler <- c(handler, sprintf("%s_%s", i, j))
          auc <- c(auc, PharmacoGx::computeAUC(concentration=sensitivity$raw[xp , tt, "Dose"], viability=sensitivity$raw[xp , tt, "Viability"]))
        }
      }
    }
    
  return(list("noisy"=noisy, "ok"=acceptable))
}
