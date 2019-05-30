


get_Conc_Viability <- function(concViability){
  
  
  noOfReps <- length(table(concViability$idx))
  
  conc <- unique(concViability$conc)
  conc <- sort(conc)
  
  output <- unlist(lapply(conc, function(c){
    viab <- concViability$value[concViability$conc==c]
    return(median(viab))
  }))
  
  
  return(list("conc"=conc,"viab"=output))
}



plotIncucyteDRCSet_Modified <- function (idrc_set, grouped = FALSE) 
{
  library(plyr)
  if (grouped & is.null(idrc_set$fitted_data_grouped)) {
    stop("Need to fit growth curves first using fitGrowthCurvesGrouped")
  }
  if (!grouped & is.null(idrc_set$fitted_data_grouped)) {
    stop("Need to fit growth curves first using fitGrowthCurvesIndividual")
  }
  library(dplyr)
  data1 <- idrc_set$platemap %>% dplyr::inner_join(idrc_set$platedata$data, by = "wellid")
  out_plot <- ggplot(data1, aes(x = elapsed, y = value, colour = as.factor(round(conc, 3)))) + geom_point(alpha = 0.2) + scale_colour_discrete(name = "conc") + facet_wrap(~sampleid) + theme_bw()
  if (grouped) {
    out_plot <- out_plot + geom_line(data = idrc_set$fitted_data_grouped, 
                                     aes(group = paste0(sampleid, conc)))
  }  else {
    out_plot <- out_plot + geom_line(data = idrc_set$fitted_data_indiv, 
                                     aes(group = wellid))
  }
  if (is.numeric(data1$cuttime.AAC)) {
    out_plot <- out_plot + geom_vline(data=ddply(data1, "sampleid", summarize, cut_time = unique(cuttime.AAC)), aes(xintercept=cut_time), colour = "blue", linetype = "dashed")
  }
  return(out_plot)
 
}

