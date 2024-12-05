library(SingleCellExperiment)
library(glmGamPoi)
library(dplyr)
sce = as.SingleCellExperiment(seuratobject)
##Aggregating counts
aggr_counts <- glmGamPoi::pseudobulk(sce, group_by = vars(orig.ident), aggregation_functions = list(counts = "rowSums2", .default = "rowMeans2"))
pseudo_counts = data.frame(aggr_counts@assays@data$counts)
pseudo_counts = ADImpute::NormalizeTPM(pseudo_counts, log = T) %>%
  data.frame()
