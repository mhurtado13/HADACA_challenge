
rm(list = ls())

mixes1_SDE5_pdac <- readRDS(file = "/home/benoitcl/HADACA/starting_kit_phase1/data/mixes1_SDE5_pdac.rds")
  
reference_pdac <- readRDS("/home/benoitcl/HADACA/starting_kit_phase1/data/reference_pdac.rds")
reference_pdac$ref_scRNA$ref_sc_peng$metadata[1:5,]
reference_pdac$ref_scRNA$ref_sc_baron$metadata[1:5,]
reference_pdac$ref_scRNA$ref_sc_raghavan$metadata[1:5,]
head(reference_pdac$ref_bulkRNA)

# Charger les bibliothèques nécessaires
library(Seurat)
library(dplyr)

# Créer l'objet Seurat à partir des données
# `counts_matrix` contient les données de comptage (raw counts)
# `metadata` contient les métadonnées associées

counts_matrix <- reference_pdac$ref_scRNA$ref_sc_peng$counts
metadata <- reference_pdac$ref_scRNA$ref_sc_peng$metadata

# Assurez-vous que les colonnes des comptes correspondent aux noms des cellules dans les métadonnées
all(colnames(counts_matrix) == rownames(metadata)) # Doit être TRUE

# Créer l'objet Seurat
seurat_object <- CreateSeuratObject(counts = counts_matrix, meta.data = metadata)

# Normalisation des données
seurat_object <- NormalizeData(seurat_object)

# Trouver les gènes variables
seurat_object <- FindVariableFeatures(seurat_object)

# Échelle des données (scaling)
seurat_object <- ScaleData(seurat_object)

# Comparaison différentielle pour un type cellulaire contre les autres (par exemple "immune")
immune_vs_all <- FindMarkers(seurat_object, 
                             ident.1 = "immune", 
                             group.by = "cell_type", 
                             min.pct = 0.25, # Pour inclure les gènes exprimés dans au moins 25% des cellules
                             test.use = "wilcox") # Test de Wilcoxon par défaut

immune_vs_all_filtered <- immune_vs_all %>%
  dplyr::filter(p_val_adj < 0.05)

# Comparaison d'un autre type cellulaire (par exemple "endo")
endo_vs_all <- FindMarkers(seurat_object, 
                           ident.1 = "endo", 
                           group.by = "cell_type", 
                           min.pct = 0.25, 
                           test.use = "wilcox")
endo_vs_all_filtered <- endo_vs_all %>%
  dplyr::filter(p_val_adj < 0.05)

fibro_vs_all <- FindMarkers(seurat_object, 
                           ident.1 = "fibro", 
                           group.by = "cell_type", 
                           min.pct = 0.25, 
                           test.use = "wilcox")
fibro_vs_all_filtered <- fibro_vs_all %>%
  dplyr::filter(p_val_adj < 0.05)

classic_vs_all <- FindMarkers(seurat_object, 
                            ident.1 = "classic", 
                            group.by = "cell_type", 
                            min.pct = 0.25, 
                            test.use = "wilcox")
classic_vs_all_filtered <- classic_vs_all %>%
  dplyr::filter(p_val_adj < 0.05)


basal_vs_all <- FindMarkers(seurat_object, 
                              ident.1 = "basal", 
                              group.by = "cell_type", 
                              min.pct = 0.25, 
                              test.use = "wilcox")
basal_vs_all_filtered <- basal_vs_all %>%
  dplyr::filter(p_val_adj < 0.05)

# Affichage des résultats
head(immune_vs_all_filtered)!
head(endo_vs_all_filtered)
head(fibro_vs_all_filtered)
head(classic_vs_all_filtered)
head(basal_vs_all_filtered)

gene_names <- rownames(immune_vs_all_filtered)
# Create a new data frame
immune_genes_df <- data.frame(
  gene_names = gene_names,
  cell_type = "immune"  # Assign "immune" to all rows
)
gene_names <- rownames(endo_vs_all_filtered)
# Create a new data frame
endo_genes_df <- data.frame(
  gene_names = gene_names,
  cell_type = "endo"  # Assign "immune" to all rows
)
gene_names <- rownames(classic_vs_all_filtered)
# Create a new data frame
classic_genes_df <- data.frame(
  gene_names = gene_names,
  cell_type = "classic"  # Assign "immune" to all rows
)
gene_names <- rownames(fibro_vs_all_filtered)
fibro_genes_df <- data.frame(
  gene_names = gene_names,
  cell_type = "fibro"  # Assign "immune" to all rows
)
gene_names <- rownames(basal_vs_all_filtered)
base_genes_df <- data.frame(
  gene_names = gene_names,
  cell_type = "basal"  # Assign "immune" to all rows
)

gene_list <- rbind(immune_genes_df, endo_genes_df, classic_genes_df, fibro_genes_df, base_genes_df)

computeMCP <- function(TPM_matrix, genes_path) {
  library(MCPcounter)
  genes <- read.table(paste0(genes_path, "/MCPcounter/MCPcounter-genes.txt"), sep = "\t", stringsAsFactors = FALSE, header = TRUE, colClasses = "character", check.names = FALSE)
  mcp <- MCPcounter.estimate(TPM_matrix, genes = genes, featuresType = "HUGO_symbols", probesets = NULL) %>%
    t()
  
  colnames(mcp) = paste0("MCP_", colnames(mcp))
  colnames(mcp) <- colnames(mcp) %>%
    str_replace_all(., " ", "_")xx
  
  return(mcp)
}
