
rm(list = ls())

mixes1_SDE5_pdac <- readRDS(file = "/home/benoitcl/HADACA/starting_kit_phase1/data/mixes1_SDE5_pdac.rds")
  
reference_pdac <- readRDS("/home/benoitcl/HADACA/starting_kit_phase1/data/reference_pdac.rds")

# Charger les bibliothèques nécessaires
library(Seurat)
library(dplyr)

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

gene_names <- rownames(immune_vs_all_filtered)
# Create a new data frame
immune_genes_df <- data.frame(
  gene_names = gene_names,
  cell_type = "immune",
  p_adj = immune_vs_all_filtered$p_val_adj # Assign "immune" to all rows
)
gene_names <- rownames(endo_vs_all_filtered)
# Create a new data frame
endo_genes_df <- data.frame(
  gene_names = gene_names,
  cell_type = "endo",
  p_adj = endo_vs_all_filtered$p_val_adj# Assign "immune" to all rows
)
gene_names <- rownames(classic_vs_all_filtered)
# Create a new data frame
classic_genes_df <- data.frame(
  gene_names = gene_names,
  cell_type = "classic",  # Assign "immune" to all rows
  p_adj = classic_vs_all_filtered$p_val_adj
)
gene_names <- rownames(fibro_vs_all_filtered)
fibro_genes_df <- data.frame(
  gene_names = gene_names,
  cell_type = "fibro",
  p_adj = fibro_vs_all_filtered$p_val_adj# Assign "immune" to all rows
)
gene_names <- rownames(basal_vs_all_filtered)
base_genes_df <- data.frame(
  gene_names = gene_names,
  cell_type = "basal",
  p_adj = basal_vs_all_filtered$p_val_adj# Assign "immune" to all rows
)

gene_list <- rbind(immune_genes_df, endo_genes_df, classic_genes_df, fibro_genes_df, base_genes_df)
colnames(gene_list) <- c("HUGO symbols", "Cell population")
write.table(gene_list, file = "peng_gene_list.txt", sep = ",", quote = FALSE, row.names = FALSE)

#gene_list <- read.table(file = "peng_gene_list.txt", sep = ",", header = TRUE)
gene_list <- gene_list %>%
  group_by(cell_type) %>%
  filter(row_number() <= 50)

bulk_genes <- rownames(mixes1_SDE5_pdac$mix_rna)

# Liste des gènes dans gene_list
marker_genes <- unique(gene_list$gene_names)

# Gènes communs
common_genes <- intersect(bulk_genes, marker_genes)

# Vérification
length(common_genes) # Combien de gènes sont communs
#length(marker_genes) # Combien de gènes dans gene_list
#setdiff(marker_genes, common_genes)

# Filtrer les données bulk RNA pour conserver uniquement les gènes communs
filtered_bulkRNA <- mixes1_SDE5_pdac$mix_rna[common_genes, ]

# Filtrer gene_list pour ne conserver que les gènes communs
filtered_gene_list <- gene_list[gene_list$gene_names %in% common_genes, ]

# Construire la matrice de signature en prenant la moyenne pour chaque type cellulaire


# signature_matrix <- sapply(unique(filtered_gene_list$cell_type), function(cell_type) {
#   #genes <- filtered_gene_list$gene_names[filtered_gene_list$cell_type == cell_type]
#   rowMeans(reference_pdac$ref_bulkRNA[filtered_gene_list$gene_names , drop = FALSE], na.rm = FALSE)  # Handling NAs in rowMeans
# })

#head(t(signature_matrix))
#signature_matrix <- t(signature_matrix)
#View(filtered_bulkRNA)
#dim(filtered_bulkRNA)
#dim(signature_matrix)
proportions <- apply(filtered_bulkRNA, 2, function(sample_expr) {
  fit <- nnls(reference_pdac$ref_bulkRNA[rownames(reference_pdac$ref_bulkRNA) %in% gene_list$gene_names,], sample_expr)
  coef(fit)
})

# Normalisation des proportions
proportions_df <- as.data.frame(t(proportions))
colnames(proportions_df) <- colnames(reference_pdac$ref_bulkRNA)
rownames(proportions_df) <- colnames(filtered_bulkRNA)
proportions_df <- proportions_df / rowSums(proportions_df)
proportions_df <- t(proportions_df) 
#print(proportions_df)
return(proportions_df)



