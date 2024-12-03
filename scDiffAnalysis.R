
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

# Affichage des résultats
head(immune_vs_all)!
head(endo_vs_all)
