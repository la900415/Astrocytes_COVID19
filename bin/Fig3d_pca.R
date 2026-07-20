library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(readr)
library(pcaMethods)
library(msdap)

setwd("Q:/data_analysis/COVID_Prot_USP/celulas/astrocytes_COVID19/study_DIANN/22_reanalysis_same_expSpecLib_61samp_without_isoforms_FINAL/msdap_results/2024-10-09_23-34-50_combinedCTRL_msqrob_FDR0.01_FINAL-BEST")

load("dataset.RData")
################################################################################
## Peptide intensity matrix for PCA
################################################################################

# peptide level data matrix, using filtered+normalized peptides across all groups
tibw_noexclude = dataset$peptides %>% 
  select(peptide_id, sample_id, intensity_all_group) %>%
  filter(!is.na(intensity_all_group)) %>%
  pivot_wider(id_cols = peptide_id, names_from = sample_id, values_from = intensity_all_group)
matrix_sample_intensities = msdap:::as_matrix_except_first_column(tibw_noexclude)

# compute PCA
PPCA = pcaMethods::pca(base::scale(t(matrix_sample_intensities), center=TRUE, scale=FALSE), method="ppca", nPcs = 10, seed = 123, maxIterations = 2000)
mat_pca = PPCA@scores # PCA coordinate matrix; rownames = sample_id and colnames = PCA dimension
pca_var = array(PPCA@R2, dimnames = list(colnames(PPCA@scores))) # variation per dimension

# if you only want a table with the PCA coordinates, just print the mat_pca variable and you're done at this point

# example plot code
dims = c(1,2) # plot dimensions. in this example, can be; 1:2, 2:3, c(1,3)
p = ggplot(data.frame(x=mat_pca[,dims[1]], y=mat_pca[,dims[2]], label=rownames(mat_pca)), 
           aes(x = x, y = y, label = label)) +
  geom_point() +
  geom_text() +
  labs(x = sprintf("dimension %d (%.1f%%)", dims[1], pca_var[dims[1]] * 100),
       y = sprintf("dimension %d (%.1f%%)", dims[2], pca_var[dims[2]] * 100) ) +
  theme_bw()
print(p)


# Metadata
pca_pep_df <-
  as.data.frame(mat_pca) %>%
  tibble::rownames_to_column("sample_id") %>%
  left_join(
    dataset$samples |> filter(exclude == FALSE),
    by = "sample_id"
  ) 

# final plot code
pep_pca <-
  ggplot(
    pca_pep_df,
    aes(
      PC1,
      PC2,
      colour = group
    )  ) +
  geom_hline(
    yintercept = 0,
    colour = "grey90",
    linewidth = 0.4 ) +
  geom_vline(
    xintercept = 0,
    colour = "grey90",
    linewidth = 0.4) +
  geom_point(
    size = 3.8,
    stroke = 1,
    alpha = .95
  ) +
  ggrepel::geom_text_repel(
    aes(label = shortname),
    size = 4,
    max.overlaps = Inf,
    box.padding = .35,
    point.padding = .25,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = c("Ctrl"="#7E6148B2", "Ga"="#E64B35B2","Gc"="#F39B7FB2", "Wa"="#3C5488B2", "Wc"="#8491B4B1")  ) +
  #scale_shape_manual(values = shape_values  ) +
  labs(
    x = sprintf(
      "PC1 (%.1f%%)",
      100 * pca_var[1]
    ),
    y = sprintf(
      "PC2 (%.1f%%)",
      100 * pca_var[2]
    )
  ) +
  coord_equal() +
  theme_bw(base_size = 14) 

pep_pca
