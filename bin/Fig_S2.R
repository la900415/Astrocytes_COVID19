library(dplyr)
library(tidyr)
library(tidyverse)
library(protti)

# Defining directory with setwd()
load("dataset.RData") # from msdap output

########## Combine dataset with the annotation ###############
dataset$peptides <- dataset$peptides %>%  left_join(y = dataset$samples, by = "key_sample")
dataset$peptides$key_group.y <- NULL
dataset$peptides$sample_id.y <- NULL
dataset$peptides$detected_proteins <- NULL
dataset$peptides$detected_peptides <- NULL
dataset$peptides <- dataset$peptides %>% 
  rename(sample_id = sample_id.x) %>% 
  rename(key_group = key_group.x)

########## Loading function to generate visualization of MSDAP #############
ggplot_sample_detect_counts_barplots = function(samples, samples_colors_long) {
  tib = samples
  
  # if there are 'all peptide' counts (detect + MBR/quant-only), include these in plot tibble
  if("all_peptides" %in% colnames(samples)) {
    tib$quantified_peptides = tib$all_peptides - tib$detected_peptides
    tib$quantified_proteins = tib$all_proteins - tib$detected_proteins
  }
  
  # from wide to long format
  tib = tib %>%
    select(!!c("shortname", grep("(detected|quantified)_(peptides|proteins)", colnames(tib), value = T))) %>%
    pivot_longer(cols = -shortname, names_to = "type", values_to = "count")
  
  # flip levels/sorting because we use coord_flip() downstream in ggplot
  tib$shortname = factor(tib$shortname, levels = rev(samples$shortname))
  tib$pep_or_prot = ifelse(grepl("peptide", tib$type), "peptides", "proteins")
  tib$detect_or_quant = ifelse(grepl("detected", tib$type), "detect", "quant")
  
  # rename labels for plot clarity
  tib$type["detected_peptides" == tib$type] = "peptides: identified & quantified"
  tib$type["quantified_peptides" == tib$type] = "peptides: only quantified"
  tib$type["detected_proteins" == tib$type] = "proteins: identified & quantified"
  tib$type["quantified_proteins" == tib$type] = "proteins: only quantified"
  
  # color-coding
  # clr = c("detected_peptides"="#0570b0", "quantified_peptides"="#74a9cf", "detected_proteins"="#6a51a3", "quantified_proteins"="#9e9ac8")
  clr = c("peptides: identified & quantified"="#0570b0", "peptides: only quantified"="#74a9cf",
          "proteins: identified & quantified"="#6a51a3", "proteins: only quantified"="#9e9ac8")
  
  # coordinates for our customized dual-barplot
  tib = tib %>% group_by(pep_or_prot, shortname) %>% mutate(total=sum(count)) %>% arrange(shortname, type)
  tib_dual = tib %>% group_by(pep_or_prot) %>% mutate(count_scaled = count / max(total) * ifelse(pep_or_prot=="proteins", -1, 1),
                                                      count_scaled_max = total / max(total) * ifelse(pep_or_prot=="proteins", -1, 1))
  tib_dual$outlier_lowside = abs(tib_dual$count_scaled_max) < 0.3
  
  # custom x-axis labels, since this dimension not ranges from -1:1
  ticks_peptides = pretty(1:max(tib %>% filter(pep_or_prot=="peptides") %>% pull(total)), n = 5, min.n = 4)
  ticks_peptides_scaled = ticks_peptides / max(tib %>% filter(pep_or_prot=="peptides") %>% pull(total))
  ticks_proteins = pretty(1:max(tib %>% filter(pep_or_prot=="proteins") %>% pull(total)), n = 5, min.n = 4)
  ticks_proteins_scaled = ticks_proteins / max(tib %>% filter(pep_or_prot=="proteins") %>% pull(total))
  # combine
  ticks_coord = c(-1 * rev(ticks_proteins_scaled), ticks_peptides_scaled[-1])
  ticks_label = c(rev(ticks_proteins), ticks_peptides[-1])
  
  
  # plot code is somewhat convoluted by careful alignment of text labels. simplest QC plot of data as-is: ggplot(tib_dual, aes(x = shortname, y = count_scaled, fill = type)) + geom_bar(stat = "identity", position = position_stack(reverse = T)) + coord_flip()
  text_size = ifelse(nrow(tib) > 50, 2, 3.5)
  p = ggplot(tib_dual, aes(x = shortname, y = count_scaled, fill = type)) +
    geom_bar(stat = "identity", position = position_stack(reverse = T)) +
    geom_text(aes(label = count,
                  y = ifelse(detect_or_quant=="detect", ifelse(pep_or_prot=="proteins", -0.025, 0.025), count_scaled_max),
                  hjust = ifelse((type %in% c("peptides: identified & quantified", "proteins: only quantified") & !c(type=="proteins: only quantified" & outlier_lowside)) |
                                   (type %in% c("peptides: only quantified") & outlier_lowside), -0.1, 1.1)),
              colour = ifelse(tib_dual$outlier_lowside & tib_dual$detect_or_quant=="quant", "darkgrey", "white"),
              size = text_size, # 4 is default
              check_overlap = F) +
    scale_y_continuous(breaks=ticks_coord, labels=ticks_label, ) +
    scale_fill_manual(values = clr, guide = guide_legend(nrow = 2)) +
    coord_flip() +
    labs(x = "Samples", y = "Protein / peptide counts") +
    ggpubr::theme_pubr(base_size = 10) +
    theme(legend.position = "bottom", legend.title = element_blank(), legend.key.size = unit(0.75,"line"), legend.text = element_text(size=9),
          panel.grid = element_blank(), axis.line.y.right = element_blank(),
          axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1))
  
  if(nrow(tib) > 50)
    p = p + theme(axis.text.y.left = element_text(size=6))
  
  return(p)
  ###### some reference code for separate plots
  # tib_pep = tib %>% filter(pep_or_prot == "peptides") %>% droplevels()
  # tib_prot = tib %>% filter(pep_or_prot == "proteins") %>% droplevels()
  #
  # plot_text_lim = 0.1 * max(tib_pep$count)
  # p_pep = ggplot(tib_pep, aes(x = shortname, y = count, fill = type)) +
  #   geom_bar(stat = "identity", position = position_stack(reverse = T)) +
  #   geom_text(aes(label = count, hjust = ifelse(count>plot_text_lim, 1.25, 0)), position = position_stack(reverse = T), colour = "white", check_overlap = T) + # v1
  #   scale_fill_manual(values = clr) +
  #   coord_flip() +
  #   labs(x = "", y = "") +
  #   ggpubr::theme_pubr() +
  #   theme(legend.position = "bottom", legend.title = element_blank(), panel.grid = element_blank(), axis.line.y.right = element_blank())
  #
  # plot_text_lim = 0.1 * max(tib_prot$count)
  # p_prot = ggplot(tib_prot, aes(x = shortname, y = count, fill = type)) +
  #   geom_bar(stat = "identity", position = position_stack(reverse = T)) +
  #   geom_text(aes(label = count, hjust = ifelse(count>plot_text_lim, 1.25, 0)), position = position_stack(reverse = T), colour = "white", check_overlap = T) + # v1
  #   scale_fill_manual(values = clr) +
  #   coord_flip() +
  #   labs(x = "", y = "") +
  #   ggpubr::theme_pubr() +
  #   theme(legend.position = "bottom", legend.title = element_blank(), panel.grid = element_blank(), axis.line.y.right = element_blank())
}

########## Figure S2A. Barplots of IDs from MSDAP ##########################
ggplot_sample_detect_counts_barplots(
  dataset$samples %>% filter(group != "pool"), 
  samples_colors_long=c("peptides: identified & quantified"="#0570b0", 
                        "peptides: only quantified"="#74a9cf",
                        "proteins: identified & quantified"="#6a51a3", 
                        "proteins: only quantified"="#9e9ac8")) +
  labs(title = "PGs & peptides identification") +
  theme(plot.title = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 12) )

########## Figure S2F. Pears corr using protti ####################
qc_sample_correlation(
  data = dataset$peptides %>% filter(group != "pool"),
  #dataset$protein_level$global_tidy,
  sample = shortname,
  grouping = peptide_id,
  #protein_id,
  intensity_log2 = intensity_all_group,
  condition = group,
  digestion = age,
  run_order = virus_strain,
  method = "pearson",
  interactive = F
)
