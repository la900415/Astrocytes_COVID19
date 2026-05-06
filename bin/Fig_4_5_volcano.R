######################## Fig 3. Volcano plots #########################################
setwd("C:/Users/Laura/Downloads/astrocitos_COVID19/2nd_reanalysis/msdap_results/2024-10-09_23-34-50_combinedCTRL_msqrob_FDR0.01")
load("dataset.RData")
dataset$de_proteins <- dataset$de_proteins %>% left_join(dataset$proteins, by=c("protein_id"="protein_id"))
########################################################################################

                         # Volcano plot for each contrast
cbbPalette <- c("#7E6148B2", "#E64B35B2", "#F39B7FB2", "#3C5488B2", "#8491B4B1") #palette for Ctrl, Ga, Gc, Wa, Wc groups
                                                       "blue2", "cornflowerblue" # other colors for Wuhan 

###### Gc_vs_Ga ########################################################################
tib = dataset$de_proteins %>% dplyr::filter(contrast=="contrast: Ga vs Gc # condition_variable: group") %>%
  filter(dea_algorithm == "msqrob") %>%
  drop_na(foldchange.log2, pvalue, qvalue) %>%
  # minlog10 conversion must be performed within this loop! scales zero's to max-value, ONLY valid within same statistical test
  mutate(minlog10qval = -log10(qvalue),
         foldchange.log2_abs = abs(foldchange.log2)) %>%
  # order by p-value so we draw the most significant proteins last (thus their symbols/PCH are on top)
  arrange(desc(minlog10qval)) %>%
  # reduce tibble size
  select(contrast, protein_id, gene_symbols, 
         signif_threshold_log2fc, signif_threshold_qvalue,
         foldchange.log2, foldchange.log2_abs, minlog10qval, signif)

# which proteins should get a text label?
tib$flag_plot_label = FALSE
tib$flag_plot_label[1:25] = TRUE

# classify up/down regulated
tib$updown = ifelse(tib$foldchange.log2 < 0, "down-regulated", "up-regulated")
tib$updown[tib$signif != TRUE] = "not significant"

tib %>% 
  ggplot(aes(foldchange.log2, minlog10qval, colour = updown, fill = updown, shape = updown, label=gene_symbols) ) +
  geom_point(alpha = 0.5, size=3) + # Alpha sets the transparency of the points
  # Add dotted lines to indicate the threshold, semi-transparent
  geom_hline(yintercept = -log10(tib$signif_threshold_qvalue),      colour = "darkgrey", linetype = "dashed", linewidth = 1) + 
  geom_vline(xintercept =   tib$signif_threshold_log2fc,  colour = "darkgrey", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = -(tib$signif_threshold_log2fc), colour = "darkgrey", linetype = "dashed", linewidth = 1) +
  # Set the colour of the points
  #scale_colour_manual(values = c("up"="#d55e00aa", "down"="#56b4e9aa", "unchanged"="black")) +
  scale_discrete_manual(aesthetics = c("colour", "fill"), values = c("up-regulated" = "#F39B7FB2", "down-regulated" = "#E64B35B2", "not significant" = "#22222299") ) +
  scale_shape_manual( values = c("up-regulated" = 24, "down-regulated" = 25, "not significant" = 19) ) +
  xlab("log2 fold-change (Gc/Ga)") + ylab("-log10 FDR adj. p-value") + ggtitle("Gc_vs_Ga") +
  theme_bw() +
  scale_y_continuous(limits=c(0, 2.5)) +
  theme(plot.title = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 14)) +
  ggrepel::geom_text_repel(data = tib %>% filter(flag_plot_label == TRUE), alpha=1, # vjust = 0.6,
                           show.legend = FALSE, size = 3, segment.alpha = .3, min.segment.length = 0, 
                           na.rm = TRUE, max.time = 1, max.iter = 1e5, max.overlaps = Inf, point.padding = 0, 
                           box.padding = 0.2, seed = 123) + # min.segment.length = unit(0.2, 'lines')
  annotate("text", x = -2.4, y = 2.4, label = "Ga (0)", hjust = 0.8, vjust = 0, size = 6, fontface="bold", colour="#E64B35B2") +
  annotate("text", x = 1.5, y = 2.4, label = "Gc (1)",   hjust = 0.5, vjust = 0, size = 6, fontface="bold", colour="#F39B7FB2") +
  #annotate("text", x = -3, y = 2.2, label = "",   hjust = 0.5, vjust = 0, size = 5, fontface="bold", colour="#d55e00aa") +
  #annotate("text", x = 2.5, y = 2.2, label = "",   hjust = 0.5, vjust = 0, size = 5, fontface="bold", colour="#56b4e9aa") +
  annotate("text", x = 0, y = 2.4, label = "(4562)",   hjust = 0.5, vjust = 0, size = 6, fontface="bold", colour="#22222299") 


###### Wc_vs_Wa ########################################################################
dataset$de_proteins <- dataset$de_proteins %>% left_join(dataset$proteins, by=c("protein_id"="protein_id"))

tib = dataset$de_proteins %>% dplyr::filter(contrast=="contrast: Wa vs Wc # condition_variable: group") %>%
  filter(dea_algorithm == "msqrob") %>%
  drop_na(foldchange.log2, pvalue, qvalue) %>%
  # minlog10 conversion must be performed within this loop! scales zero's to max-value, ONLY valid within same statistical test
  mutate(minlog10qval = -log10(qvalue),
         foldchange.log2_abs = abs(foldchange.log2)) %>%
  # order by p-value so we draw the most significant proteins last (thus their symbols/PCH are on top)
  arrange(desc(minlog10qval)) %>%
  # reduce tibble size
  select(contrast, protein_id, gene_symbols, 
         signif_threshold_log2fc, signif_threshold_qvalue,
         foldchange.log2, foldchange.log2_abs, minlog10qval, signif)

# which proteins should get a text label?
tib$flag_plot_label = FALSE
tib$flag_plot_label[1:25] = TRUE

# classify up/down regulated
tib$updown = ifelse(tib$foldchange.log2 < 0, "down-regulated", "up-regulated")
tib$updown[tib$signif != TRUE] = "not significant"

tib %>% 
  ggplot(aes(foldchange.log2, minlog10qval, colour = updown, fill = updown, shape = updown, label=gene_symbols) ) +
  geom_point(alpha = 0.5, size=3) + # Alpha sets the transparency of the points
  # Add dotted lines to indicate the threshold, semi-transparent
  geom_hline(yintercept = -log10(tib$signif_threshold_qvalue),colour = "darkgrey", linetype = "dashed", linewidth = 1) + 
  geom_vline(xintercept = tib$signif_threshold_log2fc, colour = "darkgrey", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = -(tib$signif_threshold_log2fc), colour = "darkgrey", linetype = "dashed", linewidth = 1) +
  # Set the colour of the points
  #scale_colour_manual(values = c("up"="#d55e00aa", "down"="#56b4e9aa", "unchanged"="black")) +
  scale_discrete_manual(aesthetics = c("colour", "fill"), values = c("up-regulated" = "cornflowerblue", "down-regulated" = "blue2", "not significant" = "#22222299") ) +
  scale_shape_manual( values = c("up-regulated" = 24, "down-regulated" = 25, "not significant" = 19) ) +
  xlab("log2 fold-change (Wc/Wa)") + ylab("-log10 FDR adj. p-value") + ggtitle("Wc_vs_Wa") +
  theme_bw() +
  scale_y_continuous(limits=c(0, 5)) +
  theme(plot.title = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 14)) +
  ggrepel::geom_text_repel(data = tib %>% filter(flag_plot_label == TRUE), alpha=1, # vjust = 0.6,
                           show.legend = FALSE, size = 3, segment.alpha = .3, min.segment.length = 0, 
                           na.rm = TRUE, max.time = 1, max.iter = 1e5, max.overlaps = Inf, point.padding = 0, 
                           box.padding = 0.2, seed = 123) + # min.segment.length = unit(0.2, 'lines')
  annotate("text", x = -4, y = 4.8, label = "Wa (0)", hjust = 0.8, vjust = 0, size = 6, fontface="bold", colour="blue2") +
  annotate("text", x = 1.5, y = 4.8, label = "Wc (1)",   hjust = 0.5, vjust = 0, size = 6, fontface="bold", colour="cornflowerblue") +
  #annotate("text", x = -3, y = 2.2, label = "",   hjust = 0.5, vjust = 0, size = 5, fontface="bold", colour="#d55e00aa") +
  #annotate("text", x = 2.5, y = 2.2, label = "",   hjust = 0.5, vjust = 0, size = 5, fontface="bold", colour="#56b4e9aa") +
  annotate("text", x = 0, y = 4.8, label = "(4962)",   hjust = 0.5, vjust = 0, size = 5, fontface="bold", colour="#22222299") 

###### Wa_vs_Ga ########################################################################
dataset$de_proteins <- dataset$de_proteins %>% left_join(dataset$proteins, by=c("protein_id"="protein_id"))

tib = dataset$de_proteins %>% dplyr::filter(contrast=="contrast: Ga vs Wa # condition_variable: group") %>%
  filter(dea_algorithm == "msqrob") %>%
  drop_na(foldchange.log2, pvalue, qvalue) %>%
  # minlog10 conversion must be performed within this loop! scales zero's to max-value, ONLY valid within same statistical test
  mutate(minlog10qval = -log10(qvalue),
         foldchange.log2_abs = abs(foldchange.log2)) %>%
  # order by p-value so we draw the most significant proteins last (thus their symbols/PCH are on top)
  arrange(desc(minlog10qval)) %>%
  # reduce tibble size
  select(contrast, protein_id, gene_symbols, 
         signif_threshold_log2fc, signif_threshold_qvalue,
         foldchange.log2, foldchange.log2_abs, minlog10qval, signif)

# which proteins should get a text label?
tib$flag_plot_label = FALSE
tib$flag_plot_label[1:25] = TRUE

# classify up/down regulated
tib$updown = ifelse(tib$foldchange.log2 < 0, "down-regulated", "up-regulated")
tib$updown[tib$signif != TRUE] = "not significant"

tib %>% 
  ggplot(aes(foldchange.log2, minlog10qval, colour = updown, fill = updown, shape = updown, label=gene_symbols) ) +
  geom_point(alpha = 0.5, size=3) + # Alpha sets the transparency of the points
  # Add dotted lines to indicate the threshold, semi-transparent
  geom_hline(yintercept = -log10(tib$signif_threshold_qvalue),colour = "darkgrey", linetype = "dashed", linewidth = 1) + 
  geom_vline(xintercept = tib$signif_threshold_log2fc, colour = "darkgrey", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = -(tib$signif_threshold_log2fc), colour = "darkgrey", linetype = "dashed", linewidth = 1) +
  # Set the colour of the points
  #scale_colour_manual(values = c("up"="#d55e00aa", "down"="#56b4e9aa", "unchanged"="black")) +
  scale_discrete_manual(aesthetics = c("colour", "fill"), values = c("up-regulated" = "blue2", "down-regulated" = "#E64B35B2", "not significant" = "#22222299") ) +
  scale_shape_manual( values = c("up-regulated" = 24, "down-regulated" = 25, "not significant" = 19) ) +
  xlab("log2 fold-change (Wa/Ga)") + ylab("-log10 FDR adj. p-value") + ggtitle("Wa_vs_Ga") +
  theme_bw() +
  scale_y_continuous(limits=c(0, 7.2)) +
  theme(plot.title = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 14)) +
  ggrepel::geom_text_repel(data = tib %>% filter(flag_plot_label == TRUE), alpha=1, # vjust = 0.6,
                           show.legend = FALSE, size = 3, segment.alpha = .3, min.segment.length = 0, 
                           na.rm = TRUE, max.time = 1, max.iter = 1e5, max.overlaps = Inf, point.padding = 0, 
                           box.padding = 0.2, seed = 123) + # min.segment.length = unit(0.2, 'lines')
  annotate("text", x = -1.5, y = 7, label = "Ga (0)", hjust = 0.8, vjust = 0, size = 6, fontface="bold", colour="#E64B35B2") +
  annotate("text", x = 3, y = 7, label = "Wa (1)",   hjust = 0.5, vjust = 0, size = 6, fontface="bold", colour="blue2") +
  #annotate("text", x = -3, y = 2.2, label = "",   hjust = 0.5, vjust = 0, size = 5, fontface="bold", colour="#d55e00aa") +
  #annotate("text", x = 2.5, y = 2.2, label = "",   hjust = 0.5, vjust = 0, size = 5, fontface="bold", colour="#56b4e9aa") +
  annotate("text", x = 0, y = 7, label = "(4368)",   hjust = 0.5, vjust = 0, size = 6, fontface="bold", colour="#22222299") 

###### Wc_vs_Gc ########################################################################
dataset$de_proteins <- dataset$de_proteins %>% left_join(dataset$proteins, by=c("protein_id"="protein_id"))

tib = dataset$de_proteins %>% dplyr::filter(contrast=="contrast: Gc vs Wc # condition_variable: group") %>%
  filter(dea_algorithm == "msqrob") %>%
  drop_na(foldchange.log2, pvalue, qvalue) %>%
  # minlog10 conversion must be performed within this loop! scales zero's to max-value, ONLY valid within same statistical test
  mutate(minlog10qval = -log10(qvalue),
         foldchange.log2_abs = abs(foldchange.log2)) %>%
  # order by p-value so we draw the most significant proteins last (thus their symbols/PCH are on top)
  arrange(desc(minlog10qval)) %>%
  # reduce tibble size
  select(contrast, protein_id, gene_symbols, 
         signif_threshold_log2fc, signif_threshold_qvalue,
         foldchange.log2, foldchange.log2_abs, minlog10qval, signif)

# which proteins should get a text label?
tib$flag_plot_label = FALSE
tib$flag_plot_label[1:25] = TRUE

# classify up/down regulated
tib$updown = ifelse(tib$foldchange.log2 < 0, "down-regulated", "up-regulated")
tib$updown[tib$signif != TRUE] = "not significant"

tib %>% 
  ggplot(aes(foldchange.log2, minlog10qval, colour = updown, fill = updown, shape = updown, label=gene_symbols) ) +
  geom_point(alpha = 0.5, size=3) + # Alpha sets the transparency of the points
  # Add dotted lines to indicate the threshold, semi-transparent
  geom_hline(yintercept = -log10(tib$signif_threshold_qvalue),colour = "darkgrey", linetype = "dashed", linewidth = 1) + 
  geom_vline(xintercept = tib$signif_threshold_log2fc, colour = "darkgrey", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = -(tib$signif_threshold_log2fc), colour = "darkgrey", linetype = "dashed", linewidth = 1) +
  # Set the colour of the points
  #scale_colour_manual(values = c("up"="#d55e00aa", "down"="#56b4e9aa", "unchanged"="black")) +
  scale_discrete_manual(aesthetics = c("colour", "fill"), values = c("up-regulated" = "cornflowerblue", "down-regulated" = "#F39B7FB2", "not significant" = "#22222299") ) +
  scale_shape_manual( values = c("up-regulated" = 24, "down-regulated" = 25, "not significant" = 19) ) +
  xlab("log2 fold-change (Wc/Gc)") + ylab("-log10 FDR adj. p-value") + ggtitle("Wc_vs_Gc") +
  theme_bw() +
  scale_y_continuous(limits=c(0, 10)) +
  theme(plot.title = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 14)) +
  ggrepel::geom_text_repel(data = tib %>% filter(flag_plot_label == TRUE), alpha=1, # vjust = 0.6,
                           show.legend = FALSE, size = 3, segment.alpha = .3, min.segment.length = 0, 
                           na.rm = TRUE, max.time = 1, max.iter = 1e5, max.overlaps = Inf, point.padding = 0, 
                           box.padding = 0.2, seed = 123) + # min.segment.length = unit(0.2, 'lines')
  annotate("text", x = -1.5, y = 9.8, label = "Gc (56)", hjust = 0.8, vjust = 0, size = 6, fontface="bold", colour="#F39B7FB2") +
  annotate("text", x = 1.5, y = 9.8, label = "Wc (9)",   hjust = 0.5, vjust = 0, size = 6, fontface="bold", colour="cornflowerblue") +
  #annotate("text", x = -3, y = 2.2, label = "",   hjust = 0.5, vjust = 0, size = 5, fontface="bold", colour="#d55e00aa") +
  #annotate("text", x = 2.5, y = 2.2, label = "",   hjust = 0.5, vjust = 0, size = 5, fontface="bold", colour="#56b4e9aa") +
  annotate("text", x = 0, y = 9.8, label = "(5144)",   hjust = 0.5, vjust = 0, size = 6, fontface="bold", colour="#22222299") 


###### Wc_vs_Ctrl ########################################################################
dataset$de_proteins <- dataset$de_proteins %>% left_join(dataset$proteins, by=c("protein_id"="protein_id"))

tib = dataset$de_proteins %>% dplyr::filter(contrast=="contrast: Ctrl vs Wc # condition_variable: group") %>%
  filter(dea_algorithm == "msqrob") %>%
  drop_na(foldchange.log2, pvalue, qvalue) %>%
  # minlog10 conversion must be performed within this loop! scales zero's to max-value, ONLY valid within same statistical test
  mutate(minlog10qval = -log10(qvalue),
         foldchange.log2_abs = abs(foldchange.log2)) %>%
  # order by p-value so we draw the most significant proteins last (thus their symbols/PCH are on top)
  arrange(desc(minlog10qval)) %>%
  # reduce tibble size
  select(contrast, protein_id, gene_symbols, 
         signif_threshold_log2fc, signif_threshold_qvalue,
         foldchange.log2, foldchange.log2_abs, minlog10qval, signif)

# which proteins should get a text label?
tib$flag_plot_label = FALSE
tib$flag_plot_label[1:25] = TRUE

# classify up/down regulated
tib$updown = ifelse(tib$foldchange.log2 < 0, "down-regulated", "up-regulated")
tib$updown[tib$signif != TRUE] = "not significant"

tib %>% 
  ggplot(aes(foldchange.log2, minlog10qval, colour = updown, fill = updown, shape = updown, label=gene_symbols) ) +
  geom_point(alpha = 0.5, size=3) + # Alpha sets the transparency of the points
  # Add dotted lines to indicate the threshold, semi-transparent
  geom_hline(yintercept = -log10(tib$signif_threshold_qvalue),colour = "darkgrey", linetype = "dashed", linewidth = 1) + 
  geom_vline(xintercept = tib$signif_threshold_log2fc, colour = "darkgrey", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = -(tib$signif_threshold_log2fc), colour = "darkgrey", linetype = "dashed", linewidth = 1) +
  # Set the colour of the points
  #scale_colour_manual(values = c("up"="#d55e00aa", "down"="#56b4e9aa", "unchanged"="black")) +
  scale_discrete_manual(aesthetics = c("colour", "fill"), values = c("up-regulated" = "cornflowerblue", "down-regulated" = "#7E6148B2", "not significant" = "#22222299") ) +
  scale_shape_manual( values = c("up-regulated" = 24, "down-regulated" = 25, "not significant" = 19) ) +
  xlab("log2 fold-change (Wc/Ctrl)") + ylab("-log10 FDR adj. p-value") + ggtitle("Wc_vs_Ctrl") +
  theme_bw() +
  scale_y_continuous(limits=c(0, 55)) +
  theme(plot.title = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 14)) +
  ggrepel::geom_text_repel(data = tib %>% filter(flag_plot_label == TRUE), alpha=1, # vjust = 0.6,
                           show.legend = FALSE, size = 3, segment.alpha = .3, min.segment.length = 0, 
                           na.rm = TRUE, max.time = 1, max.iter = 1e5, max.overlaps = Inf, point.padding = 0, 
                           box.padding = 0.2, seed = 123) + # min.segment.length = unit(0.2, 'lines')
  annotate("text", x = -2.5, y = 52.5, label = "Ctrl (532)", hjust = 0.8, vjust = 0, size = 6, fontface="bold", colour="#7E6148B2") +
  annotate("text", x = 3.5, y = 52.5, label = "Wc (737)",   hjust = 0.5, vjust = 0, size = 6, fontface="bold", colour="cornflowerblue") +
  #annotate("text", x = -3, y = 2.2, label = "",   hjust = 0.5, vjust = 0, size = 5, fontface="bold", colour="#d55e00aa") +
  #annotate("text", x = 2.5, y = 2.2, label = "",   hjust = 0.5, vjust = 0, size = 5, fontface="bold", colour="#56b4e9aa") +
  annotate("text", x = 0, y = 52.5, label = "(2166)",   hjust = 0.5, vjust = 0, size = 6, fontface="bold", colour="black") 

###### Wa_vs_Ctrl ########################################################################
dataset$de_proteins <- dataset$de_proteins %>% left_join(dataset$proteins, by=c("protein_id"="protein_id"))

tib = dataset$de_proteins %>% dplyr::filter(contrast=="contrast: Ctrl vs Wa # condition_variable: group") %>%
  filter(dea_algorithm == "msqrob") %>%
  drop_na(foldchange.log2, pvalue, qvalue) %>%
  # minlog10 conversion must be performed within this loop! scales zero's to max-value, ONLY valid within same statistical test
  mutate(minlog10qval = -log10(qvalue),
         foldchange.log2_abs = abs(foldchange.log2)) %>%
  # order by p-value so we draw the most significant proteins last (thus their symbols/PCH are on top)
  arrange(desc(minlog10qval)) %>%
  # reduce tibble size
  select(contrast, protein_id, gene_symbols, 
         signif_threshold_log2fc, signif_threshold_qvalue,
         foldchange.log2, foldchange.log2_abs, minlog10qval, signif)

# which proteins should get a text label?
tib$flag_plot_label = FALSE
tib$flag_plot_label[1:25] = TRUE

# classify up/down regulated
tib$updown = ifelse(tib$foldchange.log2 < 0, "down-regulated", "up-regulated")
tib$updown[tib$signif != TRUE] = "not significant"

tib %>% 
  ggplot(aes(foldchange.log2, minlog10qval, colour = updown, fill = updown, shape = updown, label=gene_symbols) ) +
  geom_point(alpha = 0.5, size=3) + # Alpha sets the transparency of the points
  # Add dotted lines to indicate the threshold, semi-transparent
  geom_hline(yintercept = -log10(tib$signif_threshold_qvalue),colour = "darkgrey", linetype = "dashed", linewidth = 1) + 
  geom_vline(xintercept = tib$signif_threshold_log2fc, colour = "darkgrey", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = -(tib$signif_threshold_log2fc), colour = "darkgrey", linetype = "dashed", linewidth = 1) +
  # Set the colour of the points
  #scale_colour_manual(values = c("up"="#d55e00aa", "down"="#56b4e9aa", "unchanged"="black")) +
  scale_discrete_manual(aesthetics = c("colour", "fill"), values = c("up-regulated" = "blue3", "down-regulated" = "#7E6148B2", "not significant" = "#22222299") ) +
  scale_shape_manual( values = c("up-regulated" = 24, "down-regulated" = 25, "not significant" = 19) ) +
  xlab("log2 fold-change (Wa/Ctrl)") + ylab("-log10 FDR adj. p-value") + ggtitle("Wa_vs_Ctrl") +
  theme_bw() +
  scale_y_continuous(limits=c(0, 37)) +
  theme(plot.title = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 14)) +
  ggrepel::geom_text_repel(data = tib %>% filter(flag_plot_label == TRUE), alpha=1, # vjust = 0.6,
                           show.legend = FALSE, size = 3, segment.alpha = .3, min.segment.length = 0, 
                           na.rm = TRUE, max.time = 1, max.iter = 1e5, max.overlaps = Inf, point.padding = 0, 
                           box.padding = 0.2, seed = 123) + # min.segment.length = unit(0.2, 'lines')
  annotate("text", x = -1.5, y=35, label = "Ctrl (419)", hjust = 0.8, vjust = 0, size = 6, fontface="bold", colour="#7E6148B2") +
  annotate("text", x = 2.7, y=35, label = "Wa (546)",   hjust = 0.5, vjust = 0, size = 6, fontface="bold", colour="blue3") +
  #annotate("text", x = -3, y = 2.2, label = "",   hjust = 0.5, vjust = 0, size = 5, fontface="bold", colour="#d55e00aa") +
  #annotate("text", x = 2.5, y = 2.2, label = "",   hjust = 0.5, vjust = 0, size = 5, fontface="bold", colour="#56b4e9aa") +
  annotate("text", x = 0, y=35, label = "(2162)",   hjust = 0.5, vjust = 0, size = 6, fontface="bold", colour="black") 

###### Ga_vs_Ctrl ########################################################################
dataset$de_proteins <- dataset$de_proteins %>% left_join(dataset$proteins, by=c("protein_id"="protein_id"))

tib = dataset$de_proteins %>% dplyr::filter(contrast=="contrast: Ctrl vs Ga # condition_variable: group") %>%
  filter(dea_algorithm == "msqrob") %>%
  drop_na(foldchange.log2, pvalue, qvalue) %>%
  # minlog10 conversion must be performed within this loop! scales zero's to max-value, ONLY valid within same statistical test
  mutate(minlog10qval = -log10(qvalue),
         foldchange.log2_abs = abs(foldchange.log2)) %>%
  # order by p-value so we draw the most significant proteins last (thus their symbols/PCH are on top)
  arrange(desc(minlog10qval)) %>%
  # reduce tibble size
  select(contrast, protein_id, gene_symbols, 
         signif_threshold_log2fc, signif_threshold_qvalue,
         foldchange.log2, foldchange.log2_abs, minlog10qval, signif)

# which proteins should get a text label?
tib$flag_plot_label = FALSE
tib$flag_plot_label[1:25] = TRUE

# classify up/down regulated
tib$updown = ifelse(tib$foldchange.log2 < 0, "down-regulated", "up-regulated")
tib$updown[tib$signif != TRUE] = "not significant"

tib %>% 
  ggplot(aes(foldchange.log2, minlog10qval, colour = updown, fill = updown, shape = updown, label=gene_symbols) ) +
  geom_point(alpha = 0.5, size=3) + # Alpha sets the transparency of the points
  # Add dotted lines to indicate the threshold, semi-transparent
  geom_hline(yintercept = -log10(tib$signif_threshold_qvalue),colour = "darkgrey", linetype = "dashed", linewidth = 1) + 
  geom_vline(xintercept = tib$signif_threshold_log2fc, colour = "darkgrey", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = -(tib$signif_threshold_log2fc), colour = "darkgrey", linetype = "dashed", linewidth = 1) +
  # Set the colour of the points
  #scale_colour_manual(values = c("up"="#d55e00aa", "down"="#56b4e9aa", "unchanged"="black")) +
  scale_discrete_manual(aesthetics = c("colour", "fill"), values = c("up-regulated" = "#E64B35B2", "down-regulated" = "#7E6148B2", "not significant" = "#22222299") ) +
  scale_shape_manual( values = c("up-regulated" = 24, "down-regulated" = 25, "not significant" = 19) ) +
  xlab("log2 fold-change (Ga/Ctrl)") + ylab("-log10 FDR adj. p-value") + ggtitle("Ga_vs_Ctrl") +
  theme_bw() +
  scale_y_continuous(limits=c(0, 40)) +
  theme(plot.title = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 14)) +
  ggrepel::geom_text_repel(data = tib %>% filter(flag_plot_label == TRUE), alpha=1, # vjust = 0.6,
                           show.legend = FALSE, size = 3, segment.alpha = .3, min.segment.length = 0, 
                           na.rm = TRUE, max.time = 1, max.iter = 1e5, max.overlaps = Inf, point.padding = 0, 
                           box.padding = 0.2, seed = 123) + # min.segment.length = unit(0.2, 'lines')
  annotate("text", x = -1.5, y=38, label = "Ctrl (447)", hjust = 0.8, vjust = 0, size = 6, fontface="bold", colour="#7E6148B2") +
  annotate("text", x = 3, y=38, label = "Ga (582)",   hjust = 0.5, vjust = 0, size = 6, fontface="bold", colour="#E64B35B2") +
  #annotate("text", x = -3, y = 2.2, label = "",   hjust = 0.5, vjust = 0, size = 5, fontface="bold", colour="#d55e00aa") +
  #annotate("text", x = 2.5, y = 2.2, label = "",   hjust = 0.5, vjust = 0, size = 5, fontface="bold", colour="#56b4e9aa") +
  annotate("text", x = 0, y=38, label = "(1863)",   hjust = 0.5, vjust = 0, size = 6, fontface="bold", colour="black") 

###### Gc_vs_Ctrl ########################################################################
dataset$de_proteins <- dataset$de_proteins %>% left_join(dataset$proteins, by=c("protein_id"="protein_id"))

tib = dataset$de_proteins %>% dplyr::filter(contrast=="contrast: Ctrl vs Gc # condition_variable: group") %>%
  filter(dea_algorithm == "msqrob") %>%
  drop_na(foldchange.log2, pvalue, qvalue) %>%
  # minlog10 conversion must be performed within this loop! scales zero's to max-value, ONLY valid within same statistical test
  mutate(minlog10qval = -log10(qvalue),
         foldchange.log2_abs = abs(foldchange.log2)) %>%
  # order by p-value so we draw the most significant proteins last (thus their symbols/PCH are on top)
  arrange(desc(minlog10qval)) %>%
  # reduce tibble size
  select(contrast, protein_id, gene_symbols, 
         signif_threshold_log2fc, signif_threshold_qvalue,
         foldchange.log2, foldchange.log2_abs, minlog10qval, signif)

# which proteins should get a text label?
tib$flag_plot_label = FALSE
tib$flag_plot_label[1:25] = TRUE

# classify up/down regulated
tib$updown = ifelse(tib$foldchange.log2 < 0, "down-regulated", "up-regulated")
tib$updown[tib$signif != TRUE] = "not significant"

tib %>% 
  ggplot(aes(foldchange.log2, minlog10qval, colour = updown, fill = updown, shape = updown, label=gene_symbols) ) +
  geom_point(alpha = 0.5, size=3) + # Alpha sets the transparency of the points
  # Add dotted lines to indicate the threshold, semi-transparent
  geom_hline(yintercept = -log10(tib$signif_threshold_qvalue),colour = "darkgrey", linetype = "dashed", linewidth = 1) + 
  geom_vline(xintercept = tib$signif_threshold_log2fc, colour = "darkgrey", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = -(tib$signif_threshold_log2fc), colour = "darkgrey", linetype = "dashed", linewidth = 1) +
  # Set the colour of the points
  #scale_colour_manual(values = c("up"="#d55e00aa", "down"="#56b4e9aa", "unchanged"="black")) +
  scale_discrete_manual(aesthetics = c("colour", "fill"), values = c("up-regulated" = "#F39B7FB2", "down-regulated" = "#7E6148B2", "not significant" = "#22222299") ) +
  scale_shape_manual( values = c("up-regulated" = 24, "down-regulated" = 25, "not significant" = 19) ) +
  xlab("log2 fold-change (Gc/Ctrl)") + ylab("-log10 FDR adj. p-value") + ggtitle("Gc_vs_Ctrl") +
  theme_bw() +
  scale_y_continuous(limits=c(0, 40)) +
  theme(plot.title = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 14)) +
  ggrepel::geom_text_repel(data = tib %>% filter(flag_plot_label == TRUE), alpha=1, # vjust = 0.6,
                           show.legend = FALSE, size = 3, segment.alpha = .3, min.segment.length = 0, 
                           na.rm = TRUE, max.time = 1, max.iter = 1e5, max.overlaps = Inf, point.padding = 0, 
                           box.padding = 0.2, seed = 123) + # min.segment.length = unit(0.2, 'lines')
  annotate("text", x = -1.5, y=38, label = "Ctrl (617)", hjust = 0.8, vjust = 0, size = 6, fontface="bold", colour="#7E6148B2") +
  annotate("text", x = 3, y=38, label = "Gc (601)",   hjust = 0.5, vjust = 0, size = 6, fontface="bold", colour="#F39B7FB2") +
  #annotate("text", x = -3, y = 2.2, label = "",   hjust = 0.5, vjust = 0, size = 5, fontface="bold", colour="#d55e00aa") +
  #annotate("text", x = 2.5, y = 2.2, label = "",   hjust = 0.5, vjust = 0, size = 5, fontface="bold", colour="#56b4e9aa") +
  annotate("text", x = 0, y=38, label = "(2057)",   hjust = 0.5, vjust = 0, size = 6, fontface="bold", colour="black") 
