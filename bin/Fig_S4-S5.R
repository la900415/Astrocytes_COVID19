##### Markers #######
astrocyte_marker <- read_excel("Fig2_rank_plot/marker_astrocytes_Jurga2021.xlsx", sheet = "table2_astrocytes", skip = 1)
astrocyte_marker$specif_astrocyte_markers <-  as.logical(astrocyte_marker$specif_astrocyte_markers)
astrocyte_marker$general_astrocyte_markers <- as.logical(astrocyte_marker$general_astrocyte_markers)
#astrocyte_marker <- separate(astrocyte_marker, Unique_Gene_Marker, into = c("Unique_Gene_Marker", "U2"), sep = ";")
astrogliosis_marker <- read_excel("Fig2_rank_plot/marker_astrocytes_Jurga2021.xlsx", sheet = "table3_astrogliosis", skip = 1)
astrogliosis_marker$markers_reactive_astrogliosis <- as.logical(astrogliosis_marker$markers_reactive_astrogliosis)

##### Wc_vs_Ctrl #######
Wc_vs_Ctrl <- read_delim("protein_abundance__filter by contrast; Ctrl vs Wc # condition_variable_ group.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
Wc_vs_Ctrl$avg_log2_Intens <- rowMeans(Wc_vs_Ctrl[, 4:11], na.rm = TRUE)
Wc_vs_Ctrl$avg_Intens <- 2^(Wc_vs_Ctrl$avg_log2_Intens)
Wc_vs_Ctrl$log10_Intens <- log10(Wc_vs_Ctrl$avg_Intens)
Wc_vs_Ctrl <- Wc_vs_Ctrl %>% arrange( desc(log10_Intens))
Wc_vs_Ctrl$rank <- 1:3435
Wc_vs_Ctrl_astrocyte    <- Wc_vs_Ctrl %>%  left_join(astrocyte_marker, by=c("gene_symbols_or_id"="Unique_Gene_Marker"))
Wc_vs_Ctrl_astrogliosis <- Wc_vs_Ctrl %>%  left_join(astrogliosis_marker, by=c("gene_symbols_or_id"="GeneSymbol_Marker_unique"))

#A1 <- 
Wc_vs_Ctrl_astrocyte  %>%  
  ggplot( aes(rank, log10_Intens, color=Function, fill=Function) ) +
  geom_point(alpha=0.5, size=3 ) + 
  scale_color_manual(values=c("Metabolic"="blue4", 
                              "Struct & Memb"="red4", 
                              "Transcrip Fact & intracell"="black", 
                              "NA"="grey40") ) + 
  geom_text_repel(data=Wc_vs_Ctrl_astrocyte %>% filter(general_astrocyte_markers==TRUE), 
                  aes(label=gene_symbols_or_id),
                  fontface='bold', segment.color='grey30',
                  size = 4, max.overlaps = (getOption("ggrepel.max.overlaps", default=30)), min.segment.length=.2, inherit.aes=T, show.legend=F, na.rm=F) + 
  geom_point(data=Wc_vs_Ctrl_astrocyte %>% filter(Function=="Metabolic"), 
             aes(x=rank, y=log10_Intens), color="blue4", size=3) +
  geom_point(data=Wc_vs_Ctrl_astrocyte %>% filter(Function=="Struct & Memb"), 
             aes(x=rank, y=log10_Intens), color="red4", size=3) +
  geom_point(data=Wc_vs_Ctrl_astrocyte %>% filter(Function=="Transcrip Fact & intracell"), 
             aes(x=rank, y=log10_Intens), color="black", size=3) +
  theme_bw( ) + labs(x= "Proteins rank", y= "log10 (LFQ mean)", title = "") +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "inside", legend.justification = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        strip.text = element_text(size=14) ) +
  geom_text(aes(x=1500, y=9, label="Wc_vs_Ctrl"), size=6, color="black") 

#A2 <- 
Wc_vs_Ctrl_astrogliosis %>%
  ggplot( aes(rank, log10_Intens, color=Localization, fill=Localization, label=gene_symbols_or_id) ) +
  geom_point(alpha=0.5, size=3 ) + 
  scale_color_manual(values=c("Intracell"="black", 
                              "Struct & Memb"="red4", 
                              "Secreted"="orange4", 
                              "NA"="grey40") ) + 
  geom_text_repel(data=Wc_vs_Ctrl_astrogliosis %>% filter(markers_reactive_astrogliosis==TRUE), 
                  aes(label=gene_symbols_or_id),
                  fontface='bold', segment.color='grey30',
                  size = 4, max.overlaps = (getOption("ggrepel.max.overlaps", default=30)), min.segment.length=.2, inherit.aes=T, show.legend=F, na.rm=F) +
  geom_point(data=Wc_vs_Ctrl_astrogliosis %>% filter(Localization=="Intracell"), 
             aes(x=rank, y=log10_Intens), color="black", size=3) +
  geom_point(data=Wc_vs_Ctrl_astrogliosis %>% filter(Localization=="Struct & Memb"), 
             aes(x=rank, y=log10_Intens), color="red4", size=3) +
  geom_point(data=Wc_vs_Ctrl_astrogliosis %>% filter(Localization=="Secreted"), 
             aes(x=rank, y=log10_Intens), color="orange4", size=3) +
  theme_bw( ) + labs(x= "Proteins rank", y= "log10 (LFQ mean)", title = " ") +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "inside", legend.justification = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        strip.text = element_text(size=14) ) +
  geom_text(aes(x=1500, y=9, label="Wc_vs_Ctrl"), size=6, color="black")

##### Wa_vs_Ctrl #######
Wa_vs_Ctrl <- read_delim("protein_abundance__filter by contrast; Ctrl vs Wa # condition_variable_ group.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
Wa_vs_Ctrl$avg_log2_Intens <- rowMeans(Wa_vs_Ctrl[, 4:10], na.rm = TRUE)
Wa_vs_Ctrl$avg_Intens <- 2^(Wa_vs_Ctrl$avg_log2_Intens)
Wa_vs_Ctrl$log10_Intens <- log10(Wa_vs_Ctrl$avg_Intens)
Wa_vs_Ctrl <- Wa_vs_Ctrl %>% arrange( desc(log10_Intens))
Wa_vs_Ctrl$rank <- 1:3127
Wa_vs_Ctrl_astrocyte    <- Wa_vs_Ctrl %>%  left_join(astrocyte_marker, by=c("gene_symbols_or_id"="Unique_Gene_Marker"))
Wa_vs_Ctrl_astrogliosis <- Wa_vs_Ctrl %>%  left_join(astrogliosis_marker, by=c("gene_symbols_or_id"="GeneSymbol_Marker_unique"))

#B1 <- 
Wa_vs_Ctrl_astrocyte  %>%  
  ggplot( aes(rank, log10_Intens, color=Function, fill=Function) ) +
  geom_point(alpha=0.5, size=3 ) + 
  scale_color_manual(values=c("Metabolic"="blue4", 
                              "Struct & Memb"="red4", 
                              "Transcrip Fact & intracell"="black", 
                              "NA"="grey40") ) + 
  geom_text_repel(data=Wa_vs_Ctrl_astrocyte %>% filter(general_astrocyte_markers==TRUE), 
                  aes(label=gene_symbols_or_id),
                  fontface='bold', segment.color='grey30',
                  size = 4, max.overlaps = (getOption("ggrepel.max.overlaps", default=30)), min.segment.length=.2, inherit.aes=T, show.legend=F, na.rm=F) + 
  geom_point(data=Wa_vs_Ctrl_astrocyte %>% filter(Function=="Metabolic"), 
             aes(x=rank, y=log10_Intens), color="blue4", size=3) +
  geom_point(data=Wa_vs_Ctrl_astrocyte %>% filter(Function=="Struct & Memb"), 
             aes(x=rank, y=log10_Intens), color="red4", size=3) +
  geom_point(data=Wa_vs_Ctrl_astrocyte %>% filter(Function=="Transcrip Fact & intracell"), 
             aes(x=rank, y=log10_Intens), color="black", size=3) +
  theme_bw( ) + labs(x= "Proteins rank", y= "log10 (LFQ mean)", title = "") +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "inside", legend.justification = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        strip.text = element_text(size=14) ) +
  geom_text(aes(x=1500, y=9, label="Wa_vs_Ctrl"), size=6, color="black")

#B2 <- 
Wa_vs_Ctrl_astrogliosis %>%
  ggplot( aes(rank, log10_Intens, color=Localization, fill=Localization, label=gene_symbols_or_id) ) +
  geom_point(alpha=0.5, size=3 ) + 
  scale_color_manual(values=c("Intracell"="black", 
                              "Struct & Memb"="red4", 
                              "Secreted"="orange4", 
                              "NA"="grey40") ) + 
  geom_text_repel(data=Wa_vs_Ctrl_astrogliosis %>% filter(markers_reactive_astrogliosis==TRUE), 
                  aes(label=gene_symbols_or_id),
                  fontface='bold', segment.color='grey30',
                  size = 4, max.overlaps = (getOption("ggrepel.max.overlaps", default=30)), min.segment.length=.2, inherit.aes=T, show.legend=F, na.rm=F) +
  geom_point(data=Wa_vs_Ctrl_astrogliosis %>% filter(Localization=="Intracell"), 
             aes(x=rank, y=log10_Intens), color="black", size=3) +
  geom_point(data=Wa_vs_Ctrl_astrogliosis %>% filter(Localization=="Struct & Memb"), 
             aes(x=rank, y=log10_Intens), color="red4", size=3) +
  geom_point(data=Wa_vs_Ctrl_astrogliosis %>% filter(Localization=="Secreted"), 
             aes(x=rank, y=log10_Intens), color="orange4", size=3) +
  theme_bw( ) + labs(x= "Proteins rank", y= "log10 (LFQ mean)", title = " ") +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "inside", legend.justification = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        strip.text = element_text(size=14) ) +
  geom_text(aes(x=1500, y=9, label="Wa_vs_Ctrl"), size=6, color="black")


##### Gc_vs_Ctrl #######
Gc_vs_Ctrl <- read_delim("protein_abundance__filter by contrast; Ctrl vs Gc # condition_variable_ group.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
Gc_vs_Ctrl$avg_log2_Intens <- rowMeans(Gc_vs_Ctrl[, 4:11], na.rm = TRUE)
Gc_vs_Ctrl$avg_Intens <- 2^(Gc_vs_Ctrl$avg_log2_Intens)
Gc_vs_Ctrl$log10_Intens <- log10(Gc_vs_Ctrl$avg_Intens)
Gc_vs_Ctrl <- Gc_vs_Ctrl %>% arrange( desc(log10_Intens))
Gc_vs_Ctrl$rank <- 1:3275
Gc_vs_Ctrl_astrocyte    <- Gc_vs_Ctrl %>%  left_join(astrocyte_marker, by=c("gene_symbols_or_id"="Unique_Gene_Marker"))
Gc_vs_Ctrl_astrogliosis <- Gc_vs_Ctrl %>%  left_join(astrogliosis_marker, by=c("gene_symbols_or_id"="GeneSymbol_Marker_unique"))

#C1 <- 
Gc_vs_Ctrl_astrocyte  %>%  
  ggplot( aes(rank, log10_Intens, color=Function, fill=Function) ) +
  geom_point(alpha=0.5, size=3 ) + 
  scale_color_manual(values=c("Metabolic"="blue4", 
                              "Struct & Memb"="red4", 
                              "Transcrip Fact & intracell"="black", 
                              "NA"="grey40") ) + 
  geom_text_repel(data=Gc_vs_Ctrl_astrocyte %>% filter(general_astrocyte_markers==TRUE), 
                  aes(label=gene_symbols_or_id),
                  fontface='bold', segment.color='grey30',
                  size = 4, max.overlaps = (getOption("ggrepel.max.overlaps", default=30)), min.segment.length=.2, inherit.aes=T, show.legend=F, na.rm=F) + 
  geom_point(data=Gc_vs_Ctrl_astrocyte %>% filter(Function=="Metabolic"), 
             aes(x=rank, y=log10_Intens), color="blue4", size=3) +
  geom_point(data=Gc_vs_Ctrl_astrocyte %>% filter(Function=="Struct & Memb"), 
             aes(x=rank, y=log10_Intens), color="red4", size=3) +
  geom_point(data=Gc_vs_Ctrl_astrocyte %>% filter(Function=="Transcrip Fact & intracell"), 
             aes(x=rank, y=log10_Intens), color="black", size=3) +
  theme_bw( ) + labs(x= "Proteins rank", y= "log10 (LFQ mean)", title = "") +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "inside", legend.justification = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        strip.text = element_text(size=14) ) +
  geom_text(aes(x=1500, y=9, label="Gc_vs_Ctrl"), size=6, color="black")

#C2 <- 
Gc_vs_Ctrl_astrogliosis %>%
  ggplot( aes(rank, log10_Intens, color=Localization, fill=Localization, label=gene_symbols_or_id) ) +
  geom_point(alpha=0.5, size=3 ) + 
  scale_color_manual(values=c("Intracell"="black", 
                              "Struct & Memb"="red4", 
                              "Secreted"="orange4", 
                              "NA"="grey40") ) + 
  geom_text_repel(data=Gc_vs_Ctrl_astrogliosis %>% filter(markers_reactive_astrogliosis==TRUE), 
                  aes(label=gene_symbols_or_id),
                  fontface='bold', segment.color='grey30',
                  size = 4, max.overlaps = (getOption("ggrepel.max.overlaps", default=30)), min.segment.length=.2, inherit.aes=T, show.legend=F, na.rm=F) +
  geom_point(data=Gc_vs_Ctrl_astrogliosis %>% filter(Localization=="Intracell"), 
             aes(x=rank, y=log10_Intens), color="black", size=3) +
  geom_point(data=Gc_vs_Ctrl_astrogliosis %>% filter(Localization=="Struct & Memb"), 
             aes(x=rank, y=log10_Intens), color="red4", size=3) +
  geom_point(data=Gc_vs_Ctrl_astrogliosis %>% filter(Localization=="Secreted"), 
             aes(x=rank, y=log10_Intens), color="orange4", size=3) +
  theme_bw( ) + labs(x= "Proteins rank", y= "log10 (LFQ mean)", title = " ") +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "inside", legend.justification = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        strip.text = element_text(size=14) ) +
  geom_text(aes(x=1500, y=9, label="Gc_vs_Ctrl"), size=6, color="black")

##### Ga_vs_Ctrl #######
Ga_vs_Ctrl <- read_delim("protein_abundance__filter by contrast; Ctrl vs Ga # condition_variable_ group.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
Ga_vs_Ctrl$avg_log2_Intens <- rowMeans(Ga_vs_Ctrl[, 4:10], na.rm = TRUE)
Ga_vs_Ctrl$avg_Intens <- 2^(Ga_vs_Ctrl$avg_log2_Intens)
Ga_vs_Ctrl$log10_Intens <- log10(Ga_vs_Ctrl$avg_Intens)
Ga_vs_Ctrl <- Ga_vs_Ctrl %>% arrange( desc(log10_Intens))
Ga_vs_Ctrl$rank <- 1:2892
Ga_vs_Ctrl_astrocyte    <- Ga_vs_Ctrl %>%  left_join(astrocyte_marker, by=c("gene_symbols_or_id"="Unique_Gene_Marker"))
Ga_vs_Ctrl_astrogliosis <- Ga_vs_Ctrl %>%  left_join(astrogliosis_marker, by=c("gene_symbols_or_id"="GeneSymbol_Marker_unique"))

#D1 <- 
Ga_vs_Ctrl_astrocyte  %>%  
  ggplot( aes(rank, log10_Intens, color=Function, fill=Function) ) +
  geom_point(alpha=0.5, size=3 ) + 
  scale_color_manual(values=c("Metabolic"="blue4", 
                              "Struct & Memb"="red4", 
                              "Transcrip Fact & intracell"="black", 
                              "NA"="grey40") ) + 
  geom_text_repel(data=Ga_vs_Ctrl_astrocyte %>% filter(general_astrocyte_markers==TRUE), 
                  aes(label=gene_symbols_or_id),
                  fontface='bold', segment.color='grey30',
                  size = 4, max.overlaps = (getOption("ggrepel.max.overlaps", default=30)), min.segment.length=.2, inherit.aes=T, show.legend=F, na.rm=F) + 
  geom_point(data=Ga_vs_Ctrl_astrocyte %>% filter(Function=="Metabolic"), 
             aes(x=rank, y=log10_Intens), color="blue4", size=3) +
  geom_point(data=Ga_vs_Ctrl_astrocyte %>% filter(Function=="Struct & Memb"), 
             aes(x=rank, y=log10_Intens), color="red4", size=3) +
  geom_point(data=Ga_vs_Ctrl_astrocyte %>% filter(Function=="Transcrip Fact & intracell"), 
             aes(x=rank, y=log10_Intens), color="black", size=3) +
  theme_bw( ) + labs(x= "Proteins rank", y= "log10 (LFQ mean)", title = "") +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "inside", 
        legend.justification = "bottom",
        legend.justification.inside = c(.01, .01),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        strip.text = element_text(size=14) ) +
  geom_text(aes(x=1500, y=9, label="Ga_vs_Ctrl"), size=6, color="black") +
  scale_x_continuous(limits=c(0, 3020))

#D2 <- 
Ga_vs_Ctrl_astrogliosis %>%
  ggplot( aes(rank, log10_Intens, color=Localization, fill=Localization, label=gene_symbols_or_id) ) +
  geom_point(alpha=0.5, size=3 ) + 
  scale_color_manual(values=c("Intracell"="black", 
                              "Struct & Memb"="red4", 
                              "Secreted"="orange4", 
                              "NA"="grey40") ) + 
  geom_text_repel(data=Ga_vs_Ctrl_astrogliosis %>% filter(markers_reactive_astrogliosis==TRUE), 
                  aes(label=gene_symbols_or_id),
                  fontface='bold', segment.color='grey30',
                  size = 4, max.overlaps = (getOption("ggrepel.max.overlaps", default=30)), min.segment.length=.2, inherit.aes=T, show.legend=F, na.rm=F) +
  geom_point(data=Ga_vs_Ctrl_astrogliosis %>% filter(Localization=="Intracell"), 
             aes(x=rank, y=log10_Intens), color="black", size=3) +
  geom_point(data=Ga_vs_Ctrl_astrogliosis %>% filter(Localization=="Struct & Memb"), 
             aes(x=rank, y=log10_Intens), color="red4", size=3) +
  geom_point(data=Ga_vs_Ctrl_astrogliosis %>% filter(Localization=="Secreted"), 
             aes(x=rank, y=log10_Intens), color="orange4", size=3) +
  theme_bw( ) + labs(x= "Proteins rank", y= "log10 (LFQ mean)", title = " ") +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "inside", 
        legend.justification = "bottom",
        legend.justification.inside = c(.01, .01),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        strip.text = element_text(size=14) ) +
  geom_text(aes(x=1500, y=9, label="Ga_vs_Ctrl"), size=6, color="black") +
  scale_x_continuous(limits=c(0, 3020))



##### Wc_vs_Wa #######
Wc_vs_Wa <- read_delim("protein_abundance__filter by contrast; Wa vs Wc # condition_variable_ group.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
Wc_vs_Wa$avg_log2_Intens <- rowMeans(Wc_vs_Wa[, 4:10], na.rm = TRUE)
Wc_vs_Wa$avg_Intens <- 2^(Wc_vs_Wa$avg_log2_Intens)
Wc_vs_Wa$log10_Intens <- log10(Wc_vs_Wa$avg_Intens)
Wc_vs_Wa <- Wc_vs_Wa %>% arrange( desc(log10_Intens))
Wc_vs_Wa$rank <- 1:4963
Wc_vs_Wa_astrocyte    <- Wc_vs_Wa %>%  left_join(astrocyte_marker, by=c("gene_symbols_or_id"="Unique_Gene_Marker"))
Wc_vs_Wa_astrogliosis <- Wc_vs_Wa %>%  left_join(astrogliosis_marker, by=c("gene_symbols_or_id"="GeneSymbol_Marker_unique"))

Wc_vs_Wa_astrocyte  %>%  
  ggplot( aes(rank, log10_Intens, color=Function, fill=Function) ) +
  geom_point(alpha=0.5, size=3 ) + 
  scale_color_manual(values=c("Metabolic"="blue4", 
                              "Struct & Memb"="red4", 
                              "Transcrip Fact & intracell"="black", 
                              "NA"="grey40") ) + 
  geom_text_repel(data=Wc_vs_Wa_astrocyte %>% filter(general_astrocyte_markers==TRUE), 
                  aes(label=gene_symbols_or_id),
                  fontface='bold', segment.color='grey30',
                  size = 4, max.overlaps = (getOption("ggrepel.max.overlaps", default=30)), min.segment.length=.2, inherit.aes=T, show.legend=F, na.rm=F) + 
  geom_point(data=Wc_vs_Wa_astrocyte %>% filter(Function=="Metabolic"), 
             aes(x=rank, y=log10_Intens), color="blue4", size=3) +
  geom_point(data=Wc_vs_Wa_astrocyte %>% filter(Function=="Struct & Memb"), 
             aes(x=rank, y=log10_Intens), color="red4", size=3) +
  geom_point(data=Wc_vs_Wa_astrocyte %>% filter(Function=="Transcrip Fact & intracell"), 
             aes(x=rank, y=log10_Intens), color="black", size=3) +
  theme_bw( ) + labs(x= "Proteins rank", y= "log10 (LFQ mean)", title = "") +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "inside", 
        legend.position.inside = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        strip.text = element_text(size=14) ) +
  geom_text(aes(x=1000, y=9.5, label="Wc_vs_Wa"), size=6, color="black")

Wc_vs_Wa_astrogliosis %>%
  ggplot( aes(rank, log10_Intens, color=Localization, fill=Localization, label=gene_symbols_or_id) ) +
  geom_point(alpha=0.5, size=3 ) + 
  scale_color_manual(values=c("Intracell"="black", 
                              "Struct & Memb"="red4", 
                              "Secreted"="orange4", 
                              "NA"="grey40") ) + 
  geom_text_repel(data=Wc_vs_Wa_astrogliosis %>% filter(markers_reactive_astrogliosis==TRUE), 
                  aes(label=gene_symbols_or_id),
                  fontface='bold', segment.color='grey30',
                  size = 4, max.overlaps = (getOption("ggrepel.max.overlaps", default=30)), min.segment.length=.2, inherit.aes=T, show.legend=F, na.rm=F) +
  geom_point(data=Wc_vs_Wa_astrogliosis %>% filter(Localization=="Intracell"), 
             aes(x=rank, y=log10_Intens), color="black", size=3) +
  geom_point(data=Wc_vs_Wa_astrogliosis %>% filter(Localization=="Struct & Memb"), 
             aes(x=rank, y=log10_Intens), color="red4", size=3) +
  geom_point(data=Wc_vs_Wa_astrogliosis %>% filter(Localization=="Secreted"), 
             aes(x=rank, y=log10_Intens), color="orange4", size=3) +
  theme_bw( ) + labs(x= "Proteins rank", y= "log10 (LFQ mean)", title = " ") +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "inside", 
        legend.position.inside = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        strip.text = element_text(size=14) ) +
  geom_text(aes(x=1000, y=9.5, label="Wc_vs_Wa"), size=6, color="black")

##### Wc_vs_Gc #######
Wc_vs_Gc <- read_delim("protein_abundance__filter by contrast; Gc vs Wc # condition_variable_ group.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
Wc_vs_Gc$avg_log2_Intens <- rowMeans(Wc_vs_Gc[, 4:11], na.rm = TRUE)
Wc_vs_Gc$avg_Intens <- 2^(Wc_vs_Gc$avg_log2_Intens)
Wc_vs_Gc$log10_Intens <- log10(Wc_vs_Gc$avg_Intens)
Wc_vs_Gc <- Wc_vs_Gc %>% arrange( desc(log10_Intens))
Wc_vs_Gc$rank <- 1:5209
Wc_vs_Gc_astrocyte    <- Wc_vs_Gc %>%  left_join(astrocyte_marker, by=c("gene_symbols_or_id"="Unique_Gene_Marker"))
Wc_vs_Gc_astrogliosis <- Wc_vs_Gc %>%  left_join(astrogliosis_marker, by=c("gene_symbols_or_id"="GeneSymbol_Marker_unique"))

Wc_vs_Gc_astrocyte  %>%  
  ggplot( aes(rank, log10_Intens, color=Function, fill=Function) ) +
  geom_point(alpha=0.5, size=3 ) + 
  scale_color_manual(values=c("Metabolic"="blue4", 
                              "Struct & Memb"="red4", 
                              "Transcrip Fact & intracell"="black", 
                              "NA"="grey40") ) + 
  geom_text_repel(data=Wc_vs_Gc_astrocyte %>% filter(general_astrocyte_markers==TRUE), 
                  aes(label=gene_symbols_or_id),
                  fontface='bold', segment.color='grey30',
                  size = 4, max.overlaps = (getOption("ggrepel.max.overlaps", default=30)), min.segment.length=.2, inherit.aes=T, show.legend=F, na.rm=F) + 
  geom_point(data=Wc_vs_Gc_astrocyte %>% filter(Function=="Metabolic"), 
             aes(x=rank, y=log10_Intens), color="blue4", size=3) +
  geom_point(data=Wc_vs_Gc_astrocyte %>% filter(Function=="Struct & Memb"), 
             aes(x=rank, y=log10_Intens), color="red4", size=3) +
  geom_point(data=Wc_vs_Gc_astrocyte %>% filter(Function=="Transcrip Fact & intracell"), 
             aes(x=rank, y=log10_Intens), color="black", size=3) +
  theme_bw( ) + labs(x= "Proteins rank", y= "log10 (LFQ mean)", title = "") +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "inside", 
        legend.position.inside = c(.99, .99),
        legend.justification = c("right", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        strip.text = element_text(size=14) ) +
  geom_text(aes(x=1000, y=9, label="Wc_vs_Gc"), size=6, color="black")

Wc_vs_Gc_astrogliosis %>%
  ggplot( aes(rank, log10_Intens, color=Localization, fill=Localization, label=gene_symbols_or_id) ) +
  geom_point(alpha=0.5, size=3 ) + 
  scale_color_manual(values=c("Intracell"="black", 
                              "Struct & Memb"="red4", 
                              "Secreted"="orange4", 
                              "NA"="grey40") ) + 
  geom_text_repel(data=Wc_vs_Gc_astrogliosis %>% filter(markers_reactive_astrogliosis==TRUE), 
                  aes(label=gene_symbols_or_id),
                  fontface='bold', segment.color='grey30',
                  size = 4, max.overlaps = (getOption("ggrepel.max.overlaps", default=30)), min.segment.length=.2, inherit.aes=T, show.legend=F, na.rm=F) +
  geom_point(data=Wc_vs_Gc_astrogliosis %>% filter(Localization=="Intracell"), 
             aes(x=rank, y=log10_Intens), color="black", size=3) +
  geom_point(data=Wc_vs_Gc_astrogliosis %>% filter(Localization=="Struct & Memb"), 
             aes(x=rank, y=log10_Intens), color="red4", size=3) +
  geom_point(data=Wc_vs_Gc_astrogliosis %>% filter(Localization=="Secreted"), 
             aes(x=rank, y=log10_Intens), color="orange4", size=3) +
  theme_bw( ) + labs(x= "Proteins rank", y= "log10 (LFQ mean)", title = " ") +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "inside", 
        legend.position.inside = c(.99, .99),
        legend.justification = c("right", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        strip.text = element_text(size=14) ) +
  geom_text(aes(x=1000, y=9, label="Wc_vs_Gc"), size=6, color="black")

##### Wa_vs_Ga #######
Wa_vs_Ga <- read_delim("protein_abundance__filter by contrast; Ga vs Wa # condition_variable_ group.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
Wa_vs_Ga$avg_log2_Intens <- rowMeans(Wa_vs_Ga[, 4:9], na.rm = TRUE)
Wa_vs_Ga$avg_Intens <- 2^(Wa_vs_Ga$avg_log2_Intens)
Wa_vs_Ga$log10_Intens <- log10(Wa_vs_Ga$avg_Intens)
Wa_vs_Ga <- Wa_vs_Ga %>% arrange( desc(log10_Intens))
Wa_vs_Ga$rank <- 1:4369
Wa_vs_Ga_astrocyte    <- Wa_vs_Ga %>%  left_join(astrocyte_marker, by=c("gene_symbols_or_id"="Unique_Gene_Marker"))
Wa_vs_Ga_astrogliosis <- Wa_vs_Ga %>%  left_join(astrogliosis_marker, by=c("gene_symbols_or_id"="GeneSymbol_Marker_unique"))

Wa_vs_Ga_astrocyte  %>%  
  ggplot( aes(rank, log10_Intens, color=Function, fill=Function) ) +
  geom_point(alpha=0.5, size=3 ) + 
  scale_color_manual(values=c("Metabolic"="blue4", 
                              "Struct & Memb"="red4", 
                              "Transcrip Fact & intracell"="black", 
                              "NA"="grey40") ) + 
  geom_text_repel(data=Wa_vs_Ga_astrocyte %>% filter(general_astrocyte_markers==TRUE), 
                  aes(label=gene_symbols_or_id),
                  fontface='bold', segment.color='grey30',
                  size = 4, max.overlaps = (getOption("ggrepel.max.overlaps", default=30)), min.segment.length=.2, inherit.aes=T, show.legend=F, na.rm=F) + 
  geom_point(data=Wa_vs_Ga_astrocyte %>% filter(Function=="Metabolic"), 
             aes(x=rank, y=log10_Intens), color="blue4", size=3) +
  geom_point(data=Wa_vs_Ga_astrocyte %>% filter(Function=="Struct & Memb"), 
             aes(x=rank, y=log10_Intens), color="red4", size=3) +
  geom_point(data=Wa_vs_Ga_astrocyte %>% filter(Function=="Transcrip Fact & intracell"), 
             aes(x=rank, y=log10_Intens), color="black", size=3) +
  theme_bw( ) + labs(x= "Proteins rank", y= "log10 (LFQ mean)", title = "") +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "inside", 
        legend.position.inside = c(.99, .99),
        legend.justification = c("right", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        strip.text = element_text(size=14) ) +
  geom_text(aes(x=1000, y=9.5, label="Wa_vs_Ga"), size=6, color="black")

Wa_vs_Ga_astrogliosis %>%
  ggplot( aes(rank, log10_Intens, color=Localization, fill=Localization, label=gene_symbols_or_id) ) +
  geom_point(alpha=0.5, size=3 ) + 
  scale_color_manual(values=c("Intracell"="black", 
                              "Struct & Memb"="red4", 
                              "Secreted"="orange4", 
                              "NA"="grey40") ) + 
  geom_text_repel(data=Wa_vs_Ga_astrogliosis %>% filter(markers_reactive_astrogliosis==TRUE), 
                  aes(label=gene_symbols_or_id),
                  fontface='bold', segment.color='grey30',
                  size = 4, max.overlaps = (getOption("ggrepel.max.overlaps", default=30)), min.segment.length=.2, inherit.aes=T, show.legend=F, na.rm=F) +
  geom_point(data=Wa_vs_Ga_astrogliosis %>% filter(Localization=="Intracell"), 
             aes(x=rank, y=log10_Intens), color="black", size=3) +
  geom_point(data=Wa_vs_Ga_astrogliosis %>% filter(Localization=="Struct & Memb"), 
             aes(x=rank, y=log10_Intens), color="red4", size=3) +
  geom_point(data=Wa_vs_Ga_astrogliosis %>% filter(Localization=="Secreted"), 
             aes(x=rank, y=log10_Intens), color="orange4", size=3) +
  theme_bw( ) + labs(x= "Proteins rank", y= "log10 (LFQ mean)", title = " ") +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "inside", 
        legend.position.inside = c(.99, .99),
        legend.justification = c("right", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        strip.text = element_text(size=14) ) +
  geom_text(aes(x=1000, y=9.5, label="Wa_vs_Ga"), size=6, color="black")

##### Gc_vs_Ga #######
Gc_vs_Ga <- read_delim("protein_abundance__filter by contrast; Ga vs Gc # condition_variable_ group.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
Gc_vs_Ga$avg_log2_Intens <- rowMeans(Gc_vs_Ga[, 4:10], na.rm = TRUE)
Gc_vs_Ga$avg_Intens <- 2^(Gc_vs_Ga$avg_log2_Intens)
Gc_vs_Ga$log10_Intens <- log10(Gc_vs_Ga$avg_Intens)
Gc_vs_Ga <- Gc_vs_Ga %>% arrange( desc(log10_Intens))
Gc_vs_Ga$rank <- 1:4563
Gc_vs_Ga_astrocyte    <- Gc_vs_Ga %>%  left_join(astrocyte_marker, by=c("gene_symbols_or_id"="Unique_Gene_Marker"))
Gc_vs_Ga_astrogliosis <- Gc_vs_Ga %>%  left_join(astrogliosis_marker, by=c("gene_symbols_or_id"="GeneSymbol_Marker_unique"))

Gc_vs_Ga_astrocyte  %>%  
  ggplot( aes(rank, log10_Intens, color=Function, fill=Function) ) +
  geom_point(alpha=0.5, size=3 ) + 
  scale_color_manual(values=c("Metabolic"="blue4", 
                              "Struct & Memb"="red4", 
                              "Transcrip Fact & intracell"="black", 
                              "NA"="grey40") ) + 
  geom_text_repel(data=Gc_vs_Ga_astrocyte %>% filter(general_astrocyte_markers==TRUE), 
                  aes(label=gene_symbols_or_id),
                  fontface='bold', segment.color='grey30',
                  size = 4, max.overlaps = (getOption("ggrepel.max.overlaps", default=30)), min.segment.length=.2, inherit.aes=T, show.legend=F, na.rm=F) + 
  geom_point(data=Gc_vs_Ga_astrocyte %>% filter(Function=="Metabolic"), 
             aes(x=rank, y=log10_Intens), color="blue4", size=3) +
  geom_point(data=Gc_vs_Ga_astrocyte %>% filter(Function=="Struct & Memb"), 
             aes(x=rank, y=log10_Intens), color="red4", size=3) +
  geom_point(data=Gc_vs_Ga_astrocyte %>% filter(Function=="Transcrip Fact & intracell"), 
             aes(x=rank, y=log10_Intens), color="black", size=3) +
  theme_bw( ) + labs(x= "Proteins rank", y= "log10 (LFQ mean)", title = "") +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "inside", 
        legend.position.inside = c(.99, .99),
        legend.justification = c("right", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        strip.text = element_text(size=14) ) +
  geom_text(aes(x=1000, y=9.5, label="Gc_vs_Ga"), size=6, color="black")

Gc_vs_Ga_astrogliosis %>%
  ggplot( aes(rank, log10_Intens, color=Localization, fill=Localization, label=gene_symbols_or_id) ) +
  geom_point(alpha=0.5, size=3 ) + 
  scale_color_manual(values=c("Intracell"="black", 
                              "Struct & Memb"="red4", 
                              "Secreted"="orange4", 
                              "NA"="grey40") ) + 
  geom_text_repel(data=Gc_vs_Ga_astrogliosis %>% filter(markers_reactive_astrogliosis==TRUE), 
                  aes(label=gene_symbols_or_id),
                  fontface='bold', segment.color='grey30',
                  size = 4, max.overlaps = (getOption("ggrepel.max.overlaps", default=30)), min.segment.length=.2, inherit.aes=T, show.legend=F, na.rm=F) +
  geom_point(data=Gc_vs_Ga_astrogliosis %>% filter(Localization=="Intracell"), 
             aes(x=rank, y=log10_Intens), color="black", size=3) +
  geom_point(data=Gc_vs_Ga_astrogliosis %>% filter(Localization=="Struct & Memb"), 
             aes(x=rank, y=log10_Intens), color="red4", size=3) +
  geom_point(data=Gc_vs_Ga_astrogliosis %>% filter(Localization=="Secreted"), 
             aes(x=rank, y=log10_Intens), color="orange4", size=3) +
  theme_bw( ) + labs(x= "Proteins rank", y= "log10 (LFQ mean)", title = " ") +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "inside", 
        legend.position.inside = c(.99, .99),
        legend.justification = c("right", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        strip.text = element_text(size=14) ) +
  geom_text(aes(x=1000, y=9.5, label="Gc_vs_Ga"), size=6, color="black")
