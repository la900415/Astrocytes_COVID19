library(tidyverse)
library(readxl)
library(ggrepel)

# import table
setwd("Q:/data_analysis/COVID_Prot_USP/celulas/astrocytes_COVID19/study_DIANN/22_reanalysis_same_expSpecLib_61samp_without_isoforms_FINAL/msdap_results/2024-10-09_23-34-50_combinedCTRL_msqrob_FDR0.01_FINAL-BEST")
global_2724_astrocyte_long <- read_excel("Fig2_rank_plot/global_2724_astrocyte_long.xlsx")

#Avg_all <- 
global_2724_astrocyte_long  %>%  
  dplyr::filter(sample == "Avg_all") %>%
  ggplot( aes(rank, log10_Intens, color=Function, fill=Function) ) +
  geom_point(alpha=0.5, size=3 ) + 
  scale_color_manual(values=c("Metabolic"="blue4", 
                              "Struct & Memb"="red4", 
                              "Transcrip Fact & intracell"="black", 
                              "non-markers"="grey40") ) + 
  geom_text_repel(data=global_2724_astrocyte_long %>% filter(general_astrocyte_markers==TRUE) %>% filter(sample=="Avg_all"), 
                  aes(label=gene_symbols_or_id),
                  fontface='bold', segment.color='grey30',
                  size = 4, max.overlaps = (getOption("ggrepel.max.overlaps", default=30)), min.segment.length=.2, inherit.aes=T, show.legend=F, na.rm=F) + 
  geom_point(data=global_2724_astrocyte_long %>% filter(Function=="Metabolic") %>% filter(sample=="Avg_all"), 
             aes(x=rank, y=log10_Intens), color="blue4", size=3) +
  geom_point(data=global_2724_astrocyte_long %>% filter(Function=="Struct & Memb") %>% filter(sample=="Avg_all"), 
             aes(x=rank, y=log10_Intens), color="red4", size=3) +
  geom_point(data=global_2724_astrocyte_long %>% filter(Function=="Transcrip Fact & intracell") %>% filter(sample=="Avg_all"), 
             aes(x=rank, y=log10_Intens), color="black", size=3) +
  scale_y_continuous(limits=c(4, 10)) + scale_x_continuous(limits=c(0, 3000)) +
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
        strip.text = element_text(size=14) ) #
#geom_text(aes(x=1500, y=10, label="Astrocyte markers"), size=6, color="black")