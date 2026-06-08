setwd("C:/Users/Laura/Downloads/astrocitos_COVID19/2nd_reanalysis/msdap_results/2024-10-09_23-34-50_combinedCTRL_msqrob_FDR0.01")
load("dataset.RData")

################ Fig 3B-3C, S3B-S3C. Identif & quantif #########################################
library(ggsci)
library(RColorBrewer)
library(tidyverse)
cbbPalette <- c("#7E6148B2", "#E64B35B2", "#F39B7FB2", "#3C5488B2", "#8491B4B1") #palette for Ctrl, Ga, Gc, Wa, Wc groups
dataset$samples$shortname <-  dataset$samples$shortname %>% str_replace_all ( c ("Ca1"="Ctrl1", "Ca3"="Ctrl3", "Cc1"="Ctrl4", "Cc2"="Ctrl2" ) )

# Identified proteins
#A <- 
dataset$samples %>% dplyr::filter(exclude=="FALSE") %>% 
  ggplot( aes(x=shortname, y=detected_proteins, fill=group) ) +
  geom_bar(position = position_dodge(width = 0.9), stat="identity", width=.8, colour="black") +
  geom_hline(yintercept=5320, color = "#E64B35B2", linetype = "dashed", size=0.8) +
  geom_hline(yintercept=3980, color = "#7E6148B2", linetype = "dashed", size=0.8) +
  geom_text(aes(label=detected_proteins),angle=90,vjust=.25,hjust=1.2,size=2.5,color="black") +
  theme_bw() +
  xlab(label = "Samples") + ylab(label = "# identif. proteins") + 
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = "top",
        legend.text = element_text(size = 11),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 11, angle = 90, vjust=0.5),
        axis.text.y = element_text(size = 12) ) +
  scale_fill_manual(values = c("Ctrl"="#7E6148B2", "Ga"="#E64B35B2","Gc"="#F39B7FB2", "Wa"="#3C5488B2", "Wc"="#8491B4B1") ) 
#scale_fill_jama() + # library(ggsci)
#scale_fill_brewer(type="seq",palette="Blues",direction=1,guide="coloursteps",aesthetics = "fill")

# Identified peptides
#B <- 
dataset$samples %>% filter(exclude=="FALSE") %>% 
  ggplot( aes(x=shortname, y=detected_peptides, fill=group) ) +
  geom_bar(position = position_dodge2(), stat="identity", width=.8, colour="black") +
  geom_hline(yintercept=33768, color = "#E64B35B2", linetype = "dashed", size=0.8) +
  geom_hline(yintercept=18944, color = "#7E6148B2", linetype = "dashed", size=0.8) +
  geom_text(aes(label=detected_proteins),angle=90,vjust=.25,hjust=1.2,size=2.5,color="black") +
  theme_bw() +
  xlab(label = "Samples") + ylab(label = "# identif. peptides") + 
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = "right",
        legend.text = element_text(size = 11),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 11, angle = 90, vjust=0.5),
        axis.text.y = element_text(size = 12) ) +
  scale_fill_manual(values = c("Ctrl"="#7E6148B2", "Ga"="#E64B35B2","Gc"="#F39B7FB2", "Wa"="#3C5488B2", "Wc"="#8491B4B1") ) 

# Summary of Quantified proteins and peptides
library(RColorBrewer)
RColorBrewer::display.brewer.all()

quantif <- data.frame(contrast=c("Gc/Ctrl", "Ga/Ctrl", "Wc/Ctrl", "Wa/Ctrl", "Gc/Ga", "Wc/Wa", "Wa/Ga", "Wc/Gc", "global"),
                      comparison=c("Inf_vs_Ctrl", "Inf_vs_Ctrl", "Inf_vs_Ctrl", "Inf_vs_Ctrl", "Inf_vs_Inf", "Inf_vs_Inf", "Inf_vs_Inf", "Inf_vs_Inf", "All_groups"),
                      proteins_dea=c(3275,  2892,  3435,  3127,  4563,  4963,  4369,  5209, 2724),
                      peptides_dea=c(14004, 11857, 15220, 13358, 26539, 32083, 25008, 33984, 11026),
                      qvalue_dea = c( 0.01, 0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01, NA),
                      log2FC_dea = c(0.579, 0.752, 0.561, 0.744, 0.572, 0.588, 0.662, 0.535, NA),
                      proteins_dd=c(  4363, 3861,  4685,  4201,  4278,  4655,  4114, 4708, NA),
                      proteins_zsc_dd=c(0,   0,     0,    0,     61,    59,    17,  40, NA)
)
saveRDS(object=quantif, file = "Fig2_identif_quantif/quantif.RDS")

# quantified proteins
library(RColorBrewer)
RColorBrewer::display.brewer.all()
#C <- 
quantif %>%
  ggplot( aes(x=reorder(contrast, proteins_dea), y=proteins_dea, fill=comparison)) + 
  geom_bar(position = position_dodge(width = 0.9), stat="identity", width=.8, colour="black") +
  geom_hline(yintercept = 2892, color = "#7E6148B2", linetype = "dashed", size=0.8) +
  geom_hline(yintercept = 4369, color = "#00A087B2", linetype = "dashed", size=0.8) +
  geom_text(aes(label=proteins_dea),angle=90,vjust=.3,hjust=1.2,size=3.5,color="black") +
  theme_bw() +
  xlab(label = "contrast") + ylab(label = "# quantif. proteins") +
  scale_fill_manual(values = c("Inf_vs_Inf"="#00A087B2", "Inf_vs_Ctrl"="#7E6148B2", "global"="gray", "All_groups"="gray") ) +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, vjust=1.3, hjust=1.2),
        axis.text.y = element_text(size = 12) )

# quantified peptides
#D <- 
quantif %>% 
  ggplot( aes(x=reorder(contrast,peptides_dea), y=peptides_dea, fill=comparison) ) +
  geom_bar(position = position_dodge(width = 0.9), stat="identity", width=.8, colour="black") +
  geom_hline(yintercept = 11857, color = "#7E6148B2", linetype = "dashed", size=0.8) +
  geom_hline(yintercept = 25008, color = "#00A087B2", linetype = "dashed", size=0.8) +
  geom_text(aes(label=peptides_dea),angle=90,vjust=.3,hjust=1.2,size=3.5,color="black") +
  theme_bw() +
  xlab(label = "Contrast") + ylab(label = "# quantif. peptides") + 
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, vjust=1.3, hjust=1.2),
        axis.text.y = element_text(size = 12) ) +
  scale_fill_manual(values = c("Inf_vs_Inf"="#00A087B2", "Inf_vs_Ctrl"="#7E6148B2", "global"="gray", "All_groups"="gray") ) 
#scale_fill_grey() + # library(ggsci)

