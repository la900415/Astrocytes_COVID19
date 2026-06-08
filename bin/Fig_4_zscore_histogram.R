
setwd("C:/Users/Laura/Downloads/astrocitos_COVID19/2nd_reanalysis/msdap_results/2024-10-09_23-34-50_combinedCTRL_msqrob_FDR0.01")
load("dataset.RData")
dataset$dd_proteins <- dataset$dd_proteins %>% left_join(dataset$proteins, by=c("protein_id"="protein_id"))

########### zscore of DDA 

###### Gc _vs_Ga ######################################################################
dataset$dd_proteins %>% filter(contrast=="contrast: Ga vs Gc # condition_variable: group") %>% 
  filter(pass_filters=="TRUE") %>%
  ggplot( aes(x=zscore) ) +
  geom_histogram(binwidth = 3, fill = "gray", color = "black") +
  theme(panel.grid.major.y = element_line(color = "gray", size = 0.7, linetype = "dotted"),
        panel.grid.major.x = element_line(color = "gray", size = 0.7, linetype = "dotted"),
        panel.background   = element_rect(fill = "white", colour = "grey50"),
        legend.box.background = element_blank(),
        legend.box.margin = margin(6,6,6,6),
        legend.key = element_blank() ) + #insert vertical line x=6 and x=-6 with red color
  geom_vline(xintercept = c(6, -6), color = "red", linetype = "dashed") +
  xlab(label = "Differential-detected z−score") +
  ylab(label = "Protein count") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 10) ) +
  labs(title = "#proteins: 4278, #|z-score| >=6: 61", subtitle = NULL, caption = NULL, tag = "") +
  annotate("text", x = -12.5, y = 3000, label = "Ga (9)", hjust = 0.8, vjust = 0, size = 6, fontface="bold", colour="#E64B35B2") +
  annotate("text", x = 12.5, y = 3000, label = "Gc (52)",   hjust = 0.5, vjust = 0, size = 6, fontface="bold", colour="#F39B7FB2")

###### Wc _vs_Wa ######################################################################
dataset$dd_proteins %>% filter(contrast=="contrast: Wa vs Wc # condition_variable: group") %>% 
  filter(pass_filters=="TRUE") %>%
  ggplot( aes(x=zscore) ) +
  geom_histogram(binwidth = 2.5, fill = "gray", color = "black") +
  theme(panel.grid.major.y = element_line(color = "gray", size = 0.7, linetype = "dotted"),
        panel.grid.major.x = element_line(color = "gray", size = 0.7, linetype = "dotted"),
        panel.background   = element_rect(fill = "white", colour = "grey50"),
        legend.box.background = element_blank(),
        legend.box.margin = margin(6, 6, 6, 6),
        legend.key = element_blank() ) + #insert vertical line x=6 and x=-6 with red color
  geom_vline(xintercept = c(6, -6), color = "red", linetype = "dashed") +
  xlab(label = "Differential-detected z−score") +
  ylab(label = "Protein count") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 10) ) +
  labs(title = "#proteins: 4655, #|z-score| >=6: 59", subtitle = NULL, caption = NULL, tag = "") +
  annotate("text", x = -8, y = 3000, label = "Wa (13)", hjust = 0.8, vjust = 0, size = 6, fontface="bold", colour="blue3") +
  annotate("text", x = 12.5, y = 3000, label = "Wc (46)",   hjust = 0.5, vjust = 0, size = 6, fontface="bold", colour="cornflowerblue")

###### Wa _vs_Ga ######################################################################
dataset$dd_proteins %>% filter(contrast=="contrast: Ga vs Wa # condition_variable: group") %>% 
  filter(pass_filters=="TRUE") %>%
  ggplot( aes(x=zscore) ) +
  geom_histogram(binwidth = 2.8, fill = "gray", color = "black") +
  theme(panel.grid.major.y = element_line(color = "gray", size = 0.7, linetype = "dotted"),
        panel.grid.major.x = element_line(color = "gray", size = 0.7, linetype = "dotted"),
        panel.background   = element_rect(fill = "white", colour = "grey50"),
        legend.box.background = element_blank(),
        legend.box.margin = margin(6, 6, 6, 6),
        legend.key = element_blank() ) + #insert vertical line x=6 and x=-6 with red color
  geom_vline(xintercept = c(6, -6), color = "red", linetype = "dashed") +
  xlab(label = "Differential-detected z−score") +
  ylab(label = "Protein count") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 10) ) +
  labs(title = "#proteins: 4114, #|z-score| >=6: 17", subtitle = NULL, caption = NULL, tag = "") +
  annotate("text", x = -7, y = 3000, label = "Ga (1)", hjust = 0.8, vjust = 0, size = 6, fontface="bold", colour="#E64B35B2") +
  annotate("text", x = 9, y = 3000, label = "Wa (16)", hjust = 0.5, vjust = 0, size = 6, fontface="bold", colour="blue3")

###### Wc _vs_Gc ######################################################################
dataset$dd_proteins %>% filter(contrast=="contrast: Gc vs Wc # condition_variable: group") %>% 
  filter(pass_filters=="TRUE") %>%
  ggplot( aes(x=zscore) ) +
  geom_histogram(binwidth = 2.8, fill = "gray", color = "black") +
  theme(panel.grid.major.y = element_line(color = "gray", size = 0.7, linetype = "dotted"),
        panel.grid.major.x = element_line(color = "gray", size = 0.7, linetype = "dotted"),
        panel.background   = element_rect(fill = "white", colour = "grey50"),
        legend.box.background = element_blank(),
        legend.box.margin = margin(6, 6, 6, 6),
        legend.key = element_blank() ) + #insert vertical line x=6 and x=-6 with red color
  geom_vline(xintercept = c(6, -6), color = "red", linetype = "dashed") +
  xlab(label = "Differential-detected z−score") +
  ylab(label = "Protein count") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 10) ) +
  labs(title = "#proteins: 4708, #|z-score| >=6: 40", subtitle = NULL, caption = NULL, tag = "") +
  annotate("text", x = -8, y = 3000, label = "Gc (14)", hjust = 0.8, vjust = 0, size = 6, fontface="bold", colour="#F39B7FB2") +
  annotate("text", x = 10, y = 3000, label = "Wc (26)", hjust = 0.5, vjust = 0, size = 6, fontface="bold", colour="cornflowerblue")
