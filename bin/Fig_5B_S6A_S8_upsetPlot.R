library(tidyverse)
library(ComplexHeatmap)

#Importing DEA results from MSDAP for Upset_plot All
dea_GaCtrl <- readRDS("export_stats_genesummary/dea_GaCtrl.RDS")
dea_GcCtrl <- readRDS("export_stats_genesummary/dea_GcCtrl.RDS")
dea_WaCtrl <- readRDS("export_stats_genesummary/dea_WaCtrl.RDS")
dea_WcCtrl <- readRDS("export_stats_genesummary/dea_WcCtrl.RDS")

GaCtrl <- dea_GaCtrl$symbol
GcCtrl <- dea_GcCtrl$symbol
WaCtrl <- dea_WaCtrl$symbol
WcCtrl <- dea_WcCtrl$symbol

GaCtrl <- as.data.frame(GaCtrl)
GcCtrl <- as.data.frame(GcCtrl)
WaCtrl <- as.data.frame(WaCtrl)
WcCtrl <- as.data.frame(WcCtrl)

colnames(GaCtrl)[1] <- "Ga/Ctrl"
colnames(GcCtrl)[1] <- "Gc/Ctrl"
colnames(WaCtrl)[1] <- "Wa/Ctrl"
colnames(WcCtrl)[1] <- "Wc/Ctrl"

GaCtrl <- as.matrix(GaCtrl)
GcCtrl <- as.matrix(GcCtrl)
WaCtrl <- as.matrix(WaCtrl)
WcCtrl <- as.matrix(WcCtrl)

All <- list("Ga/Ctrl"=GaCtrl, 
            "Gc/Ctrl" = GcCtrl, 
            "Wa/Ctrl" = WaCtrl, 
            "Wc/Ctrl" = WcCtrl)

# Using ComplexHeatmap for upset plot (Fig 5B)
library(ComplexHeatmap)
All2 = make_comb_mat(All)
set_name(All2)
comb_name(All2)
comb_degree(All2)
set_size(All2)

UpSet(All2,
      set_order = c("Ga/Ctrl", "Gc/Ctrl", "Wa/Ctrl", "Wc/Ctrl"),
      pt_size = unit(5, "mm"), lwd = 2,
      #top_annotation = upset_top_annotation(All2, gp = gpar(col = comb_degree(All2))),
      #comb_col = c("red", "blue", "black")[comb_degree(All2)],
      bg_col = c("#F0F0FF", "#FFF0F0"), bg_pt_col = "#CCCCFF",
      top_annotation = upset_top_annotation(All2, 
                                            add_numbers = TRUE, 
                                            annotation_name_rot = 90,
                                            annotation_name_side = "left",
                                            axis_param = list(side = "left")),
      left_annotation = upset_left_annotation(All2, 
                                              add_numbers = TRUE, 
                                              ylim = c(0, 1300), 
                                              gp = gpar(fill = c("#E64B35", "#F39B7F", "#1F78B4", "#A6CEE3"), col = "black"), 
      )

################# Fig S6A Upset plot of less stringent Diff Detection identif-proteins of mean groups ####################
library(readxl)
identif <- read_excel("stringdb/identif_venn_3samp_per_group_perseus/matrix23.xlsx", 
                      sheet = "venn")
Ctrl <- as.data.frame(identif$Ctrl)
Ga <- as.data.frame(identif$Ga)
Gc <- as.data.frame(identif$Gc)
Wa <- as.data.frame(identif$Wa)
Wc <- as.data.frame(identif$Wc)

Ctrl <- na.omit(Ctrl)
Ga <- na.omit(Ga)
Gc <- na.omit(Gc)
Wa <- na.omit(Wa)
Wc <- na.omit(Wc)

colnames(Ctrl)[1] <- "Ctrl"
colnames(Ga)[1] <- "Ga"
colnames(Gc)[1] <- "Gc"
colnames(Wa)[1] <- "Wa"
colnames(Wc)[1] <- "Wc"

Ctrl <- as.data.frame(Ctrl)
Ga <- as.data.frame(Ga)
Gc <- as.data.frame(Gc)
Wa <- as.data.frame(Wa)
Wc <- as.data.frame(Wc)

Ctrl <- as.matrix(Ctrl)
Ga <- as.matrix(Ga)
Gc <- as.matrix(Gc)
Wa <- as.matrix(Wa)
Wc <- as.matrix(Wc)

identif <- list("Ctrl"=Ctrl, "Ga" = Ga, "Gc" = Gc, "Wa" = Wa, "Wc" = Wc)
a_vs_c <- list("Ga" = Ga, "Gc" = Gc, "Wa" = Wa, "Wc" = Wc)

library(ComplexHeatmap)
identif2 = make_comb_mat(identif)
UpSet(identif2,
      set_order = c("Ctrl", "Ga", "Gc", "Wa", "Wc"),
      pt_size = unit(5, "mm"), lwd = 2,
      #top_annotation = upset_top_annotation(identif2, gp = gpar(col = comb_degree(identif2))),
      #comb_col = c("red", "blue", "black")[comb_degree(identif2)],
      bg_col = c("#F0F0FF", "#FFF0F0"), bg_pt_col = "#CCCCFF",
      top_annotation = upset_top_annotation(identif2, 
                                            add_numbers = TRUE, 
                                            annotation_name_rot = 90,
                                            annotation_name_side = "left",
                                            axis_param = list(side = "left")),
      left_annotation = upset_left_annotation(identif2, 
                                              add_numbers = TRUE, 
                                              ylim = c(0, 6000), 
                                              gp = gpar(fill = c("#7E6148B2", "#E64B35", "#F39B7F", "#1F78B4", "#A6CEE3"), col = "black"), 
      )
)      

            
# Upsetplot_Upreg (Fig S8A)
dea_GaCtrl <- readRDS("export_stats_genesummary/dea_GaCtrl.RDS")
GaCtrl <- dea_GaCtrl$symbol[dea_GaCtrl$dysreg=="up"]
GcCtrl <- dea_GcCtrl$symbol[dea_GcCtrl$dysreg=="up"]
WaCtrl <- dea_WaCtrl$symbol[dea_WaCtrl$dysreg=="up"]
WcCtrl <- dea_WcCtrl$symbol[dea_WcCtrl$dysreg=="up"]
GaCtrl <- as.data.frame(GaCtrl)
GcCtrl <- as.data.frame(GcCtrl)
WaCtrl <- as.data.frame(WaCtrl)
WcCtrl <- as.data.frame(WcCtrl)
colnames(GaCtrl)[1] <- "Ga/Ctrl"
colnames(GcCtrl)[1] <- "Gc/Ctrl"
colnames(WaCtrl)[1] <- "Wa/Ctrl"
colnames(WcCtrl)[1] <- "Wc/Ctrl"
GaCtrl <- as.matrix(GaCtrl)
GcCtrl <- as.matrix(GcCtrl)
WaCtrl <- as.matrix(WaCtrl)
WcCtrl <- as.matrix(WcCtrl)
upreg <- list("Ga/Ctrl" = GaCtrl, 
              "Gc/Ctrl" = GcCtrl, 
              "Wa/Ctrl" = WaCtrl, 
              "Wc/Ctrl" = WcCtrl)
upreg2 = make_comb_mat(upreg)


library(ComplexHeatmap)
UpSet(upreg2,
      set_order = c("Ga/Ctrl", "Gc/Ctrl", "Wa/Ctrl", "Wc/Ctrl"),
      pt_size = unit(5, "mm"), lwd = 2,
      #top_annotation = upset_top_annotation(upreg2, gp = gpar(col = comb_degree(upreg2))),
      #comb_col = c("red", "blue", "black")[comb_degree(upreg2)],
      bg_col = c("#F0F0FF", "#FFF0F0"), bg_pt_col = "#CCCCFF",
      top_annotation = upset_top_annotation(upreg2, 
                                            add_numbers = TRUE, 
                                            annotation_name_rot = 90,
                                            annotation_name_side = "left",
                                            axis_param = list(side = "left")),
      left_annotation = upset_left_annotation(upreg2, 
                                              add_numbers = TRUE, 
                                              ylim = c(0, 1300), 
                                              gp = gpar(fill = c("#E64B35", "#F39B7F", "#1F78B4", "#A6CEE3"), col = "black"), 
      )
)

#Downregulated DAPs (Fig S8B)
GaCtrl <- dea_GaCtrl$symbol[dea_GaCtrl$dysreg=="down"]
GcCtrl <- dea_GcCtrl$symbol[dea_GcCtrl$dysreg=="down"]
WaCtrl <- dea_WaCtrl$symbol[dea_WaCtrl$dysreg=="down"]
WcCtrl <- dea_WcCtrl$symbol[dea_WcCtrl$dysreg=="down"]
GaCtrl <- as.data.frame(GaCtrl)
GcCtrl <- as.data.frame(GcCtrl)
WaCtrl <- as.data.frame(WaCtrl)
WcCtrl <- as.data.frame(WcCtrl)
colnames(GaCtrl)[1] <- "Ga/Ctrl"
colnames(GcCtrl)[1] <- "Gc/Ctrl"
colnames(WaCtrl)[1] <- "Wa/Ctrl"
colnames(WcCtrl)[1] <- "Wc/Ctrl"
GaCtrl <- as.matrix(GaCtrl)
GcCtrl <- as.matrix(GcCtrl)
WaCtrl <- as.matrix(WaCtrl)
WcCtrl <- as.matrix(WcCtrl)
downreg <- list("Ga/Ctrl" = GaCtrl, 
              "Gc/Ctrl" = GcCtrl, 
              "Wa/Ctrl" = WaCtrl, 
              "Wc/Ctrl" = WcCtrl)
downreg2 = make_comb_mat(downreg)

library(ComplexHeatmap)
UpSet(downreg2,
      set_order = c("Ga/Ctrl", "Gc/Ctrl", "Wa/Ctrl", "Wc/Ctrl"),
      pt_size = unit(5, "mm"), lwd = 2,
      #top_annotation = upset_top_annotation(downreg2, gp = gpar(col = comb_degree(downreg2))),
      #comb_col = c("red", "blue", "black")[comb_degree(downreg2)],
      bg_col = c("#F0F0FF", "#FFF0F0"), bg_pt_col = "#CCCCFF",
      top_annotation = upset_top_annotation(downreg2, 
                                            add_numbers = TRUE, 
                                            annotation_name_rot = 90,
                                            annotation_name_side = "left",
                                            axis_param = list(side = "left")),
      left_annotation = upset_left_annotation(downreg2, 
                                              add_numbers = TRUE, 
                                              ylim = c(0, 1300), 
                                              gp = gpar(fill = c("#7E6130", "#7E6148B2"), col = "black"), 
      )
)