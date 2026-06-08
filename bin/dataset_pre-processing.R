setwd("C:/Users/Laura/Downloads/astrocitos_COVID19/2nd_reanalysis/msdap_results/2024-10-09_23-34-50_combinedCTRL_msqrob_FDR0.01")
load("dataset.RData")
################ loading DEA & DDA results #################################################################
library(hgnc)
list_archives()
last_update()
latest_archive_url()
hgnc = import_hgnc_dataset(file = "C:/Users/Laura/Downloads/astrocitos_COVID19/2nd_reanalysis/msdap_results/2024-10-09_23-34-50_combinedCTRL_msqrob_FDR0.01/hgnc_complete_set_2024_10_09.txt")

library(msdap)
hgnc_table = hgnc_lookuptable("C:/Users/Laura/Downloads/astrocitos_COVID19/2nd_reanalysis/msdap_results/2024-10-09_23-34-50_combinedCTRL_msqrob_FDR0.01/hgnc_complete_set_2024_10_09.txt")
export_stats_genesummary(dataset,
                         gene_ambiguity = "prio_specific",
                         diffdetect_zscore_threshold = 6,
                         diffdetect_type = "auto",
                         dea_logfc_instead_of_effectsize = FALSE,
                         output_dir = "C:/Users/Laura/Downloads/astrocitos_COVID19/2nd_reanalysis/msdap_results/2024-10-09_23-34-50_combinedCTRL_msqrob_FDR0.01/export_stats_genesummary",
                         hgnc = hgnc_table,
                         xref = NULL,
                         remove_nohgnc = FALSE
)

#Replace contrast names in dataset
dataset$de_proteins$contrast <-  str_replace_all (dataset$de_proteins$contrast, 
                                                  c ("contrast: Ctrl vs Gc # condition_variable: group" = "Gc/Ctrl", 
                                                     "contrast: Ctrl vs Ga # condition_variable: group" = "Ga/Ctrl",
                                                     "contrast: Ctrl vs Wc # condition_variable: group" = "Wc/Ctrl",
                                                     "contrast: Ctrl vs Wa # condition_variable: group" = "Wa/Ctrl",
                                                     "contrast: Ga vs Gc # condition_variable: group" = "Gc/Ga",
                                                     "contrast: Wa vs Wc # condition_variable: group" = "Wc/Wa",
                                                     "contrast: Gc vs Wc # condition_variable: group" = "Wc/Gc",
                                                     "contrast: Ga vs Wa # condition_variable: group" = "Wa/Ga") )

dataset$dd_proteins$contrast <-  str_replace_all (dataset$dd_proteins$contrast, 
                                                  c ("contrast: Ctrl vs Gc # condition_variable: group" = "Gc/Ctrl", 
                                                     "contrast: Ctrl vs Ga # condition_variable: group" = "Ga/Ctrl",
                                                     "contrast: Ctrl vs Wc # condition_variable: group" = "Wc/Ctrl",
                                                     "contrast: Ctrl vs Wa # condition_variable: group" = "Wa/Ctrl",
                                                     "contrast: Ga vs Gc # condition_variable: group" = "Gc/Ga",
                                                     "contrast: Wa vs Wc # condition_variable: group" = "Wc/Wa",
                                                     "contrast: Gc vs Wc # condition_variable: group" = "Wc/Gc",
                                                     "contrast: Ga vs Wa # condition_variable: group" = "Wa/Ga") )

#Read results of export_stats_genesummary
dea_GaCtrl <- read_excel("export_stats_genesummary/Ctrl_vs_Ga___condition_variable__group__msqrob.xlsx")
dea_GcCtrl <- read_excel("export_stats_genesummary/Ctrl_vs_Gc___condition_variable__group__msqrob.xlsx")
dea_WaCtrl <- read_excel("export_stats_genesummary/Ctrl_vs_Wa___condition_variable__group__msqrob.xlsx")
dea_WcCtrl <- read_excel("export_stats_genesummary/Ctrl_vs_Wc___condition_variable__group__msqrob.xlsx")

#Replace contrast names in DEA results
dea_GaCtrl$contrast <- str_replace(dea_GaCtrl$contrast, "contrast: Ctrl vs Ga # condition_variable: group", "Ga/Ctrl")
dea_GcCtrl$contrast <- str_replace(dea_GcCtrl$contrast, "contrast: Ctrl vs Gc # condition_variable: group", "Gc/Ctrl")
dea_WaCtrl$contrast <- str_replace(dea_WaCtrl$contrast, "contrast: Ctrl vs Wa # condition_variable: group", "Wa/Ctrl")
dea_WcCtrl$contrast <- str_replace(dea_WcCtrl$contrast, "contrast: Ctrl vs Wc # condition_variable: group", "Wc/Ctrl")

#up/down-regulated proteins
dea_GaCtrl <- dea_GaCtrl %>% dplyr::filter(signif=="TRUE" & score_type=="dea")
dea_GcCtrl <- dea_GcCtrl %>% dplyr::filter(signif=="TRUE" & score_type=="dea")
dea_WaCtrl <- dea_WaCtrl %>% dplyr::filter(signif=="TRUE" & score_type=="dea")
dea_WcCtrl <- dea_WcCtrl %>% dplyr::filter(signif=="TRUE" & score_type=="dea")

#Add dysreg column: up/down regulated, saveRDS, write excell
dea_GaCtrl$dysreg[dea_GaCtrl$log2fc > 0] <- "up"
dea_GaCtrl$dysreg[dea_GaCtrl$log2fc < 0] <- "down"
saveRDS(object=dea_GaCtrl, file = "export_stats_genesummary/dea_GaCtrl.RDS")
write_xlsx(dea_GaCtrl, "export_stats_genesummary/dea_GaCtrl.xlsx")

dea_GcCtrl$dysreg[dea_GcCtrl$log2fc > 0] <- "up"
dea_GcCtrl$dysreg[dea_GcCtrl$log2fc < 0] <- "down"
saveRDS(object=dea_GcCtrl, file = "export_stats_genesummary/dea_GcCtrl.RDS")
write_xlsx(dea_GcCtrl, "export_stats_genesummary/dea_GcCtrl.xlsx")

dea_WaCtrl$dysreg[dea_WaCtrl$log2fc > 0] <- "up"
dea_WaCtrl$dysreg[dea_WaCtrl$log2fc < 0] <- "down"
saveRDS(object=dea_WaCtrl, file = "export_stats_genesummary/dea_WaCtrl.RDS")
write_xlsx(dea_WaCtrl, "export_stats_genesummary/dea_WaCtrl.xlsx")

dea_WcCtrl$dysreg[dea_WcCtrl$log2fc > 0] <- "up"
dea_WcCtrl$dysreg[dea_WcCtrl$log2fc < 0] <- "down"
saveRDS(object=dea_WcCtrl, file = "export_stats_genesummary/dea_WcCtrl.RDS")
write_xlsx(dea_WcCtrl, "export_stats_genesummary/dea_WcCtrl.xlsx")

dataset$de_proteins <- dataset$de_proteins %>% left_join(dataset$proteins, by=c("protein_id"="protein_id"))
dataset$de_proteins$updown = ifelse(dataset$de_proteins$foldchange.log2 < 0, "down regulated", "up regulated")
dataset$de_proteins$updown[dataset$de_proteins$signif != TRUE] = "not significant"