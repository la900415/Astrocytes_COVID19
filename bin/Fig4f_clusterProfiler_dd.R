####### Fig 4f. ClusterProfiler of diff detection ########################################################
#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
#BiocManager::install("clusterProfiler")
library(clusterProfiler) 
#install.packages("GOplot")
library(GOplot)
#BiocManager::install("DOSE")
library(DOSE)
#BiocManager::install("enrichplot")
library(enrichplot)

ck_ego_dd <- compareCluster(gene_symbols~contrast+updown, data=dd_GcGa_WcWa, fun="enrichGO", ont="BP",
                            pvalueCutoff=0.05, keyType="SYMBOL", OrgDb=org.Hs.eg.db, pAdjustMethod="BH", 
                            readable=TRUE)
ck_ego_dd <- setReadable(ck_ego_dd, OrgDb = org.Hs.eg.db, keyType="SYMBOL")
saveRDS(object=ck_ego_dd, file = "Fig3_clusterprofiler/ck_ego_dd.RDS")

ck_ego_dd_red <- simplify(ck_ego_dd, cutoff=0.7, by="p.adjust", select_fun=min, measure="Wang", semData=NULL)
saveRDS(object=ck_ego_dd_red, file = "Fig3_clusterprofiler/ck_ego_dd_red.RDS")

dotplot(ck_ego_dd_red, x="contrast", showCategory=5,includeAll=T,label_format=40)+facet_grid(~updown)+
  theme(axis.text.x = element_text(size = 14, angle=45, hjust=1.1, vjust=1.2),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        strip.text = element_text(size = 16) )

cnetplot(ck_ego_dd, showCategory=5, label_format=30, cex_label_gene = 1.5, cex_label_category=1.5) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14) )
