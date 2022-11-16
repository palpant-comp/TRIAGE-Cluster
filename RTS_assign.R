library(Seurat)
library(ggplot2)
library(Matrix)
library(patchwork)
library(dplyr)
library(tidyr)


v1p <- readRDS("/clusterdata/uqysun19/data/manuscript_TRIAGEclustering/Priority_epimap_rts.RDS")
df <- readRDS(paste0(dir, "Gottgens_norm_data.RDS"))
v1 <- v1p[v1p[,"Gene"] %in% rownames(df),]
saveRDS(v1, paste0(dir, "Gottgens_v1genelist.RDS"))

#for each cell, get the RTS priority gene expressed and assign the corresponding RTS as a quantatitive feature for cell-type-specific cells
gene_assign <- function(v1, tdf, g.meta){
    rownames(v1) <- v1$Gene
    v1gene.list <- rownames(v1)
    tdf.sub <- tdf[v1gene.list,]
    wu.rts <- tdf.sub

    for (i in 1:length(v1gene.list)){
            if (class(wu.rts) != 'data.frame') {
                if (wu.rts[i] > 0){
                    g.meta$v1gene <- g.meta$v1gene %>% replace_na(v1gene.list[i])
                }
            } else {
                z <- wu.rts[v1gene.list[i],] > 0
                c <- colnames(wu.rts)[z]
                g.meta[which(rownames(g.meta) %in% c), "v1gene"] <- v1gene.list[i]
                wu.rts <- wu.rts[, !z]
            }
        }

    g.meta$RTS <- v1[g.meta$v1gene, "RTS"]
    return(g.meta)
}


tdf <- readRDS(paste0(dir, fn[n], "_data/", fn[n], "_v1TRIAGE_matrix.RDS"))
g.meta <- readRDS(paste0(dir, fn[n], "_data/", fn[n], "_v1TRIAGE_seu_metadata.RDS"))

gn <- gene_assign(v1, tdf, g.meta)
saveRDS(gn, paste0(dir, "Gottgens_geneassign.RDS"))
#csv file as input for python script for TRIAGE-Cluster (density estimation and peak selection)
write.csv(gn, paste0(dir, "Gottgens_geneassign.csv"), rownames = TRUE)
