library(dplyr)
library(tidyr)
library(parallel)

'%!in%' <- function(x,y)!('%in%'(x,y))

dir <- "/home/uqysun19/90days/manuscript_TRIAGEclustering/bottom2top_RTS/Gottgens_data/"



ari <- function(j){
    
    v1gene.list <- head(v1gene.list1, -j)

    df.rts <- df.new

    meta.new <- meta[,c("cell_name","v1gene")]

    for (i in 1:length(v1gene.list)){
        if (class(df.rts) != 'data.frame') {
            if (df.rts[i] > 0){
                meta.new$v1gene_new <- meta.new$v1gene_new %>% replace_na(v1gene.list[i])
                break
            }
            
        } else {
            z <- df.rts[v1gene.list[i],] > 0
            a <- colnames(df.rts)[z]
            meta.new[which(meta.new$cell_name %in% a), "v1gene_new"] <- v1gene.list[i]
            df.rts <- df.rts[, !z]
        }
    }
        if (sum(is.na(meta.new[, "v1gene_new"]) == TRUE) >= 1){
            meta.new[, "v1gene_new"][is.na(meta.new[, "v1gene_new"])] <- "GeneX"
        }

        ARI <- mclust::adjustedRandIndex(meta.new[, "v1gene_new"], meta.new[, "v1gene"])
}


meta <- readRDS(paste0(dir, "Gottgens_geneassign.RDS"))
v1 <- readRDS(paste0(dir, "Gottgens_v1genelist_new.RDS"))
df <- readRDS(paste0(dir, "Gottgens_v1TRIAGE_matrix.RDS"))

v1gene.list1 <- rownames(v1)
df.new <- df[rownames(v1),]
meta <- meta[order(meta$RTS, decreasing = T), ]

a <- c(1:14)*50
ari.t <- mclapply(a, ari, mc.cores=10)
ARI.table <- data.frame(ARI = as.character(ari.t), sep = a, type = "removeBOTTOM", row.names = paste0("remove_", a))
saveRDS(ARI.table, paste0(dir, "Gottgens_rmBOTTOMsep50_ARItable.RDS"))















