library(dplyr)
library(tidyr)
library(parallel)

'%!in%' <- function(x,y)!('%in%'(x,y))

dir <- "/home/uqysun19/90days/manuscript_TRIAGEclustering/bottom2top_RTS/Gottgens_data/"

args <- commandArgs(trailingOnly = TRUE)

n <- as.numeric(args[1]) * 50

ari <- function(n){
    whole.pick <- rownames(df)[-sample(1:length(rownames(df)),n)]
    
    v1gene.list <- intersect(v1gene.list1, whole.pick)

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


args <- commandArgs(trailingOnly = TRUE)

n <- as.numeric(args[1]) * 50

meta <- readRDS(paste0(dir, "Gottgens_geneassign.RDS"))
v1 <- readRDS(paste0(dir, "Gottgens_v1genelist_new.RDS"))
df <- readRDS(paste0(dir, "Gottgens_v1TRIAGE_matrix.RDS"))

v1gene.list1 <- rownames(v1)
df.new <- df[rownames(v1),]
meta <- meta[order(meta$RTS, decreasing = T), ]

iteration <- c(1:100)
ARI.table <- matrix(nrow = 1, ncol = length(iteration), dimnames = list("ARI", paste0("iteration", iteration)))

for (x in iteration){
  ari.t <- mclapply(n, ari, mc.cores=10)
  ARI.table[1,  x] <- ari.t[[1]]
  saveRDS(ARI.table, paste0(dir, "wholegenome_random_remove", n, "_ARItable.RDS"))
}


