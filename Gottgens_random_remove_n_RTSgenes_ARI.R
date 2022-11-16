library(dplyr)
library(tidyr)
library(parallel)

'%!in%' <- function(x,y)!('%in%'(x,y))

dir <- "/home/uqysun19/90days/manuscript_TRIAGEclustering/bottom2top_RTS/Gottgens_data/"



ari <- function(n){
    
    v1gene.list <- v1gene.list1[-sample(1:length(v1gene.list1),n)]

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
  saveRDS(ARI.table, paste0(dir, "random_remove", n, "_ARItable.RDS"))
}










# library(Seurat)
# library(ggplot2)
# library(dplyr)
# library(mclust)
# library(tidyr)

# args <- commandArgs(trailingOnly = TRUE)

# n <- as.numeric(args[1]) * 50

# g.meta <- readRDS("/home/uqysun19/90days/manuscript_TRIAGEclustering/Gottgens_geneassign.RDS")
# v1.all <- readRDS("/home/uqysun19/90days/manuscript_TRIAGEclustering/epimap_v1TRIAGE_genelist_expressed_in_Gottgens.RDS")
# gottgensTRIAGE <- readRDS("/home/uqysun19/90days/manuscript_TRIAGEclustering/Gottgens_v1TRIAGE.RDS")

# #get the common genes
# a <- intersect(rownames(gottgensTRIAGE), rownames(v1.all))

# #remove genes not detected in expression matrix from the gene list 
# v1.new <- v1.all[a,]
# #order based on RTS score from high to low
# v1.new <- v1.new[order(v1.new$RTS, decreasing = T),]
# #only choose priority genes
# v1.p <- v1.new[v1.new$Priority == "Y",]
# #make the rownames order in expression matrix same as the gene list
# gottgens.new <- gottgensTRIAGE[rownames(v1.p),]


# iteration <- 100

# ARI.table <- matrix(nrow = 1, ncol = iteration, dimnames = list("ARI", paste0("iteration", c(1:iteration))))

# # v1.newsub <- v1.new[which(v1.new$Priority == "Y"),]

# v1gene.list1 <- rownames(v1.p)

# g.meta <- g.meta[order(g.meta$RTS, decreasing = T), ]

# #v1gene.list2 <- unique(g.meta$v1gene)


# set.seed(0)
# for (j in 1:iteration){
#   v1gene.list <- v1gene.list1[-sample(1:length(v1gene.list1),n)]
  
#   gottgens.rts <- gottgens.new
  
#   g.meta.new <- g.meta[,c("cell_name","v1gene")]
  
#   for (i in 1:length(v1gene.list)){
#     if (class(gottgens.rts) != 'data.frame') {
#       if (gottgens.rts[i] > 0){
#         g.meta.new$v1gene_new <- g.meta.new$v1gene_new %>% replace_na(v1gene.list[i])
#         break
#       }
#     } else {
#       z <- gottgens.rts[v1gene.list[i],] > 0
#       a <- colnames(gottgens.rts)[z]
#       g.meta.new[which(g.meta.new$cell_name %in% a), "v1gene_new"] <- v1gene.list[i]
#       gottgens.rts <- gottgens.rts[, !z]
#       print(z)
#     }
#   }
  
#     ARI <- mclust::adjustedRandIndex(g.meta.new[, "v1gene_new"], g.meta.new[, "v1gene"])
#     ARI.table[1,  j] <- ARI
# }

# saveRDS(ARI.table, paste0("/home/uqysun19/90days/manuscript_TRIAGEclustering/random_RTSgene_test/random_remove_", n, "_ARI.RDS"))

