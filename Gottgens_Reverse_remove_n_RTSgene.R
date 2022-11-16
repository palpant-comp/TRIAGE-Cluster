library(Seurat)
library(ggplot2)
library(dplyr)
library(mclust)
library(tidyr)

g.meta <- readRDS("/home/uqysun19/90days/manuscript_TRIAGEclustering/Gottgens_geneassign.RDS")
v1.all <- readRDS("/home/uqysun19/90days/manuscript_TRIAGEclustering/epimap_v1TRIAGE_genelist_expressed_in_Gottgens.RDS")
gottgensTRIAGE <- readRDS("/home/uqysun19/90days/manuscript_TRIAGEclustering/Gottgens_v1TRIAGE.RDS")

#get the common genes
a <- intersect(rownames(gottgensTRIAGE), rownames(v1.all))

#remove genes not detected in expression matrix from the gene list 
v1.new <- v1.all[a,]
#order based on RTS score from high to low
v1.new <- v1.new[order(v1.new$RTS, decreasing = T),]
#only choose priority genes
v1.p <- v1.new[v1.new$Priority == "Y",]
#make the rownames order in expression matrix same as the gene list
gottgens.new <- gottgensTRIAGE[rownames(v1.p),]

#remove top n number of genes from RTS gene list
n<-c(1:15)*50

# iteration <- 100

ARI.table <- matrix(nrow = 1, ncol = length(n), 
                    dimnames = list("ARI", paste0("remove_", n)))

# v1.newsub <- v1.new[which(v1.new$Priority == "Y"),]
 
v1gene.list1 <- rownames(v1.p)

#make the order same as the v1gene.list1
g.meta <- g.meta[order(g.meta$RTS, decreasing = F), ]

# v1gene.list2 <- unique(g.meta$v1gene)

set.seed(0)
for (j in 1:length(n)){
    
    v1gene.list <- head(v1gene.list1, -n[j])

    gottgens.rts <- gottgens.new

    g.meta.new <- g.meta[,c("cell_name","v1gene")]

    for (i in 1:length(v1gene.list)){
        if (class(gottgens.rts) != 'data.frame') {
            if (gottgens.rts[i] > 0){
                g.meta.new$v1gene_new <- g.meta.new$v1gene_new %>% replace_na(v1gene.list[i])
                break
            }
            
        } else{
            z <- gottgens.rts[v1gene.list[i],] > 0
            a <- colnames(gottgens.rts)[z]
            g.meta.new[which(g.meta.new$cell_name %in% a), "v1gene_new"] <- v1gene.list[i]
            gottgens.rts <- gottgens.rts[, !z]
        }
    }
        if (sum(is.na(g.meta.new[, "v1gene_new"]) == TRUE) >= 1){
        g.meta.new[, "v1gene_new"][is.na(g.meta.new[, "v1gene_new"])] <- "GeneX"
    } else{
        ARI <- mclust::adjustedRandIndex(g.meta.new[, "v1gene_new"], g.meta.new[, "v1gene"])
        ARI.table[1,  j] <- ARI
    }
}

saveRDS(ARI.table, "/home/uqysun19/90days/manuscript_TRIAGEclustering/reverse_RTSgene_test/Gottgens_ARI_reverse_remove_top_RTS.RDS")



