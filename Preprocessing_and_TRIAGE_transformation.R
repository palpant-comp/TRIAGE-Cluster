library(data.table)
library(dplyr)
library(tidyr)
library(scran)
library(scater)
library(Seurat)
library(ggplot2)
library(Matrix)
library(patchwork)

'%!in%' <- function(x,y)!('%in%'(x,y))

#convert mouse gene name to human gene name use database from Ensembl
#input is raw count matrix (row is gene, column is cell.)
dir <- "/home/uqysun19/90days/manuscript_TRIAGEclustering/bottom2top_RTS/"

idcon <- function(m, m_names = NULL, species = NULL, geneid = NULL, dir = NULL){
    
    gene_sum <- function(gn = NULL){
    db1 <- hm.db[hm.db[,gn] %in% m_names,]
    m[,gn] <- rownames(m)
    mdf <- merge(db1[c(gn, "Human.gene.name")], m, by = gn)
    mdf[,gn] <- NULL
    mdf1 <-mdf[duplicated(mdf["Human.gene.name"])|duplicated(mdf["Human.gene.name"], fromLast=TRUE),]
    mdf2 <- as.data.frame(mdf1 %>% group_by(Human.gene.name) %>% summarise_all(sum))
    mdf3 <- mdf[mdf[,"Human.gene.name"] %!in% mdf1[,"Human.gene.name"],]
    mdf4 <- rbind(mdf2, mdf3)
    rownames(mdf4) <- mdf4[,"Human.gene.name"]
    mdf4[,"Human.gene.name"] <- NULL    
    return(mdf4)
}
    
	if (length(grep('^w$', ls(.GlobalEnv))) == 0){
        hm.db <<- readRDS('/home/uqysun19/90days/manuscript_TRIAGEclustering/HumanMouse_geneEnsembl_DB.RDS')
	}
	if (length(m_names) == 0){m_names <- rownames(m)}

	#check species of input data
	if (species == "Mouse") {
		if (geneid == TRUE){
            out <- gene_sum(gn = "Gene.name")

			} else {
            out <- gene_sum(gn = "Gene.stable.ID")

			}
	} else {
		if (geneid == TRUE){
			gene_name <- gsub("\\--.*", "", rownames(m))
            m$gene <- gene_name
            out <- as.data.frame(m %>% group_by(gene) %>% summarise_all(sum))
            rownames(out) <- out[,"gene"]
            out[,"gene"] <- NULL
            
        } else {
            out <- gene_sum(gn = "Human.gene.stable.ID")

        }
	  }
    return (out)
    }
 
#remove genes that do not express in any cells
#use scran to perform cell filtering and normalisation
#the output normalised matrix is not log transformed because scran usng log2 
#and we perform log10 transformation later in TRIAGE transformation

gene_filter <- function(df){
    Toremove <- which(rowSums(df)==0)
    if (length(Toremove) == 0){
        df.filter <- df    
    } else{
        df.filter <- df[-Toremove,]
    }
    return(df.filter)
} 

cell_filter <- function(df){
   sce <- SingleCellExperiment(assays=list(counts=df))
   qcstats <- perCellQCMetrics(sce)
   qcfilter <- quickPerCellQC(qcstats)
   sce.filter <- sce[,!qcfilter$discard] 
   return(sce.filter)
}


scran_norm <- function(df){
    sce <- computeSumFactors(df)
    #log2 in scran
    sce.norm <- logNormCounts(sce, log = FALSE)
#     sce.norm <- as.data.frame(assay(sce,"normcounts"))
    return(sce.norm)
}

#TRIAGE transformation generate discordance matrix using # original expression matrix (normalised matrix) 
#The normalised matrix could be log transformed or not, # if not, put log = FALSE)
dir <- "/clusterdata/uqysun19/data/manuscript_TRIAGEclustering/bottom2top_RTS/"

TRIAGE1 <- function(m, m_names = NULL, species = NULL, geneid = NULL, 
	log = NULL, dir = NULL){
	if (length(grep('^w$', ls(.GlobalEnv))) == 0){
		w <<- as.data.frame(fread(paste0(dir, 'repressive_hg19_cont_prop.txt')))
        g <<- as.data.frame(fread(paste0(dir, 'ensembl_id_mapping.txt')))
	}
	if (length(m_names) == 0){m_names <- row.names(m)}

	#check species of input data
	if (species == "Human") {
		s.col <- "Gene"
	} else{
		s.col <- paste(species, "gene")
	}

    #check ensembl id or gene names of input data gene 
	if (geneid == TRUE){
    g.col <- "name"
    } else{
    g.col <- "stable ID"
    }
    
    #get the colnames should be used to select column from reference data "ensembl_id_mapping.txt"
    select.col <- paste(s.col, g.col)

    colfn <- function(x){return(g[which(x == g[, select.col]), ][1,2])}
    hgname <- sapply(m_names, colfn)

	mappable <- which(hgname %in% w$V1)
	l <- rep(0, nrow(m))

	RTSfn <- function(gene){return(as.numeric(w[which(w$V1 == gene),][1,2]))}
	l[mappable] <- sapply(hgname[mappable], RTSfn)
    
    #check if input normalised data is log transformed or not
    if (log == TRUE){
    	out <- log(m+1)*l

    } else {out <- m*l}

	# row.names(out) <- m_names
	return(out)
}

#create seurat object and perform dimension reduction
seuob <- function(df, me, dfn = NULL){
	seu <- CreateSeuratObject(df, meta.data = me)
    seu$stim <- dfn
    seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = nrow(df), verbose = F)
    seu <- ScaleData(seu, verbose = F)
    seu <- RunPCA(seu, npcs = 50, verbose = F)
    seu <- RunUMAP(seu, reduction = "pca", dims = 1:50, 
                    umap.method = "umap-learn", metric = "correlation", n.components = 2L)

    return(seu)
}

#extract metadata and add UMAP coordinates
metadf <- function(seu){
	g.meta <- cbind(seu@meta.data, as.data.frame(Embeddings(seu[["umap"]])))
    g.meta$cell_name <- rownames(g.meta) 
    return(g.meta)
}






#example
df <- readRDS(paste0(dir, "Gottgens_raw_data.RDS"))
tdf <- idcon(df, species = "Mouse", geneid = TRUE, dir = dir)
df_filter <- gene_filter(tdf)
df_filter1 <- cell_filter(df_filter)
df.norm <- scran_norm(df_filter1)
sce.norm <- as.data.frame(assay(df.norm, "normcounts"))
#species choose "Human" here because the gene name is converted to human gene name at the first gene conversion step
Tdf <- TRIAGE1(sce.norm, species = "Human", geneid = TRUE, log = TRUE, dir = dir)

#make sure the cell number is same as the cell filtered matrix
me <- readRDS(paste0(dir, "Gottgens_metadata.RDS"))

seu <- seuob(Tdf, me)
seu.me <- metadf(seu)

saveRDS(tdf, paste0(dir, "Gottgens_geneconsum_raw_data.RDS"))
saveRDS(df.norm, paste0(dir, "Gottgens_norm_SCEob.RDS"))
saveRDS(sce.norm, paste0(dir, "Gottgens_norm_data.RDS")) 
saveRDS(df_filter, paste0(dir, "Gottgens_genefiltered_data.RDS"))    
saveRDS(df_filter1, paste0(dir, "Gottgens_cellfiltered_SCEob.RDS"))    
saveRDS(Tdf, paste0(dir, "Gottgens_v1TRIAGE_matrix.RDS"))
saveRDS(seu, paste0(dir, "Gottgens_v1TRIAGE_seu.RDS"))
saveRDS(seu.me, paste0(dir, "Gottgens_v1TRIAGE_seu_metadata.RDS"))


















