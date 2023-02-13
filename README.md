# TRIAGE-Cluster
TRIAGE-Cluster is designed to analyse single-cell data to distinguish cell identity using epigenetic information as biological anchor.

KEY FILES:

ensembl_id_mapping.txt - ensembl id and gene name for different species
repressive_hg19_cont_prop.txt - Repressive tendency score (RTS) file
Priority_epimap_rts.txt - The priority RTS gene generated from EpiMap data
HumanMouse_geneEnsembl_DB.txt - Mapping between mouse and human gene name and ensembl id


CODE DESCRIPTION:

Preprocessing_and_TRIAGE_transformation.R - Pre-Processing, Normalisation and transform expression matrix to discordance matrix
RTS_assign.R - Identify the RTS priority gene with highest expression and corresponding RTS for each cell
TRIAGE_Cluster.py - Peaks generation
