# TRIAGE-Cluster

TRIAGE-Cluster is designed to analyse single-cell data to distinguish cell identity using epigenetic information as biological anchor.

KEY FILES:

1. ensembl_id_mapping.txt - ensembl id and gene name for different species
2. repressive_hg19_cont_prop.txt - Repressive tendency score (RTS) file
3. Priority_epimap_rts.txt - The priority RTS gene generated from EpiMap data
4. HumanMouse_geneEnsembl_DB.txt - Mapping between mouse and human gene name and ensembl id


CODE DESCRIPTION:

1. Preprocessing_and_TRIAGE_transformation.R - Pre-Processing, Normalisation and transform expression matrix to discordance matrix
2. RTS_assign.R - Identify the RTS priority gene with highest expression and corresponding RTS for each cell
3. TRIAGE_Cluster.py - Peaks generation

If you use TRIAGE-Cluster in your research, please cite our paper: 
Sun, Y., Shim, W., et al. "Inferring cell diversity in single cell data using consortium-scale epigenetic data as a biological anchor for cell identity" *Nucleic Acids Research*, Volume 51, Issue 11, 23 June 2023, https://doi.org/10.1093/nar/gkad307.

See also https://github.com/palpant-comp/TRIAGE_R_Package for the R package version of this method. 
