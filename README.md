# CCLE_integrative

Integrating CCLE RNA-seq and RPPA data. 

Data using:
  1) raw RNA-seq gene counts*
  2) TPM normalized RNA-seq*
  3) RPKM normalized RNA-seq
  3) normalized RPPA*
* Used these datasets for analysis

Data was obtained from CCLE: https://portals.broadinstitute.org/ccle

Names of datasets used (based on v3.1 script):
  1) CCLE_RNAseq_genes_counts_20180929.gct (found in CCLE site)
      - Raw RNA-seq counts
      - used for PCA

  2) CCLE_RNAseq_rsem_genes_tpm_20180929.txt (found in CCLE site)
      - RNA-seq counts TPM normalized

  3) breast_cell_lines_V2.csv (found in metadata folder in github)
      - breast cancer cell lines only
      - ER positive vs. negative metadata

  4) CCLE_RPPA_20181003.csv (found in CCLE site)
      - RPPA data
      - log2 median-centered
      - for more info: https://www.mdanderson.org/content/dam/mdanderson/documents/core-facilities/Functional%20Proteomics%20RPPA%20Core%20Facility/RPPA%20Materials_Methods_2016.pdf

  5) CCLE_RPPA_Ab_info_20181226.csv (found in CCLE site)
      - RPPA antibody/gene data
  
Folders in github that are relevant:
  1) comparison_plots_081420
      - plots that compare integrative vs. RNA-seq only pathways
  2) int_080320
      - plots for integrative pathways
  3) rna_only_080320
      - plots for RNA-seq only pathways
  4) scripts
      - R code (use the most current one: v3.1)
  5) metadata
      - extra metadata sets
  6) presentations
      - presentations on methods and results
  7) lm_FC
      - plots for linear model transformation


