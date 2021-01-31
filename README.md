# CCLE_integrative

Integrating CCLE RNA-seq and RPPA data. 

Data using:
  1) raw RNA-seq gene counts*
  2) TPM normalized RNA-seq*
  3) RPKM normalized RNA-seq
  3) normalized RPPA*
* Used these datasets for analysis

Names of datasets used (based on v3.1 script):
1) CCLE_RNAseq_genes_counts_20180929.gct
  - Raw RNA-seq counts
  - used for PCA

2) CCLE_RNAseq_rsem_genes_tpm_20180929.txt
  - RNA-seq counts TPM normalized

3) breast_cell_lines_V2.csv
  - breast cancer cell lines only
  - ER positive vs. negative metadata

4) CCLE_RPPA_20181003.csv
  - RPPA data
  - log2 median-centered
  - for more info: https://www.mdanderson.org/content/dam/mdanderson/documents/core-facilities/Functional%20Proteomics%20RPPA%20Core%20Facility/RPPA%20Materials_Methods_2016.pdf

5) CCLE_RPPA_Ab_info_20181226.csv
  - RPPA antibody/gene data
