library(tidyverse)

## RPPA data input
CCLE_RPPA <- read.csv("CCLE_RPPA_20181003.csv", header = FALSE, colClasses = "character")
#CCLE_RPPA <- read.csv("CCLE_RPPA_20180123.csv", header = FALSE, colClasses = "character")
# This code assumes that the row names are in the first column and not the rownames of the dataframe
# and column names are in the first row

# Because we are filtering rows whose first column contains "BREAST"
# This header must contain "BREAST" because we want to keep the antibody names
CCLE_RPPA[1,1] <- "BREAST_CELL_LINES_ONLY"

CCLE_RPPA_breast <- CCLE_RPPA %>% filter(grepl("BREAST", V1))

# Get genes for RPPA
CCLE_RPPA_genes <- read.csv("CCLE_RPPA_Ab_info_20181226.csv", colClasses = "character")
CCLE_RPPA_genes <- CCLE_RPPA_genes[,1:2]
genes <- as.character(CCLE_RPPA_genes$Target_Genes)

keys <- c(genes, "Description")

## RNASeq data input
CCLE_RNA3 <- read.csv("CCLE_RNAseq_rsem_genes_tpm_20180929.csv", header = FALSE)
CCLE_RNA3$V1 <- as.character(CCLE_RNA3$V1)
CCLE_RNA3["1","V1"] <- "Name"
CCLE_RNA3 <- CCLE_RNA3[,c(1,3:1021)]
#RNA_genes <- read.csv("RNA_genes_ID.csv", header = FALSE)

# Filter the genes
RNA_genes_filter <- RNA_genes %>%  filter(grepl(paste(keys, collapse="|"), V2))

CCLE_RNA5 <- left_join(RNA_genes_filter, CCLE_RNA3)

CCLE_RNA4 <- CCLE_RNA5 %>% filter(grepl(paste(keys, collapse="|"), V2))

CCLE_RNA_gene <- as.data.frame(t(CCLE_RNA4))
CCLE_RNA_gene <- CCLE_RNA_gene %>% filter(grepl("BREAST|Name|Description", V1))


## Find cell lines that dont match (four in RNA that aren't in RPPA)
cell_line_RNA <- as.character(CCLE_RNA_gene$V1)
cell_line_RPPA <- as.character(CCLE_RPPA_breast$V1)
cell_line_RNA <- cell_line_RNA[3:53]
cell_line_RPPA <- sort(cell_line_RPPA[2:48])
cell_line_diff <- setdiff(cell_line_RNA, cell_line_RPPA)
##

CCLE_RNA_filtered <- CCLE_RNA_gene %>% filter(!grepl(paste(cell_line_diff, collapse="|"), V1))

# Transpose both RNA and RPPA so that genes are rows
RNA <- as.data.frame(t(CCLE_RNA_filtered),stringsAsFactors = FALSE)
RPPA <- as.data.frame(t(CCLE_RPPA_breast),stringsAsFactors = FALSE)
CCLE_RPPA_genes <- rbind(CCLE_RPPA_genes, c("BREAST_CELL_LINES_ONLY", "genes"))
RPPA <- left_join(RPPA, CCLE_RPPA_genes, by = c("V1" = "Antibody_Name"))
# Six genes are present in the RPPA data that aren't in the RNA
# Two because they don't have a match, and four because they are in as multiple genes
RNA_genes <- as.character(RNA$V2)
RNA_genes <- RNA_genes[2:1318]
RPPA_genes <- RPPA$Target_Genes
RPPA_genes <- RPPA_genes[2:215]

#### RPPA Analysis

## Extract ER-alpha 
# Create df with ER data for each cell line
ER <- RPPA[c(1,80),]
ER <- as.data.frame(t(ER))
ER <- ER[2:48,]
names(ER) = c("cell_line", "ER_val")
ER <- ER %>% mutate(ERsign = case_when(as.numeric(as.character(ER_val)) > 0 ~ "pos",
                                       TRUE ~ "neg")) 
# ER Positive and Negative Vectors
posER <- as.character(ER[ER$ERsign == "pos", 1])
negER <- as.character(ER[ER$ERsign == "neg", 1])

# Put row names where they should go
names(RPPA) <- lapply(RPPA[1,], as.character)
RPPA <- RPPA[-1,] 
# Extract pos and neg cell lines, average their values in new dfs
RPPApos <- select(RPPA, c("BREAST_CELL_LINES_ONLY", posER))
RPPApos[,2:37] <- lapply(RPPApos[,2:37], as.character)
RPPApos[,2:37] <- lapply(RPPApos[,2:37], as.numeric)
RPPApos <- RPPApos %>% mutate(posRPPAaverage = rowSums(.[2:37])/36)
RPPApos <- RPPApos[,c("BREAST_CELL_LINES_ONLY", "posRPPAaverage")]

RPPAneg <- select(RPPA, c("BREAST_CELL_LINES_ONLY", negER))
RPPAneg[,2:12] <- lapply(RPPAneg[,2:12], as.character)
RPPAneg[,2:12] <- lapply(RPPAneg[,2:12], as.numeric)
RPPAneg <- RPPAneg %>% mutate(negRPPAaverage = rowSums(.[2:12])/11)
RPPAneg <- RPPAneg[,c("BREAST_CELL_LINES_ONLY", "negRPPAaverage")]

RPPAfinal <- left_join(RPPApos, RPPAneg)
RPPAfinal <- left_join(RPPAfinal, CCLE_RPPA_genes, by = c("BREAST_CELL_LINES_ONLY" = "Antibody_Name"))


# Fix/remove the genes_diff
RPPAfinal$Target_Genes <- as.character(RPPAfinal$Target_Genes)
RPPAfinal$Target_Genes <- recode(RPPAfinal$Target_Genes, "ACACA ACACB" = "ACACA")
RPPAfinal$Target_Genes <- recode(RPPAfinal$Target_Genes, "GSK3A GSK3B" = "GSK3B")
RPPAfinal$Target_Genes <- recode(RPPAfinal$Target_Genes, "MAPK1 MAPK3" = "MAPK1")
RPPAfinal$Target_Genes <- recode(RPPAfinal$Target_Genes, "RPS6KA1 RPS6KA2 RPS6KA3" = "RPS6KA3")
RPPAfinal <- subset(RPPAfinal, Target_Genes != "AKT1 AKT2 AKT3" & Target_Genes != "PECAM1")


#### RNA

# Put row names where they should go
names(RNA) <- lapply(RNA[1,], as.character)
RNA <- RNA[-1,] 

# Extract pos and neg cell lines, average their values in new dfs
RNApos <- select(RNA, c("Name","Description", posER))
RNApos[,3:38] <- lapply(RNApos[,3:38], as.character)
RNApos[,3:38] <- lapply(RNApos[,3:38], as.numeric)
RNApos <- RNApos %>% mutate(posRNAaverage = rowSums(.[3:38])/36)
RNApos <- RNApos[,c("Name","Description","posRNAaverage")]

RNAneg <- select(RNA, c("Name","Description",negER))
RNAneg[,3:13] <- lapply(RNAneg[,3:13], as.character)
RNAneg[,3:13] <- lapply(RNAneg[,3:13], as.numeric)
RNAneg <- RNAneg %>% mutate(negRNAaverage = rowSums(.[3:13])/11)
RNAneg <- RNAneg[,c("Name","Description","negRNAaverage")]

RNAfinal <- left_join(RNApos, RNAneg)

# Combine RPPA and RNA data
final <- left_join(RPPAfinal, RNAfinal,  by = c("Target_Genes" = "Description"))
final <- final[,c(4,1,5,2,3,6,7)]
names(final) <- c("gene_name", "antibody_name", "RNASeq_name", 
                  "posRPPAaverage", "negRPPAaverage", "posRNAaverage", "negRNAaverage")

final$posRNAaverage <- log2(final$posRNAaverage + 1)
final$negRNAaverage <- log2(final$negRNAaverage + 1)

# Create model to transform RNA
pos_lm <- lm(data = final, posRNAaverage ~ posRPPAaverage)
pos_b <- summary(pos_lm)$coefficients[1,1]
pos_m <- summary(pos_lm)$coefficients[2,1]
final$RNA_transform_pos <- final$posRPPAaverage*pos_m + pos_b

neg_lm <- lm(data = final, negRNAaverage ~ negRPPAaverage)
neg_b <- summary(neg_lm)$coefficients[1,1]
neg_m <- summary(neg_lm)$coefficients[2,1]
final$RNA_transform_neg <- final$negRPPAaverage*neg_m + neg_b

ggplot(data = final)+
  geom_point(aes(x = posRPPAaverage, y = posRNAaverage)) +
  geom_point(aes(x = posRPPAaverage, y = RNA_transform_pos), color = "blue") +
  geom_abline(slope = pos_m, intercept = pos_b, color = "red") + 
  geom_point(aes(x = negRPPAaverage, y = negRNAaverage), color = "darkgreen") +
  geom_point(aes(x = negRPPAaverage, y = RNA_transform_neg), color = "blue") +
  geom_abline(slope = neg_m, intercept = neg_b, color = "darkred")

final_transform <- final[,c(1:5,8,9)]
write.csv(final_transform, "final_TPM_transform_V2.csv")

##write.csv(final, "final.csv")
# Manually edited phos_final and nonphos_final
# Remove phos with no non-phos, lowest FC of phos probes with same nonphos, lowest FC of nonphos probes for same gene
# Use files as guide

phos_final <- read.csv("TPM_phos_final_V2.csv")
nonphos_final <- read.csv("TPM_nonphos_final_V2.csv")
#MAKE FCS, RECOMBO
nonphos_final$log2RPPAFC <- nonphos_final$posRPPAaverage - nonphos_final$negRPPAaverage
nonphos_final$log2RNAFC <- nonphos_final$RNA_transform_pos - nonphos_final$RNA_transform_neg
nonphos_final$RPPAFC <- 2^(nonphos_final$log2RPPAFC)
nonphos_final$RNAFC <- 2^(nonphos_final$log2RNAFC)

phos_final$antibody_name <- paste(phos_final$antibody_name, phos_final$antibody_name.1, sep = ", ")
phos_final$phosRPPAFC <- phos_final$posRPPAaveragePHOS - phos_final$negRPPAaveragePHOS
phos_final$nonphosRPPAFC <- phos_final$posRPPAaverage - phos_final$negRPPAaverage
phos_final$log2RPPAFC <- phos_final$phosRPPAFC - phos_final$nonphosRPPAFC
phos_final$log2RNAFC <- phos_final$RNA_transform_pos - phos_final$RNA_transform_neg
phos_final$RPPAFC <- 2^(phos_final$log2RPPAFC)
phos_final$RNAFC <- 2^(phos_final$log2RNAFC)

# Edit the dfs
nonphos_final <- nonphos_final[,c(1,2,3,10,11,8,9)]
phos_final <- phos_final[,c(1,2,3,17,18,15,16)]

# Combine them
all_final <- rbind(nonphos_final, phos_final)

#Mark from which dataset the FC is higher
for (ind in 1:nrow(all_final)) {
  if (abs(all_final$log2RPPAFC[ind]) > abs(all_final$log2RNAFC[ind]) ) {
    all_final$max_log_FC[ind] = all_final$log2RPPAFC[ind]
    all_final$max_RPPAorRNA[ind] = "RPPA"
  }
  else {
    all_final$max_log_FC[ind] = all_final$log2RNAFC[ind]
    all_final$max_RPPAorRNA[ind] = "RNA"
  }
}

#write.csv(all_final, "final_TPM_RNA2RPPA_CCLE_analysis_V2.csv", row.names = FALSE)

# Little statistics
totRPPA <- sum(all_final$max_RPPAorRNA == "RPPA")
totRNA <- sum(all_final$max_RPPAorRNA == "RNA")     
tot <- totRPPA + totRNA
percRPPA <- totRPPA / tot
percRNA <- totRNA / tot

# This is to calculate RNA-only genes, same idea
RNA_genes_only <- read.csv("RNA_genes_ID.csv", header = FALSE)

CCLE_RNA6 <- left_join(RNA_genes_only, CCLE_RNA3)
CCLE_RNA_only <- as.data.frame(t(CCLE_RNA6),stringsAsFactors = FALSE)

CCLE_RNA_gene_only <- CCLE_RNA_only %>% filter(grepl("BREAST|Name|Description", V1))

CCLE_RNA_filtered_only <- CCLE_RNA_gene_only %>% filter(!grepl(paste(cell_line_diff, collapse="|"), V1))

RNA_only <- as.data.frame(t(CCLE_RNA_filtered_only),stringsAsFactors = FALSE)

names(RNA_only) <- lapply(RNA_only[1,], as.character)
RNA_only <- RNA_only[-1,] 

#Create dataframes with new 
RNA_onlypos <- select(RNA_only, c("Name","Description", posER))
RNA_onlypos[,3:38] <- lapply(RNA_onlypos[,3:38], as.character)
RNA_onlypos[,3:38] <- lapply(RNA_onlypos[,3:38], as.numeric)
RNA_onlypos <- RNA_onlypos %>% mutate(posRNAaverage = rowSums(.[3:38])/36)
RNA_onlypos <- RNA_onlypos[,c("Name","Description","posRNAaverage")]

RNA_onlyneg <- select(RNA_only, c("Name","Description",negER))
RNA_onlyneg[,3:13] <- lapply(RNA_onlyneg[,3:13], as.character)
RNA_onlyneg[,3:13] <- lapply(RNA_onlyneg[,3:13], as.numeric)
RNA_onlyneg <- RNA_onlyneg %>% mutate(negRNAaverage = rowSums(.[3:13])/11)
RNA_onlyneg <- RNA_onlyneg[,c("Name","Description","negRNAaverage")]

RNA_onlyfinal_TPM <- left_join(RNA_onlypos, RNA_onlyneg)

RNA_onlyfinal_TPM$posRNAaverage <- log2(RNA_onlyfinal_TPM$posRNAaverage + 1)
RNA_onlyfinal_TPM$negRNAaverage <- log2(RNA_onlyfinal_TPM$negRNAaverage + 1)


RNA_onlyfinal_TPM$RPPA_transform_pos <- (RNA_onlyfinal_TPM$posRNAaverage - pos_b) / pos_m 
RNA_onlyfinal_TPM$RPPA_transform_neg <- (RNA_onlyfinal_TPM$negRNAaverage - neg_b) / neg_m

RNA_onlyfinal_TPM$log2RPPAFC <- RNA_onlyfinal_TPM$RPPA_transform_pos - RNA_onlyfinal_TPM$RPPA_transform_neg
RNA_onlyfinal_TPM$log2RNAFC <- RNA_onlyfinal_TPM$posRNAaverage - RNA_onlyfinal_TPM$negRNAaverage
RNA_onlyfinal_TPM$RPPAFC <- 2^(RNA_onlyfinal_TPM$log2RPPAFC)
RNA_onlyfinal_TPM$RNAFC <- 2^(RNA_onlyfinal_TPM$log2RNAFC)

RNA_onlyfinal_TPM <- RNA_onlyfinal_TPM[, c(1:6,9,10,7,8)]

#write.csv(RNA_onlyfinal_TPM, "TPM_RNA_only_transform_V2.csv")



ggplot(data = RNA_onlyfinal_TPM)+
  geom_point(aes(x = RPPA_transform_pos, y = posRNAaverage), color = "blue") +
  geom_abline(slope = pos_m, intercept = pos_b, color = "red") + 
  geom_point(aes(x = RPPA_transform_neg, y = negRNAaverage), color = "darkblue") +
  geom_abline(slope = neg_m, intercept = neg_b, color = "darkred")



ggplot()+
  
  
  geom_abline(slope = pos_m, intercept = pos_b, color = "red") + 
  geom_abline(slope = neg_m, intercept = neg_b, color = "darkred") + 
  
  
  geom_point(data = final,aes(x = posRPPAaverage, y = posRNAaverage)) +
  geom_point(data = final,aes(x = posRPPAaverage, y = RNA_transform_pos), color = "lightblue") +
  geom_point(data = final,aes(x = negRPPAaverage, y = negRNAaverage), color = "darkgreen") +
  geom_point(data = final,aes(x = negRPPAaverage, y = RNA_transform_neg), color = "lightblue")+
  
  geom_point(data = RNA_onlyfinal_TPM,aes(x = RPPA_transform_pos, y = posRNAaverage), color = "pink") +
  geom_point(data = RNA_onlyfinal_TPM,aes(x = RPPA_transform_neg, y = negRNAaverage), color = "grey80")
