library(tidyverse)

CCLE_RPPA <- read.csv("CCLE_RPPA_20180123.csv", header = FALSE, colClasses = "character")

# Because we are filtering rows whose first column contains "BREAST"
# This header must contain "BREAST" because we want to keep the antibody names
CCLE_RPPA[1,1] <- "BREAST_CELL_LINES_ONLY"

CCLE_RPPA_breast <- CCLE_RPPA %>% filter(grepl("BREAST", V1))


CCLE_RPPA_genes <- read.csv("CCLE_RPPA_Ab_info_20180123.csv", colClasses = "character")
CCLE_RPPA_genes <- CCLE_RPPA_genes[,1:2]
genes <- as.character(CCLE_RPPA_genes$Target_Genes)

keys <- c(genes, "Description")


CCLE_RNA3 <- read.csv("CCLE_RNAseq_rsem_genes_tpm_20180929.csv", header = FALSE)
CCLE_RNA3$V1 <- as.character(CCLE_RNA3$V1)
CCLE_RNA3["1","V1"] <- "Name"
CCLE_RNA3 <- CCLE_RNA3[,c(1,3:1021)]
RNA_genes <- read.csv("RNA_genes_ID.csv", header = FALSE)


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

genes_diff <- setdiff(RPPA_genes, RNA_genes)
genes_diff2 <- setdiff(RNA_genes, RPPA_genes)


#### RPPA
#CCLE_RPPA_genes <- RPPA[,c("BREAST_CELL_LINES_ONLY", "genes")]


ER <- RPPA[c(1,76),]
ER <- as.data.frame(t(ER))
ER <- ER[2:48,]
names(ER) = c("cell_line", "ER_val")

# ER$num <- as.numeric(as.character(ER$ER_val))
ER <- ER %>% mutate(ERsign = case_when(as.numeric(as.character(ER_val)) > 0 ~ "pos",
                                                                       TRUE ~ "neg")) 
# ER Positive and Negative Vectors
posER <- as.character(ER[ER$ERsign == "pos", 1])
negER <- as.character(ER[ER$ERsign == "neg", 1])


names(RPPA) <- lapply(RPPA[1,], as.character)
RPPA <- RPPA[-1,] 
#Create dataframes with new 
RPPApos <- select(RPPA, c("BREAST_CELL_LINES_ONLY", posER))
RPPApos[,2:24] <- lapply(RPPApos[,2:24], as.character)
RPPApos[,2:24] <- lapply(RPPApos[,2:24], as.numeric)
RPPApos <- RPPApos %>% mutate(posRPPAaverage = rowSums(.[2:24])/23)
RPPApos <- RPPApos[,c("BREAST_CELL_LINES_ONLY", "posRPPAaverage")]

RPPAneg <- select(RPPA, c("BREAST_CELL_LINES_ONLY", negER))
RPPAneg[,2:25] <- lapply(RPPAneg[,2:25], as.character)
RPPAneg[,2:25] <- lapply(RPPAneg[,2:25], as.numeric)
RPPAneg <- RPPAneg %>% mutate(negRPPAaverage = rowSums(.[2:25])/24)
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

# Filter those out
names(RNA) <- lapply(RNA[1,], as.character)
RNA <- RNA[-1,] 

#Create dataframes with new 
RNApos <- select(RNA, c("Name","Description", posER))
RNApos[,3:25] <- lapply(RNApos[,3:25], as.character)
RNApos[,3:25] <- lapply(RNApos[,3:25], as.numeric)
RNApos <- RNApos %>% mutate(posRNAaverage = rowSums(.[3:25])/23)
RNApos <- RNApos[,c("Name","Description","posRNAaverage")]

RNAneg <- select(RNA, c("Name","Description",negER))
RNAneg[,3:26] <- lapply(RNAneg[,3:26], as.character)
RNAneg[,3:26] <- lapply(RNAneg[,3:26], as.numeric)
RNAneg <- RNAneg %>% mutate(negRNAaverage = rowSums(.[,3:26])/24)
RNAneg <- RNAneg[,c("Name","Description","negRNAaverage")]

RNAfinal <- left_join(RNApos, RNAneg)



final <- left_join(RPPAfinal, RNAfinal,  by = c("Target_Genes" = "Description"))
final <- final[,c(4,1,5,2,3,6,7)]
names(final) <- c("gene_name", "antibody_name", "RNASeq_name", 
                  "posRPPAaverage", "negRPPAaverage", "posRNAaverage", "negRNAaverage")

final$posRPPAaverage <- 2^(final$posRPPAaverage)
final$negRPPAaverage <- 2^(final$negRPPAaverage)


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
write.csv(final_transform, "final_TPM_transform.csv")

#write.csv(final, "final.csv")
# Manually edited phos_final and nonphos_final
# Remove phos with no non-phos, lowest FC of phos probes with same nonphos, lowest FC of nonphos probes for same gene

phos_final <- read.csv("TPM_phos_final.csv")
nonphos_final <- read.csv("TPM_nonphos_final.csv")
#MAKE FCS, RECOMBO
nonphos_final$RPPAFC <- nonphos_final$posRPPAaverage / nonphos_final$negRPPAaverage
nonphos_final$RNAFC <- nonphos_final$RNA_transform_pos / nonphos_final$RNA_transform_neg

phos_final$antibody_name <- paste(phos_final$antibody_name, phos_final$antibody_name.1, sep = ", ")
phos_final$phosRPPAFC <- phos_final$posRPPAaveragePHOS / phos_final$negRPPAaveragePHOS
phos_final$nonphosRPPAFC <- phos_final$posRPPAaverage / phos_final$negRPPAaverage
phos_final$RPPAFC <- phos_final$phosRPPAFC / phos_final$nonphosRPPAFC
phos_final$RNAFC <- phos_final$RNA_transform_pos / phos_final$RNA_transform_neg


nonphos_final <- nonphos_final[,c(1,2,3,9,10)]
phos_final <- phos_final[,c(1,2,3,15,16)]

all_final <- rbind(nonphos_final, phos_final)
all_final$log2RPPAFC <- log2(all_final$RPPAFC)
all_final$log2RNAFC <- log2(all_final$RNAFC)


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

write.csv(all_final, "final_RPKM_RNA2RPPA_CCLE_analysis_V1.csv")


totRPPA <- sum(all_final$max_RPPAorRNA == "RPPA")
totRNA <- sum(all_final$max_RPPAorRNA == "RNA")     
tot <- totRPPA + totRNA
percRPPA <- totRPPA / tot
percRNA <- totRNA / tot
