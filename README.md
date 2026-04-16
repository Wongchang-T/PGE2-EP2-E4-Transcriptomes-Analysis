# PGE2-EP2-E4-Transcriptomes-Analysis: scRNAsq of human cancer tissues by Wongchang T.
# 1. Packages loading
```{r}
library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(sctransform)
library(dplyr)
library(Rcpp)
library(harmony)
library(dittoSeq)
library(SingleR)
library(reticulate)
library(EnhancedVolcano)
library(ggrepel)
library(ggpubr)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(enrichplot)
library(pheatmap)
library(tibble)
library(patchwork)
library(stringr)
```
# 2. Import dataset
```{r}
# Breast_1-GSE114727 #
BRCA_1_1 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE114727_BREAST/GSM3148575_BreastTumor")
BRCA_1_2 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE114727_BREAST/GSM3148576_BreastTumor")
BRCA_1_3 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE114727_BREAST/GSM3148577_BreastTumor")
BRCA_1_4 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE114727_BREAST/GSM3148578_BreastTumor")
BRCA_1_5 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE114727_BREAST/GSM3148579_BreastTumor")

# Breast_2-GSE242271 #
BRCA_2_1 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE242271/scBRCA1_CD45")
BRCA_2_2 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE242271/scBRCA3_CD45")
BRCA_2_3 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE242271/scBRCA4_CD45")
BRCA_2_4 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE242271/scBRCA5_CD45")
BRCA_2_5 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE242271/scBRCA6_CD45")

# Colorectal_1-GSE242271 #
CRC_1_1 <-  Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE242271/scCRC1_CD45")
CRC_1_2 <-  Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE242271/scCRC2_CD45")
CRC_1_3 <-  Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE242271/scCRC3_CD45")
CRC_1_4 <-  Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE242271/scCRC4_CD45")
CRC_1_5 <-  Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE242271/scCRC5_CD45")

# Colorectal_2-GSE188711 #
CRC_2_1 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE188711_COLORECTAL/GSM5688706_CRC")
CRC_2_2 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE188711_COLORECTAL/GSM5688707_CRC")
CRC_2_3 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE188711_COLORECTAL/GSM5688708_CRC")
CRC_2_4 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE188711_COLORECTAL/GSM5688709_CRC")
CRC_2_5 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE188711_COLORECTAL/GSM5688710_CRC")
CRC_2_6 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE188711_COLORECTAL/GSM5688711_CRC")

# Colorectal_3-GSE161277 #
CRC_3_1 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE161277_COLON/GSM4904234_Carcinoma")
CRC_3_2 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE161277_COLON/GSM4904236_Carcinoma")
CRC_3_3 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE161277_COLON/GSM4904239_Carcinoma")
CRC_3_4 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE161277_COLON/GSM4904245_Carcinoma")

# Colorectal_4-GSE280318
CRC_4_1 <- Read10X_h5("D:/HumanTumor_scRNAseq/GSE280318_COLON/GSM8594533_P1L1_P2CRC_BC1_count_filtered_feature_bc_matrix.h5")
CRC_4_2 <- Read10X_h5("D:/HumanTumor_scRNAseq/GSE280318_COLON/GSM8594533_P1L1_P4CRC_BC2_count_filtered_feature_bc_matrix.h5")
CRC_4_3 <- Read10X_h5("D:/HumanTumor_scRNAseq/GSE280318_COLON/GSM8594534_P1L2_P2CRC_BC1_count_filtered_feature_bc_matrix.h5")
CRC_4_4 <- Read10X_h5("D:/HumanTumor_scRNAseq/GSE280318_COLON/GSM8594534_P1L2_P4CRC_BC2_count_filtered_feature_bc_matrix.h5")
CRC_4_5 <- Read10X_h5("D:/HumanTumor_scRNAseq/GSE280318_COLON/GSM8594535_P1L3_P2CRC_BC1_count_filtered_feature_bc_matrix.h5")
CRC_4_6 <- Read10X_h5("D:/HumanTumor_scRNAseq/GSE280318_COLON/GSM8594535_P1L3_P4CRC_BC2_count_filtered_feature_bc_matrix.h5")
CRC_4_7 <- Read10X_h5("D:/HumanTumor_scRNAseq/GSE280318_COLON/GSM8594536_P1L4_P2CRC_BC1_count_filtered_feature_bc_matrix.h5")
CRC_4_8 <- Read10X_h5("D:/HumanTumor_scRNAseq/GSE280318_COLON/GSM8594536_P1L4_P4CRC_BC2_count_filtered_feature_bc_matrix.h5")

# Gastric_1-GSE163558 #
STAD_1_1 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE163558_GASTRIC/GSM5004180_PT1")
STAD_1_2 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE163558_GASTRIC/GSM5004181_PT2")
STAD_1_3 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE163558_GASTRIC/GSM5004182_PT3")

# NSCLC_1-GSE183219 #
NSCLC_1_1 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE183219_NSCLC/GSM5553494")
NSCLC_1_2 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE183219_NSCLC/GSM5553496")
NSCLC_1_3 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE183219_NSCLC/GSM5553498")
NSCLC_1_4 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE183219_NSCLC/GSM5553500")
NSCLC_1_5 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE183219_NSCLC/GSM5553502")
NSCLC_1_6 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE183219_NSCLC/GSM5553504")
NSCLC_1_7 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE183219_NSCLC/GSM5553506")
NSCLC_1_8 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE183219_NSCLC/GSM5553507")
NSCLC_1_9 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE183219_NSCLC/GSM5553508")
NSCLC_1_10 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE183219_NSCLC/GSM5553510")
NSCLC_1_11 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE183219_NSCLC/GSM5553512")
NSCLC_1_12 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE183219_NSCLC/GSM5553513")
NSCLC_1_13 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE183219_NSCLC/GSM5553516")
NSCLC_1_14 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE183219_NSCLC/GSM5553518")
NSCLC_1_15 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE183219_NSCLC/GSM5553519")

#NSCLC_2-GSE300685
NSCLC_2_1 <- Read10X_h5("D:/HumanTumor_scRNAseq/GSE300685_NSCLC/GSM9066344_LC104_filtered_feature_bc_matrix.h5")
NSCLC_2_2 <- Read10X_h5("D:/HumanTumor_scRNAseq/GSE300685_NSCLC/GSM9066345_LC115_filtered_feature_bc_matrix.h5")
NSCLC_2_3 <- Read10X_h5("D:/HumanTumor_scRNAseq/GSE300685_NSCLC/GSM9066346_LC221_filtered_feature_bc_matrix.h5")
NSCLC_2_4 <- Read10X_h5("D:/HumanTumor_scRNAseq/GSE300685_NSCLC/GSM9066347_LC38_filtered_feature_bc_matrix.h5")
NSCLC_2_5 <- Read10X_h5("D:/HumanTumor_scRNAseq/GSE300685_NSCLC/GSM9066348_LC52_filtered_feature_bc_matrix.h5")
NSCLC_2_6 <- Read10X_h5("D:/HumanTumor_scRNAseq/GSE300685_NSCLC/GSM9066349_LC57_filtered_feature_bc_matrix.h5")
NSCLC_2_7 <- Read10X_h5("D:/HumanTumor_scRNAseq/GSE300685_NSCLC/GSM9066350_LC63_filtered_feature_bc_matrix.h5")
NSCLC_2_8 <- Read10X_h5("D:/HumanTumor_scRNAseq/GSE300685_NSCLC/GSM9066351_LC71_filtered_feature_bc_matrix.h5")

# Pancrease_1-GSE197177 # 
PAAD_1_1 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE197177_PDAC/GSM5910784")
PAAD_1_2 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE197177_PDAC/GSM5910785")
PAAD_1_3 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE197177_PDAC/GSM5910786")
PAAD_1_4 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE197177_PDAC/GSM5910787")
PAAD_1_5 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE197177_PDAC/GSM5910788")
PAAD_1_6 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE197177_PDAC/GSM5910789")
PAAD_1_7 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE197177_PDAC/GSM5910790")
PAAD_1_8 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE197177_PDAC/GSM5910791")

# Liver_1-GSE235057 #
LIHC_1_1 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE235057_LIVER/GSM7494107_HCC01T")
LIHC_1_2 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE235057_LIVER/GSM7494109_HCC02T")
LIHC_1_3 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE235057_LIVER/GSM7494111_HCC03T")
LIHC_1_4 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE235057_LIVER/GSM7494113_HCC04T")
LIHC_1_5 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE235057_LIVER/GSM7494115_HCC05T")
LIHC_1_6 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE235057_LIVER/GSM7494117_HCC06T")
LIHC_1_7 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE235057_LIVER/GSM7494119_HCC07T")
LIHC_1_8 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE235057_LIVER/GSM7494121_HCC08T")
LIHC_1_9 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE235057_LIVER/GSM7494123_HCC09T")
LIHC_1_10 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE235057_LIVER/GSM7494125_HCC10T")
LIHCM_1_1 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE235057_LIVER/GSM7494127_LM02T1")
LIHCM_1_2 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE235057_LIVER/GSM7494126_LM01T")
LIHCM_1_3 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE235057_LIVER/GSM7494128_LM02T2")
LIHCM_1_4 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE235057_LIVER/GSM7494129_LM03T1")
LIHCM_1_5 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE235057_LIVER/GSM7494130_LM03T2")
LIHCM_1_6 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE235057_LIVER/GSM7494131_LM04T")
LIHCM_1_7 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE235057_LIVER/GSM7494132_LM05T")
LIHCM_1_8 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE235057_LIVER/GSM7494133_LM06T")
LIHCM_1_9 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE235057_LIVER/GSM7494134_LM07T")

# CholangioCA-GSE213452
CHOL_1_1 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE213452_CCA/GSM6586313_LFY_pCCA_C")
CHOL_1_2 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE213452_CCA/GSM6586315_SKF_pCCA_C")
CHOL_1_3 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE213452_CCA/GSM6586317_SLB_pCCA_C")
CHOL_1_4 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE213452_CCA/GSM6586318_WZG_dCCA_C")
CHOL_1_5 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE213452_CCA/GSM6586320_YP_dCCA_C")
CHOL_1_6 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE213452_CCA/GSM6586322_ZhangLZ_dCCA_C")
CHOL_1_7 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE213452_CCA/GSM6586323_ZhuLZ_dCCA_C")

# Ovarian_1-GSE114727 #
OVCA_1_1 <-  Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE242271/scOVCA1_CD45")
OVCA_1_2 <-  Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE242271/scOVCA2_CD45")
OVCA_1_3 <-  Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE242271/scOVCA3_CD45")
OVCA_1_4 <-  Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE242271/scOVCA4_CD45")
OVCA_1_5 <-  Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE242271/scOVCA5_CD45")

# Bladder_1-GSE277524 #
BLCA_1_1 <-  Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE277524_BLADDER/GSM8523962_P01_R1")
BLCA_1_2 <-  Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE277524_BLADDER/GSM8523963_P01_R2")
BLCA_1_3 <-  Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE277524_BLADDER/GSM8523964_P02_R1")
BLCA_1_4 <-  Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE277524_BLADDER/GSM8523965_P02_R2")
BLCA_1_5 <-  Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE277524_BLADDER/GSM8523966_P03_R1")
BLCA_1_7 <-  Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE277524_BLADDER/GSM8523968_P04_R1")
BLCA_1_8 <-  Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE277524_BLADDER/GSM8523969_P04_R2")
BLCA_1_9 <-  Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE277524_BLADDER/GSM8523970_P05_R1")
BLCA_1_10 <-  Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE277524_BLADDER/GSM8523971_P05_R2")

# CerviX1-GSE224327 #
CESC_1_1 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE224327_CERVIX/GSM7019487_PT1")
CESC_1_2 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE224327_CERVIX/GSM7019488_PT2")
CESC_1_3 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE224327_CERVIX/GSM7019489_PT3")
CESC_1_4 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE224327_CERVIX/GSM7019490_PT4")
CESC_1_5 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE224327_CERVIX/GSM7019491_PT5")
CESC_1_6 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE224327_CERVIX/GSM7019492_PT6")

# Head & Neck-GSE301720 #
HNSC_1_1 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE301720_HEAD&NECK/GSM9087747_CK17-5")
HNSC_1_2 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE301720_HEAD&NECK/GSM9087748_CK17-7")
HNSC_1_3 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE301720_HEAD&NECK/GSM9087749_CK17-12-E6")
HNSC_1_4 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE301720_HEAD&NECK/GSM9087750_CK17-19")
HNSC_1_5 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE301720_HEAD&NECK/GSM9087751_CK17-25")
HNSC_1_6 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE301720_HEAD&NECK/GSM9087752_CK17-129-4")
HNSC_1_7 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE301720_HEAD&NECK/GSM9087753_CK17-159-2")
HNSC_1_8 <- Read10X(data.dir = "D:/HumanTumor_scRNAseq/GSE301720_HEAD&NECK/GSM9087754_CK17-208")
```
# 3. Create suerat object
```{r}
# Get a list of all variable names in your environment that look like your data
count_matrices <- ls(pattern = "^(BRCA|CRC|STAD|NSCLC|PAAD|LIHC|CHOL|LIHCM|OVCA|BLCA|CESC|HNSC)_")

# Loop through the names, create the Seurat object, and overwrite the original variable
for (name in count_matrices) {
  # Get the actual matrix data using the name
  counts <- get(name)
  
  # Create the Seurat object
  # project = name ensures the 'orig.ident' column matches your variable name
  obj <- CreateSeuratObject(counts = counts, project = name, min.cells = 3, min.features = 50)
  
  # Assign the Seurat object back to the original name
  assign(name, obj)
}

#  Clean up temporary variables
rm(counts, obj, name, count_matrices)
```

# 4. Merge seurat object
```{r}
# Collect all the newly created Seurat objects into a list
seurat_obj_names <- ls(pattern = "^(BRCA|CRC|STAD|NSCLC|PAAD|LIHC|CHOL|LIHCM|OVCA|BLCA|CESC|HNSC)_")

# Convert the names into the actual objects
list_of_objects <- lapply(seurat_obj_names, get)

# Merge them all into one object
# We use the first object as the 'x' and the rest of the list as 'y'
Human_cancers <- merge(
  x = list_of_objects[[1]], 
  y = list_of_objects[-1], 
  add.cell.ids = seurat_obj_names, 
  project = "Human_Cancers_Combined"
)

# Verify the result
print(Human_cancers)
table(Human_cancers$orig.ident) # Check cell counts per sample
```

# 5. Quality control check
```{r}
Human_cancers[["percent.mt"]] <- PercentageFeatureSet(Human_cancers, pattern = "^MT")
VlnPlot(Human_cancers, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol=3)
VlnPlot(Human_cancers, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol=3, pt.size = 0)
```

```{r}
plot(FeatureScatter(Human_cancers, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster= FALSE)) + NoLegend()
plot(FeatureScatter(Human_cancers, feature1 = "nFeature_RNA", feature2 = "percent.mt", raster= FALSE)) + NoLegend()
```

```{r}
Human_cancers <- subset(Human_cancers,subset = nFeature_RNA > 50 & nFeature_RNA < 10000 & percent.mt < 15)
head(Human_cancers)
```

# 6. Normalizing & integration
```{r}
Human_cancers <- NormalizeData(Human_cancers)
Human_cancers <- FindVariableFeatures(Human_cancers)
Human_cancers <- ScaleData(Human_cancers)
Human_cancers <- RunPCA(Human_cancers, features = VariableFeatures(Human_cancers))
Human_cancers <- RunHarmony(Human_cancers, group.by.vars = "orig.ident")
```
# 7. Unspervised clustering
```{r}
Human_cancers <- FindNeighbors(Human_cancers, dims = 1:20)
Human_cancers <- FindClusters(Human_cancers, resolution = 0.1, algorithm = 4)
Human_cancers <- RunUMAP(Human_cancers, dims = 1:20)
DimPlot(Human_cancers, reduction = "umap", label = T)
DimPlot(Human_cancers, reduction = "umap", group.by = "orig.ident")
```
# 8. Identification of immune cell populations
```{r}
DimPlot(Human_cancers, reduction = "umap", label = T)
FeaturePlot(Human_cancers, features = c("PTPRC", "CD3E", "CD14", "C1QA", "CSF1R", "CSF3R", "VCAN", "CD79A"), ncol = 4)
FeaturePlot(Human_cancers, features = c("PTPRC", "CCR7", "CLEC9A", "CLEC10A", "CCR9", "CLEC4C"), ncol = 3)
FeaturePlot(Human_cancers, features = c("PTPRC", "EPCAM", "KRT8", "KRT18", "COL1A1", "PECAM1"), ncol = 3)
VlnPlot(Human_cancers, features = "PTPRC")
```
# 9. Immune subset and unsupervised clustering
```{r}
# Do subset
Immune_cells <- subset(Human_cancers, idents = c("1","11","14", "17","18","19","20"), invert =TRUE)
```
```{r}
# Re-PCA and harmonisation
Immune_cells <- RunPCA(Immune_cells, features = VariableFeatures(object = Immune_cells))
Immune_cells <- RunHarmony(Immune_cells, group.by.vars = "orig.ident")
```
```{r}
# Clustering and UAMP plot
Immune_cells <- FindNeighbors(Immune_cells, reduction = "harmony", dims = 1:20)
Immune_cells <- FindClusters(Immune_cells, resolution = 0.2, algorithm = 4) 
Immune_cells <- RunUMAP(Immune_cells, reduction = "harmony", dims = 1:20)
DimPlot(Immune_cells, reduction = "umap", label = TRUE) + NoLegend()
DimPlot(Immune_cells, reduction = "umap", group.by = "orig.ident", label = TRUE) + NoLegend()
```
# 10. Automated cell type annotation by SingleR
```{r}
# Collapse the layers first
Immune_cells <- JoinLayers(Immune_cells)

# Load the reference database
Reference_database <- celldex::HumanPrimaryCellAtlasData()

# Run SingleR annotation
Result_annotation <- SingleR(
  test = as.SingleCellExperiment(Immune_cells), 
  ref = Reference_database, 
  labels = Reference_database$label.main, 
  assay.type.test = 1
)

# Store the labels back into the Seurat object metadata
Immune_cells$SingleR.Label <- Result_annotation$labels
```
```{r}
#Remove non-immune cell by singleR
Idents(Immune_cells) <- "SingleR.Label"
Immune_cells <- subset(Immune_cells, idents = c("Astrocyte",
                                                "Chondrocytes",
                                                "Epithelial_cells", 
                                                "Embryonic_stem_cells",
                                                "Endothelial_cells", 
                                                "Erythroblast", 
                                                "Fibroblasts",
                                                "Gametocytes",
                                                "Hepatocytes",
                                                "Keratinocytes",
                                                "Neurons",
                                                "Neuroepithelial_cells",
                                                "Osteoblasts",
                                                "Platelets",
                                                "Smooth_muscle_cells",
                                                "Tissue_stem_cells"), invert = T) 
DimPlot(Immune_cells, group.by = "SingleR.Label", label = T)
```
# 11. Validated cell type annotation using FindAllmaerker function for DEG
```{r}
Idents(Immune_cells) <- "seurat_clusters"
Immune_Marker <- FindAllMarkers(Immune_cells, logic.threshold = 0.25, min.pct = 0.25, only.pos = T, assay = "RNA")
Immune_Marker %>%
  group_by(cluster) %>%
  top_n(n =20, wt = avg_log2FC)
```
```{r}
Idents(Immune_cells) <- "seurat_clusters"
Immune_cells <- RenameIdents(Immune_cells, 
                            '1' = "T-cells", 
                            '2' = "TAMs",
                            '3' = "T-cells",
                            '4' = "Plasma_cells",
                            '5' = "B-cells",
                            '6' = "TANs",
                            '7' = "Mast-cells")
Immune_cells$Cell_types <- Idents(Immune_cells)
head(Immune_cells)
DimPlot(Immune_cells, reduction = "umap", group.by = "seurat_clusters", label = T)
DimPlot(Immune_cells, reduction = "umap", group.by = "Cell_types", label = T)
```
```{r}
Current_levels <- levels(factor(Immune_cells$Cell_types))
New_order <- c("T-cells", "B-cells", "TAMs", "TANs","Plasma-cells", "Mast-cells",
               setdiff(Current_levels, c("T-cells", "B-cells", "TAMs", "TANs","Plasma-cells", "Mast-cells")))
Immune_cells$Cell_types <- factor(Immune_cells$Cell_types, levels = New_order)
DotPlot(Immune_cells, features = c("CD2",  "CD3E", "CD3D","IL7R", "TCF7", 
                                   "CD8A", "NKG7", "CCL5","GZMA", "GZMB", "GZMK", "PRF1", "GNLY",
                                   "CD4", "CTLA4", "IL2RA", "IL2RG","FOXP3", "TIGIT", "ICOS",
                                   "MS4A1", "CD79A", "BANK1", "HLA-DQA1", "HLA-DRA",
                                   "CD14", "CD68", "CD163","C1QA", "C1QB", "C1QC", "CSF1R", "APOE", "MRC1",
                                   "CSF3R", "CXCL8", "S100A8", "S100A9", "TREM1", "SOD2", "G0S2",
                                   "CPA", "MS4A2", "KIT", "GATA2", "IL1RL1", "HPGDS",
                                   "MKI67", "TOP2A", "STMN1", "TYMS", "CDK1",
                                   "MZB1", "JCHAIN", "DERL3", "IGKC", "IGHAI", "IGHG1"), 
        dot.scale = 10, group.by = "Cell_types") + 
        scale_colour_gradient2(low="darkblue", "red", high="yellow") + 
        theme(axis.title= element_blank(), 
              axis.text.x = element_text(size = 10, angle = 90,color="black", hjust = 1, vjust = 0.3),                                 axis.text.y = element_text(size = 15, color ="black"),legend.position = "bottom")
```
# 12. TAMs subset, unsuperivsed clustering and count abundance in various cancers
```{r}
Idents(Immune_cells) <- "Cell_types"
TAMs <- subset(Immune_cells, idents = "TAMs")
head(TAMs)
```
```{r}
# Run the clustering pipeline using the new object name
TAMs <- FindNeighbors(TAMs, dims = 1:20, reduction = "harmony")
TAMs <- FindClusters(TAMs, resolution = 0.5, reduction = "harmony", algorithm = 4)
TAMs <- RunUMAP(TAMs, dims = 1:20, reduction = "harmony", min.dist = 0.5)

# Visualize with the new object
DimPlot(TAMs, reduction = "umap", label = TRUE, group.by = "seurat_clusters")
DimPlot(TAMs, reduction = "umap", group.by = "Cancer_types")
```
```{r}
Idents(TAMs) <- "seurat_clusters"
TAM_Marker <- FindAllMarkers(TAMs , logic.threshold = 0.25, min.pct = 0.25, only.pos = T, assay = "RNA")
TAM_Marker %>%
  group_by(cluster) %>%
  top_n(n =20, wt = avg_log2FC)
```
```{r}
# 1. Create a data frame from your table
cancer_counts <- as.data.frame(table(TAMs@meta.data$Specific_organ))
colnames(cancer_counts) <- c("Specific_organ", "Count")

# 2. Plotting
ggplot(cancer_counts, aes(x = reorder(Specific_organ, Count), y = Count, fill = Specific_organ)) +
  geom_bar(stat = "identity") +
  coord_flip() +  
  theme_minimal() +
  NoLegend() +
  labs(
    title = "TAMs Distribution Across Solid Cancers",
    x = "Specific_organ",
    y = "Number of TAMs"
  ) +
  scale_fill_viridis_d(option = "turbo") 
  theme(legend.position = "none")       
```
# 13. Investigation of PTGER1-4 expressed by TAMs
```{r}
VlnPlot(TAMs, features = c("PTGER1", "PTGER2", "PTGER3", "PTGER4"), group.by = "Cell_types", ncol = 4)
VlnPlot(TAMs, 
        features = c("PTGER1", "PTGER2", "PTGER3", "PTGER4"), 
        group.by = "Specific_organ", 
        ncol = 4) & # Use & instead of +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
```
# 14. Stratification of TAMs based on PTGER2 and PTGER4 expression level
```{r}
#PTGER4 Level of Expression#
PTGER4_EXPR <- FetchData(TAMs, vars = "PTGER4")
TAMs$PTGER4_EXPR <- PTGER4_EXPR$PTGER4
TAMs$PTGER4_level <- cut(
  TAMs$PTGER4_EXPR ,
  breaks = c(-Inf, 0, 1, 2, Inf),
  labels = c("None", "Low", "Moderate", "High")
)
table(TAMs$PTGER4_level)
head(TAMs)

#PTGER2 Level of Expression#
PTGER2_EXPR <- FetchData(TAMs, vars = "PTGER2")
TAMs$PTGER2_EXPR <- PTGER2_EXPR$PTGER2
TAMs$PTGER2_level <- cut(
  TAMs$PTGER2_EXPR ,
  breaks = c(-Inf, 0, 1, 2, Inf),
  labels = c("None", "Low", "Moderate", "High")
)
table(TAMs$PTGER2_level)
head(TAMs)
```
# 15 Perform DEG analysis
# 15.1 Comparison of PTGER4 high versus low
```{r}
Idents(TAMs) <- "PTGER4_level"
TAMs_PTGER4_Hi_LowEQ <- FindMarkers(TAMs, logfc.threshold = 0,
                                  min.pct = 0.1, ident.1 = High_PTGER4, ident.2 = Low_PTGER4 , 
                                  test.use = "MAST", 
                                  latent.vars = c("orig.ident", "nCount_RNA"),
                                  max.cells.per.ident = 2218)
TAMs_PTGER4_Hi_LowEQ %>% 
  slice_max(n = Inf, order_by = avg_log2FC)
```
# 15.2 Comparison of PTGER2 high versus low
```{r}
# Comparison of PTGER2 high versus low
Idents(TAMs) <- "PTGER2_level"
TAMs_PTGER4_Hi_LowEQ <- FindMarkers(TAMs, logfc.threshold = 0,
                                  min.pct = 0.1, ident.1 = High_PTGER2, ident.2 = Low_PTGER2 , 
                                  test.use = "MAST", 
                                  latent.vars = c("orig.ident", "nCount_RNA"),
                                  max.cells.per.ident = 210)
TAMs_PTGER2_Hi_LowEQ %>% 
  slice_max(n = Inf, order_by = avg_log2FC)
```
# 16. Perform GSEA of PTGER4 high vs. low
# 16.1 GO-BP 
```{r}
# Convert rownames → column
gene_listEP4HiLowEQ <- TAMs_PTGER4_Hi_LowEQ %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::select(gene, avg_log2FC) %>%
  distinct(gene, .keep_all = TRUE) %>%
  filter(!is.na(avg_log2FC)) %>%
  arrange(desc(avg_log2FC))

# Rank the gene
gene_rankEP4HiLowEQ <- gene_listEP4HiLowEQ$avg_log2FC
names(gene_rankEP4HiLowEQ) <- gene_listEP4HiLowEQ$gene
gene_rankEP4HiLowEQ <- sort(gene_rankEP4HiLowEQ, decreasing = TRUE)
```

```{r}
# GO-BP of PTGER4 hihg versus low
gsea_GO_BP_EP4EQ <- gseGO(
  geneList     = gene_rankEP4HiLowEQ,
  OrgDb        = org.Hs.eg.db,    
  ont          = "BP",              
  keyType      = "SYMBOL",          
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE,
  eps          = 0                 
)

# View top results
head(gsea_GO_BP_EP4EQ, n = Inf)
```
# 16.2 KEGG
```{r}
# Convert gene symbol to ENTREZID
gene_dfEQ <- bitr(
  names(gene_rankEP4HiLowEQ),
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)
gene_rank_keggEP4EQ <- gene_rankEP4HiLowEQ[gene_dfEQ$SYMBOL]
names(gene_rank_keggEP4EQ) <- gene_dfEQ$ENTREZID
gene_rank_keggEP4EQ <- sort(gene_rank_keggEP4EQ, decreasing = TRUE)
```
```{r}
# KEGG PTGER4 high versus low
gsea_kegg_EP4EQ <- gseKEGG(
  geneList     = gene_rank_keggEP4EQ,
  organism     = "hsa",
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)
head(gsea_kegg_EP4EQ, Inf)
```
# 16.3 Explore efferocytosis signature gene expression
```{r}
# 1. Flatten the list into a character vector for filtering
Efferocytosis_signature <- list(c(
  # 1. Receptors & Recognition (The "Eat-me" sensors)
  "MERTK", "AXL", "TYRO3", "CD36", "SCARF1", "STAB2", "LRP1", "MRC1", "TIMD4", 
  # 2. Bridging Molecules (Opsonins)
  "GAS6", "PROS1", "MFGE8", "C1QA", "C1QB", "C1QC", "THBS1", 
  # 3. Signaling & Cytoskeletal Remodeling (The "Engulfment" machinery)
  "RAC1", "RHOA", "CDC42", "ELMO1", "DOCK1", "GULP1", "MYO1G", "PLD1", "PTK2",
  # 4. Digestion & Phagosome Maturation
  "RAB5A", "RAB7A", "LAMP1", "LAMP2", "GRN", "CTSD",
  # 5. Anti-inflammatory Response (The "Post-clearance" resolution)
  "TGFB1", "IL10", "ARG1", "NR4A1", "PPARG"
))
efferocytosis_genes <- unlist(Efferocytosis_signature)
# 2. Filter the marker table (Effer_Markers) to include ONLY these signature genes
effer_sig_data_PTGER4EQ <- TAMs_PTGER4_Hi_LowEQ %>%
  filter(rownames(.) %in% efferocytosis_genes)

# 3. Create the plot
EnhancedVolcano(effer_sig_data_PTGER4EQ,
    lab = rownames(effer_sig_data_PTGER4EQ),
    x = 'avg_log2FC',
    y = 'p_val_adj',
        selectLab = rownames(effer_sig_data_PTGER4EQ), # Selects every gene in the table
    max.overlaps = Inf,                          # Allows labels to find a spot even if crowded
    drawConnectors = TRUE,                       # Essential to link labels to non-sig dots    
    xlab = bquote(~Log[2]~ 'Fold Change'),
    ylab = bquote(~-Log[10]~ 'Adjusted P-value'),
    pCutoff = 0.05,
    FCcutoff = 0.25,
    pointSize = 4.0,           
    labSize = 3,                               # Slightly smaller to fit all 36 genes
    labFace = 'bold', 
    boxedLabels = TRUE,                         # Removing boxes helps save space
    widthConnectors = 0.5,
    colConnectors = 'grey30',
    title = 'Efferocytosis Signature: PTGER4-High vs Low (Equlized)',
    subtitle = 'All Signature Genes Labeled (Significant & Non-Significant)',
    legendPosition = 'bottom', 
    gridlines.major = FALSE, 
    gridlines.minor = FALSE,
    border = 'full'
) + 
theme_classic() + 
theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
```
# 17. Perform GSEA of PTGER4 high vs. low
# 17.1 GO-BP
```{r}
# Convert rownames → column
gene_listEP2HiLowEQ <- TAMs_PTGER2_Hi_LowEQ %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::select(gene, avg_log2FC) %>%
  distinct(gene, .keep_all = TRUE) %>%
  filter(!is.na(avg_log2FC)) %>%
  arrange(desc(avg_log2FC))
# Rank the gene
gene_rankEP2HiLowEQ <- gene_listEP2HiLowEQ$avg_log2FC
names(gene_rankEP2HiLowEQ) <- gene_listEP2HiLowEQ$gene
gene_rankEP2HiLowEQ <- sort(gene_rankEP2HiLowEQ, decreasing = TRUE)
```
```{r}
# GO-BP of PTGER2 hihg versus low
gsea_GO_BP_EP2EQ <- gseGO(
  geneList     = gene_rankEP2HiLowEQ,
  OrgDb        = org.Hs.eg.db,    
  ont          = "BP",             
  keyType      = "SYMBOL",          
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE,
  eps          = 0                 
)

# View top results
head(gsea_GO_BP_EP2EQ, n = Inf)
```
# 17.2 KEGG
```{r}
# Convert gene symbol to ENTREZID
gene_df2EQ <- bitr(
  names(gene_rankEP2HiLowEQ),
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)
gene_rank_keggEP2EQ <- gene_rankEP2HiLowEQ[gene_df2EQ$SYMBOL]
names(gene_rank_keggEP2EQ) <- gene_df2EQ$ENTREZID
gene_rank_keggEP2EQ <- sort(gene_rank_keggEP2EQ, decreasing = TRUE)
```
```{r}
# KEGG PTGER4 high versus low
gsea_kegg_EP2EQ <- gseKEGG(
  geneList     = gene_rank_keggEP2EQ,
  organism     = "hsa",
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  eps           = 0,
  verbose      = FALSE
)
# View top results
head(gsea_kegg_EP2EQ, n = Inf)
```
# 17.3 Explore efferocytosis signature gene expression
```{r}
# 1. Flatten the list into a character vector for filtering
Efferocytosis_signature <- list(c(
  # 1. Receptors & Recognition (The "Eat-me" sensors)
  "MERTK", "AXL", "TYRO3", "CD36", "SCARF1", "STAB2", "LRP1", "MRC1", "TIMD4", 
  # 2. Bridging Molecules (Opsonins)
  "GAS6", "PROS1", "MFGE8", "C1QA", "C1QB", "C1QC", "THBS1", 
  # 3. Signaling & Cytoskeletal Remodeling (The "Engulfment" machinery)
  "RAC1", "RHOA", "CDC42", "ELMO1", "DOCK1", "GULP1", "MYO1G", "PLD1", "PTK2",
  # 4. Digestion & Phagosome Maturation
  "RAB5A", "RAB7A", "LAMP1", "LAMP2", "GRN", "CTSD",
  # 5. Anti-inflammatory Response (The "Post-clearance" resolution)
  "TGFB1", "IL10", "ARG1", "NR4A1", "PPARG"
))

efferocytosis_genes <- unlist(Efferocytosis_signature)

# 2. Filter the marker table (Effer_Markers) to include ONLY these signature genes
effer_sig_data_PTGER2EQ <- TAMs_PTGER2_Hi_LowEQ %>%
  filter(rownames(.) %in% efferocytosis_genes)

# 3. Create the plot
EnhancedVolcano(effer_sig_data_PTGER2EQ,
    lab = rownames(effer_sig_data_PTGER2EQ),
    x = 'avg_log2FC',
    y = 'p_val_adj',
    selectLab = rownames(effer_sig_data_PTGER2EQ), # Selects every gene in the table
    max.overlaps = Inf,                          # Allows labels to find a spot even if crowded
    drawConnectors = TRUE,                       # Essential to link labels to non-sig dots  
    xlab = bquote(~Log[2]~ 'Fold Change'),
    ylab = bquote(~-Log[10]~ 'Adjusted P-value'),
    pCutoff = 0.05,
    FCcutoff = 0.25,
    pointSize = 4.0,           
    labSize = 3.5,                               # Slightly smaller to fit all 36 genes
    labFace = 'bold', 
    boxedLabels = TRUE,                         # Removing boxes helps save space
    widthConnectors = 0.5,
    colConnectors = 'grey30',
    title = 'Efferocytosis Signature: PTGER2-High vs Low (Equlized)',
    subtitle = 'All Signature Genes Labeled (Significant & Non-Significant)',
    legendPosition = 'bottom', 
    gridlines.major = FALSE, 
    gridlines.minor = FALSE,
    border = 'full'
) + 
theme_classic() + 
theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
```
# 18. Calculation of module score
# 18.1 Efferocytosis level
```{r}
Efferocytosis_signature <- list(c(
  # 1. Receptors & Recognition (The "Eat-me" sensors)
  "MERTK", "AXL", "TYRO3", "CD36", "SCARF1", "STAB2", "LRP1", "MRC1", "TIMD4", 
  # 2. Bridging Molecules (Opsonins)
  "GAS6", "PROS1", "MFGE8", "C1QA", "C1QB", "C1QC", "THBS1", 
  # 3. Signaling & Cytoskeletal Remodeling (The "Engulfment" machinery)
  "RAC1", "RHOA", "CDC42", "ELMO1", "DOCK1", "GULP1", "MYO1G", "PLD1", "PTK2",
  # 4. Digestion & Phagosome Maturation
  "RAB5A", "RAB7A", "LAMP1", "LAMP2", "GRN", "CTSD",
  # 5. Anti-inflammatory Response (The "Post-clearance" resolution)
  "TGFB1", "IL10", "ARG1", "NR4A1", "PPARG"))
```
```{r}
TAMs <- AddModuleScore(
  object = TAMs,
  features = Efferocytosis_signature,
  name = "Efferocytosis_Master_Score",
  ctrl = 100
)
# Rename for clarity (removing the '1' Seurat adds)
TAMs$Efferocytosis_activity <- TAMs$Efferocytosis_Master_Score1
```
# 18.2 MHC signature
```{r}
MHCII_signature <- list(c(
  # Core MHC II
  "HLA-DRA", "HLA-DRB1", "HLA-DRB5", 
  "HLA-DQA1", "HLA-DQB1", "HLA-DQA2", 
  "HLA-DQB2", "HLA-DPA1", "HLA-DPB1", 
  # Antigen processing/loading
  "CD74", "HLA-DMA", "HLA-DMB",
  # Master regulator
  "CIITA"
))
```
```{r}
TAMs <- AddModuleScore(
  object = TAMs,
  features = MHCII_signature,
  name = "MHCII_Score",
  ctrl = 100
)

# Clean column name
TAMs$MHCII_activity <- TAMs$MHCII_Score1
```
# 18.3 OXPHOS score
```{r}
# Define the OXPHOS gene list
oxphos_core <- c(
  # Complex I
  "NDUFS1","NDUFS2","NDUFS3","NDUFS6","NDUFS7",
  "NDUFA9","NDUFA10",
  "NDUFB5","NDUFB8",
  "NDUFV1",
  # Complex II
  "SDHA","SDHB","SDHC","SDHD",
  # Complex III
  "UQCRC1","UQCRC2","UQCRFS1","UQCRB","UQCRQ",
  # Complex IV
  "COX4I1","COX5A","COX5B","COX6C","COX6B1",
  # Complex V
  "ATP5F1A","ATP5F1B","ATP5F1C","ATP5MC1","ATP5O",
  # Electron transport
  "CYCS","CYC1"
)                                            

# Calculate the module score
TAMs <- AddModuleScore(
  object = TAMs,
  features = oxphos_core,
  name = "OXPHOS_Score",
  ctrl = 100
)

# Clean column name
TAMs$OXPHOS_activity <- TAMs$OXPHOS_Score1
```
# 18.4 Ribosome
```{r}
ribosome_terms <- gsea_GO_BP_EP4EQ@result %>%
  filter(str_detect(tolower(Description), "ribosom"))
```
```{r}
ribosome_genes <- ribosome_terms$core_enrichment %>%
  strsplit("/") %>%
  unlist() %>%
  unique()
ribosome_genes <- intersect(ribosome_genes, rownames(TAMs))
ribosome_signature <- list(ribosome_genes)
```
```{r}
TAMs <- AddModuleScore(
  object = TAMs,
  features = ribosome_signature,
  name = "Ribosome_Score",
  ctrl = 100
)
# Clean column name
TAMs$Ribosome_activity <- TAMs$Ribosome_Score1
```
