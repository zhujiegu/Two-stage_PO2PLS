library(readxl)
library(magrittr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(meffil)

####################################
## Methylation QC
####################################

# Load data
# Methylation data contains 485577 probes and 29 *3  samples
methyl <- data.table::fread(file = '~/Data/DOWN/GSE52588_series_matrix.txt', header = T, verbose = F, skip = 86) %>% arrange(ID_REF)
meta <- data.table::fread(file = '~/Data/DOWN/GSE52588_series_matrix.txt', header = T, verbose = F, skip = 47, nrows = 35)
probe <- data.table::fread(file = '~/Data/DOWN/GSE52588_family.soft', header = T, skip=215, verbose = F, na.strings="") %>% arrange(ID)
all.equal(probe$ID, methyl$ID_REF)

# # Quality control
# Filtering:
#   1. Samples with more than 75% of the probes with a detection p-value greater than 1e-05 (0)
# 2. Probes had a detection p-value greater than 0.05 in more than 75% samples (0, says 425 in the paper, should be removed already)
# 3. Drop CHR X and Y (11648)
# 4. Drop probes with missing (22948)
# 450981 remains (same as in the paper)

# detection pvalue
pvalue <- data.table::fread(file = '~/Data/DOWN/GSE52588_raw_data.txt', header = T, verbose = F) %>% arrange(ID_REF) %>% select(contains("Pval"))

sum(colSums(pvalue > 1e-5) > 0.75 * dim(pvalue)[1])
sum(rowSums(pvalue > 0.05) > 0.75 * dim(pvalue)[2])
# 0 sample, 0 probe removed

# drop chr X Y
methyl %<>% mutate(CHR = probe$CHR)
methyl %<>% filter(!CHR %in% c('X', 'Y'))  # 11648 removed

# drop missing
rrm_miss <- which(select(methyl, contains('GSM'))=='NULL', arr.ind = T)[,1] %>% unique  # 22948 removed
methyl <- methyl[-rrm_miss, ]
probe <- filter(probe, !CHR %in% c('X', 'Y'))[-rrm_miss, ]
rm(pvalue)

# Estimate cell count
meffil.list.cell.type.references()
cg_names <- methyl$ID_REF
methyl_meffil <- as.matrix(methyl[,-c(1,89)]) %>% apply(2,as.numeric)
rownames(methyl_meffil) <- cg_names

cellcount=meffil.estimate.cell.counts.from.betas(methyl_meffil,cell.type.reference ="blood gse35069",verbose = T )
rm(methyl_meffil)

# transpose
sample_names <- colnames(methyl)[-c(1,89)]
methyl <- methyl %>% dplyr::select(starts_with("GS")) %>% apply(1, as.numeric) %>% scale(scale = F)
colnames(methyl) <- cg_names
rownames(methyl) <- sample_names

####################################
## Glycomics
####################################

# Load data
# BASH command: awk '{print $1, $NF}' methy-glycan_name_full.txt > methy-glycan_name.txt
ds <- read_excel("/home/z/Data/DOWN/Glycan/Down_glicazione_CHITTY.xls", sheet = "down")
ma <- read_excel("/home/z/Data/DOWN/Glycan/Down_glicazione_CHITTY.xls", sheet = "mothers")
sb <- read_excel("/home/z/Data/DOWN/Glycan/Down_glicazione_CHITTY.xls", sheet = "siblings")
name_match <- read.delim("/home/z/Data/DOWN/methy-glycan_name.txt", sep = " ", header = F)
gly <- rbind(ds, ma, sb)
rm(ds,ma,sb)
gly %<>% filter(PID %in% name_match$V2)

# remove missing
rrm_miss <- which(is.na(gly %>% dplyr::select(c(P1:P10, GlycoAgeTest))), arr.ind = T)[,1] %>% unique
gly <- gly[-rrm_miss,]

# add same name as meth data
gly %<>% mutate(ID = name_match$V1[match(gly$PID, name_match$V2)])
all.equal(match(gly$PID, name_match$V2), match(gly$ID, name_match$V1))

# take log on TA data
gly[,8:17] <- gly[,8:17] %>% as.matrix() %>% log

####################################
## Overlapping
####################################

rrm_meth <- which(!rownames(methyl)%in%gly$ID)
methyl <- methyl[-rrm_meth,]
sample_names <- rownames(methyl)

# row order
gly <- gly[match(rownames(methyl), gly$ID),]
all.equal(rownames(methyl), as.character(gly$ID))

# gly data matrix ready
covt <- gly %>% dplyr::select(ID, Age, Sex)
gly %<>% dplyr::select(P1:P10) %>% as.matrix()

rownames(gly) <- rownames(methyl)
colnames(gly) <- c(paste('P',1:10, sep = ''))

####################################
## Covariate
####################################
cellcount <- cellcount[rownames(cellcount)%in%sample_names,]
all.equal(rownames(cellcount), rownames(methyl))

covt <- cbind(covt, cellcount)
covt$Sex %<>% as.factor()

# Save data
save(methyl, gly, covt, file = "~/Data/DOWN/Pre-integrate_unadj.RData")
