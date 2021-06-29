# Covariate adjust

library(magrittr)
library(dplyr)
library(ggplot2)
library(tidyr)
load('~/Data/DOWN/Pre-integrate_unadj.RData')

# Methylation
if(!all.equal(covt$ID %>% as.character(), rownames(methyl))) stop("rownames do not match")

library(parallel)

correct <- mclapply(1:ncol(methyl), mc.cores = 4, function(i){
  y <- methyl[,i]
  fit <- lm(data = covt, y~.-ID, na.action = na.exclude)
  return(residuals(fit))
})

meth <- do.call(cbind, correct)
colnames(meth) <- colnames(methyl)
all.equal(rownames(meth), rownames(methyl))

# Glycomics
if(!all.equal(covt$ID %>% as.character(), rownames(gly))) stop("rownames do not match")

correct <- lapply(1:ncol(gly), function(i){
  y <- gly[,i]
  fit <- lm(data = covt, y~Age+Sex, na.action = na.exclude)
  return(residuals(fit))
})

gly <- do.call(cbind, correct)
colnames(gly) <- c(paste('P',1:10, sep = ''))
all.equal(rownames(gly), rownames(meth))

probe <- data.table::fread(file = '~/Data/DOWN/GSE52588_family.soft', header = T, skip=215, verbose = F, na.strings="") %>% arrange(ID)
chr <- probe %>% filter(ID %in% colnames(meth)) %>% select(ID, CHR)
all.equal(chr$ID, colnames(meth))
save(meth, gly, chr, file = "~/Data/DOWN/Pre-integrate_adj.RData")

##############################################################################
# leave out DSM, don't correct for age
library(magrittr)
library(dplyr)
library(ggplot2)
library(tidyr)
load('~/Data/DOWN/Pre-integrate_unadj.RData')
methyl <- methyl[1:56,]
gly <- gly[1:56,]
covt <- covt[1:56,]

# Methylation
if(!all.equal(covt$ID %>% as.character(), rownames(methyl))) stop("rownames do not match")

library(parallel)

correct <- mclapply(1:ncol(methyl), mc.cores = 4, function(i){
  y <- methyl[,i]
  fit <- lm(data = covt, y~.-ID-Age, na.action = na.exclude)
  return(residuals(fit))
})

meth <- do.call(cbind, correct)
colnames(meth) <- colnames(methyl)
all.equal(rownames(meth), rownames(methyl))

# Glycomics
if(!all.equal(covt$ID %>% as.character(), rownames(gly))) stop("rownames do not match")

correct <- lapply(1:ncol(gly), function(i){
  y <- gly[,i]
  fit <- lm(data = covt, y~Sex, na.action = na.exclude)
  return(residuals(fit))
})

gly <- do.call(cbind, correct)
colnames(gly) <- c(paste('P',1:10, sep = ''))
all.equal(rownames(gly), rownames(meth))

probe <- data.table::fread(file = '~/Data/DOWN/GSE52588_family.soft', header = T, skip=215, verbose = F, na.strings="") %>% arrange(ID)
chr <- probe %>% filter(ID %in% colnames(meth)) %>% select(ID, CHR)
all.equal(chr$ID, colnames(meth))
save(meth, gly, chr, covt, file = "~/Data/DOWN/Pre-integrate_DSP_DSS_adj.RData")

##############################################
# Methylation only cell count
load('~/Data/DOWN/Pre-integrate_unadj.RData')
library(parallel)

# Adjust for cell count only
correct <- mclapply(1:ncol(methyl), mc.cores = 4, function(i){
  y <- methyl[,i]
  fit <- lm(data = covt, y~.-ID - Age - Sex, na.action = na.exclude)
  return(residuals(fit))
})

meth <- do.call(cbind, correct)
colnames(meth) <- colnames(methyl)
all.equal(rownames(meth), rownames(methyl))

probe <- data.table::fread(file = '~/Data/DOWN/GSE52588_family.soft', header = T, skip=215, verbose = F, na.strings="") %>% arrange(ID)
chr <- probe %>% filter(ID %in% colnames(meth)) %>% select(ID, CHR)
all.equal(chr$ID, colnames(meth))
save(meth, gly, chr, covt, file = "~/Data/DOWN/Pre-integrate_count_adj.RData")
