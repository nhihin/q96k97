library(Biobase)
library(GEOquery)
library(hgu133a.db)
library(hgu133plus2.db)
library(limma)


# ------------------------------- Attempt 1 --------------------------------


# Download normalised expression matrix from GEO
# (Authors didn't make raw data available)
gset <- getGEO("GSE3045", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# Expression matrix
exprs <- gset@assayData$exprs
exprs %>% log2 %>% boxplot()

# Sample table
samples <- gset@phenoData@data %>%
  as.data.frame %>%
  dplyr::select(title, geo_accession) %>%
  dplyr::mutate(treatment = c(rep("hypoxia", 3), rep("normoxia",3))) %>%
  set_rownames(.$geo_accession)

# Initial PCA to check if things look legit
pca <- exprs %>% log2 %>% t %>% prcomp()
pca$x %>% 
  as.data.frame %>% 
  rownames_to_column("geo_accession") %>%
  left_join(samples, by = "geo_accession") %>%
  ggplot(aes(x = PC1, y = PC2, colour = treatment)) +
  geom_point(size=3) +
  geom_text_repel(aes(label = title))

# No, things did not look legit. 


# ------------------------------- Attempt 2 --------------------------------

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE12792
## Import raw data -----------------------------------------------------------
affybatch <- ReadAffy(celfile.path = 
                        here::here("R/IronDepleted/data/GSE12792_RAW"))

# Convert from raw probe intensities to expression 
normData <- expresso(affybatch, 
                     bgcorrect.method = "rma", 
                     normalize.method = "constant", 
                     pmcorrect.method = "pmonly", 
                     summary.method = "avgdiff")

# Expression matrix
exprs <- normData@assayData$exprs

## Samples table from GEO ------------------------------------------------
geo <- getGEO("GSE12792", GSEMatrix = TRUE)
samples <- geo$GSE12792_series_matrix.txt.gz%>% 
  phenoData() %>% 
  pData %>% 
  as_tibble %>%
  dplyr::select(title, 
                geo_accession) %>% 
  dplyr::rename(
    sample = title, 
    id = geo_accession
  ) %>%
  dplyr::mutate(
    treatment = c(rep("normoxia",3), rep("hypoxia",3))
  ) 
samples <- normData@phenoData@data %>%
  rownames_to_column("filename") %>%
  mutate(id = gsub(x = filename, 
                   pattern = ".CEL.gz", 
                   replacement = "")) %>%
  dplyr::select(-sample) %>%
  left_join(samples, by="id") %>%
  dplyr::arrange(treatment) %>%
  column_to_rownames("filename") 


## Principal component analysis --------------------------------------------
# This one kind of checks out
pca <- exprs %>% log2 %>% t %>% prcomp()
summary(pca)
pca_GSE12792 <- pca$x %>% 
  as.data.frame %>% 
  rownames_to_column("filename") %>% 
  left_join(samples %>%
              rownames_to_column("filename")) %>%
  ggplot(aes(x = PC1, y = PC2, colour = treatment)) + 
  geom_point(size = 3) + 
  geom_text_repel(aes(label = filename)) +
  theme(aspect.ratio = 1) +
  labs(x = "Principal Component 1 (30.6%)", 
       y = "Principal Component 2 (23.6%)",
       colour = "Treatment")
pca_GSE12792
export::graph2ppt(pca_GSE12792, here::here("R/IronDepleted/fig/pca_GSE12792"))
saveRDS(pca_GSE12792, here::here("R/IronDepleted/fig/pca_GSE12792.rds"))

## Annotation ---------------------------------------------------------------
# Check that the probe IDs are in the annotation package 
rownames(exprs) %in% keys(hgu133a.db) %>% summary  # all TRUE

# Retrieve gene id annotations
annot <- AnnotationDbi::select(
  x = hgu133a.db, 
  keys = rownames(exprs),
  columns = c("PROBEID", "ENSEMBL", "ENTREZID", "SYMBOL"),
  keytype = "PROBEID"
)
all(rownames(exprs) %in% annot$PROBEID)  # TRUE

# Check for duplicate probe ids (assigned to different genes)
dup.ids <- annot$PROBEID[duplicated(annot$PROBEID)] %>% 
  unique %>%
  sort
length(dup.ids)  # 4430

# Example probe
annot[annot$PROBEID == dup.ids[2], ]

# Collapse down the duplicated genes so each row is 
# a unique probe ID
collapser <- function(x){
  x %>% unique %>% sort %>% paste(collapse = "|")
}
annot <- annot %>% 
  group_by(PROBEID) %>%
  summarise_each(funs(collapser)) %>% 
  ungroup

# Filter out probes which don't have ensembl IDs
annot %<>% dplyr::filter(!ENSEMBL == "")

head(annot)
dim(annot) # 43,285 probes

# We also would want to restrict the expression data to these probes
# and get the ordering to match
exprs <- normData@assayData$exprs %>% 
  as.data.frame %>%
  rownames_to_column("PROBEID") %>% 
  dplyr::filter(PROBEID %in% annot$PROBEID) %>%
  column_to_rownames("PROBEID") %>%
  magrittr::extract(annot$PROBEID, )

# Check the ordering of rownames is exactly the same
(rownames(exprs)==annot$PROBEID) %>% summary  # all TRUE


## Probe level DE analysis --------------------------------------------------
design <- model.matrix(~0 + treatment, data = samples) %>%
  set_colnames(gsub(x = colnames(.), 
                    pattern = "treatment", 
                    replacement = ""))

contrasts <- makeContrasts(
  levels = colnames(design), 
  hypoxiaEffect = hypoxia - normoxia
)

fit <- lmFit(log2(exprs), design) %>%
  contrasts.fit(contrasts)%>%
  eBayes()

results <- decideTests(
  fit, 
  p.value = 0.05, 
  adjust.method = "fdr"
)

results %>% summary()

volcanoplot(fit)

## Gene level analysis -----------------------------------------------------
# Expression matrix back to uncollapsed
exprs <- normData@assayData$exprs
annot <- AnnotationDbi::select(
  x = hgu133a.db, 
  keys = rownames(exprs),
  columns = c("PROBEID", "ENSEMBL", "ENTREZID", "SYMBOL"),
  keytype = "PROBEID"
)
all(rownames(exprs) %in% annot$PROBEID)  # TRUE

# Combine the expression data with annotations
datET <- exprs %>% 
  as.data.frame %>% 
  rownames_to_column("PROBEID") %>%
  left_join(annot, by = "PROBEID") %>% 
  dplyr::distinct(PROBEID, .keep_all = TRUE) %>%
  set_rownames(.$PROBEID)

# Collapse rows
collapsed <- WGCNA::collapseRows(datET[,2:7], 
                                 rowGroup = datET$ENSEMBL, 
                                 rowID = rownames(datET), 
                                 method ="MaxMean")

colExprs <- collapsed$datETcollapsed
dim(colExprs)  # 12830 genes

colAnnot <- data.frame(ENSEMBL = rownames(colExprs)) %>%
  left_join(annot, by = "ENSEMBL") %>%
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

newFit <- lmFit(log2(colExprs), design) %>%
  contrasts.fit(contrasts) %>%
  eBayes()

results <- decideTests(newFit, 
                       p.value = 0.05, 
                       adjust.method = "fdr")

results %>% summary()

## IRE Enrichment Analysis -----------------------------------------------
source(here::here("R/GSEA/combinedGSEA_ma.R"))
humanIreGenes<-readRDS(here::here("R/IREGenes/data/human_ireGenes.rds"))

humanIreIdx <- limma::ids2indices(humanIreGenes, 
                                  identifiers = rownames(colExprs))

ireEnr <- combinedGSEA_ma(log2(colExprs), 
                          fit = newFit, 
                          design = design, 
                          contrasts = contrasts, 
                          idx = humanIreIdx)


ireEnr %>%saveRDS(here::here("R/IronDepleted/data/GSE12792_ireEnr.rds"))




# ------------------------------- Attempt 3 --------------------------------

# Unfortunately, the following dataset only has n=2 so it's probably 
# going to be dodgy
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE4086

# Import raw data
affybatch <- ReadAffy(celfile.path = here::here("R/IronDepleted/data/GSE4086_RAW"))

# Convert from raw probe intensities to expression 
normData <- expresso(affybatch, 
                     bgcorrect.method = "rma", 
                     normalize.method = "constant", 
                     pmcorrect.method = "pmonly", 
                     summary.method = "avgdiff")

# Expression matrix
exprs <- normData@assayData$exprs

# Samples table from GEO
geo <- getGEO("GSE4086", GSEMatrix = TRUE)
samples <- geo$GSE4086_series_matrix.txt.gz%>% 
  phenoData() %>% 
  pData %>% 
  as_tibble %>%
  dplyr::select(title, 
                geo_accession) %>% 
  dplyr::rename(
    sample = title, 
    id = geo_accession
  ) %>%
  dplyr::mutate(
    treatment = c(rep("hypoxia",2), rep("normoxia",2))
  ) 
samples <- normData@phenoData@data %>%
  rownames_to_column("filename") %>%
  mutate(id = gsub(x = filename, 
                   pattern = ".CEL.gz", 
                   replacement = "")) %>%
  dplyr::select(-sample) %>%
  left_join(samples, by="id") %>%
  dplyr::arrange(treatment) %>%
  column_to_rownames("filename") 


# Principal component analysis
# This one kind of checks out
pca <- exprs %>% log2 %>% t %>% prcomp()
summary(pca)
pca_GSE4086 <-pca$x %>% as.data.frame %>% 
  rownames_to_column("filename") %>% 
  left_join(samples %>%
              rownames_to_column("filename")) %>%
  ggplot(aes(x = PC1, y = PC2, colour = treatment)) + 
  geom_point(size=3) + 
  geom_text_repel(aes(label = filename))+
  theme(aspect.ratio = 1)+
  labs(x="Principal Component 1 (59.8%)",
       y="Principal Component 2 (27.4%)",
       colour="Treatment")
export::graph2ppt(pca_GSE4086, here::here("R/IronDepleted/fig/pca_GSE4086"))

# Check that the probe IDs are in the annotation package 
rownames(exprs) %in% keys(hgu133plus2.db) %>% summary  # all TRUE

# Retrieve gene id annotations
annot <- AnnotationDbi::select(
  x = hgu133plus2.db, 
  keys = rownames(exprs),
  columns = c("PROBEID", "ENSEMBL", "ENTREZID", "SYMBOL"),
  keytype = "PROBEID"
)
all(rownames(exprs) %in% annot$PROBEID)  # TRUE

# Check for duplicate probe ids (assigned to different genes)
dup.ids <- annot$PROBEID[duplicated(annot$PROBEID)] %>% 
  unique %>%
  sort
length(dup.ids)  # 4430

# Example probe
annot[annot$PROBEID == dup.ids[2], ]



# Limma analysis
design <- model.matrix(~0 + treatment, data = samples) %>%
  set_colnames(gsub(x = colnames(.), 
                    pattern = "treatment", 
                    replacement = ""))

contrasts <- makeContrasts(
  levels = colnames(design), 
  hypoxiaEffect = hypoxia - normoxia
)

fit <- lmFit(log2(exprs), design) %>%
  contrasts.fit(contrasts)%>%
  eBayes()

hist(fit$Amean)

results <- decideTests(
  fit, 
  p.value = 0.05, 
  adjust.method = "fdr"
)

results %>% summary()

volcanoplot(fit)

## Gene level analysis -----------------------------------------------------
# Expression matrix back to uncollapsed
exprs <- normData@assayData$exprs
annot <- AnnotationDbi::select(
  x = hgu133plus2.db, 
  keys = rownames(exprs),
  columns = c("PROBEID", "ENSEMBL", "ENTREZID", "SYMBOL"),
  keytype = "PROBEID"
)
all(rownames(exprs) %in% annot$PROBEID)  # TRUE

# Combine the expression data with annotations
datET <- exprs %>% 
  as.data.frame %>% 
  rownames_to_column("PROBEID") %>%
  left_join(annot, by = "PROBEID") %>% 
  dplyr::distinct(PROBEID, .keep_all = TRUE) %>%
  set_rownames(.$PROBEID)

# Collapse rows
collapsed <- WGCNA::collapseRows(datET[,2:5], 
                                 rowGroup = datET$ENSEMBL, 
                                 rowID = rownames(datET), 
                                 method ="MaxMean")

colExprs <- collapsed$datETcollapsed
dim(colExprs)  # 20143

colAnnot <- data.frame(ENSEMBL = rownames(colExprs)) %>%
  left_join(annot, by = "ENSEMBL") %>%
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

newFit <- lmFit(log2(colExprs), design) %>%
  contrasts.fit(contrasts) %>%
  eBayes()

results <- decideTests(newFit, 
                       p.value = 0.05, 
                       adjust.method = "fdr")

results %>% summary()



# IRE enrichment analysis
source(here::here("R/GSEA/combinedGSEA_ma.R"))
humanIreGenes<-readRDS(here::here("R/IREGenes/data/human_ireGenes.rds"))
humanIreIdx <- limma::ids2indices(humanIreGenes, 
                                  identifiers = rownames(colExprs))

ireEnr <- combinedGSEA_ma(log2(colExprs), 
                          fit = newFit, 
                          design = design, 
                          contrasts = contrasts, 
                          idx = humanIreIdx)

# Attempt 4 ----------------------------------------------------------------------

# Import raw data downloaded from GEO
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22282
affybatch <- ReadAffy(celfile.path = here::here("R/IronDepleted/data/GSE22282_RAW"))

# Convert from raw probe intensities to expression 
normData <- expresso(affybatch, 
                     bgcorrect.method = "rma", 
                     normalize.method = "constant", 
                     pmcorrect.method = "pmonly", 
                     summary.method = "avgdiff")

# Expression matrix
exprs <- normData@assayData$exprs

# Samples table from GEO
geo <- getGEO("GSE22282", GSEMatrix = TRUE)
samples <- geo$GSE22282_series_matrix.txt.gz%>% 
  phenoData() %>% 
  pData %>% 
  as_tibble %>%
  dplyr::select(title, 
                geo_accession,
                `treatment group:ch1`) %>% 
  dplyr::rename(
    sample = title, 
    id = geo_accession, 
    treatment = `treatment group:ch1`
  ) %>%
  dplyr::mutate(
    treatment = c(rep("hypoxia", 3), rep("normoxia", 3))
  ) 
# Join this onto the file names (sample names in normData)
samples <- normData@phenoData@data %>%
  rownames_to_column("filename") %>%
  mutate(id = gsub(x = filename, pattern = ".CEL.gz", replacement = "")) %>%
  dplyr::select(-sample) %>%
  left_join(samples, by="id") %>%
  dplyr::arrange(treatment) %>%
  column_to_rownames("filename") 


# Principal component analysis
# This one checks out
pca <- exprs %>% log2 %>% t %>% prcomp()
summary(pca)
pca_plot_GSE22282<- pca$x %>% as.data.frame %>% 
  rownames_to_column("filename") %>% 
  left_join(samples%>%rownames_to_column("filename")) %>%
  ggplot(aes(x = PC1, y = PC2, colour = treatment)) + 
  geom_point(size=3) + 
  geom_text_repel(aes(label = filename)) +
  theme(aspect.ratio = 1) +
  labs(x = "Principal Component 1 (63.8%)", 
       y = "Principal Component 2 (19.5%)",
       colour = "Treatment")
pca_plot_GSE22282
export::graph2ppt(pca_plot_GSE22282, here::here("R/IronDepleted/fig/pca_plot_GSE22282"))
saveRDS(pca_plot_GSE22282, here::here("R/IronDepleted/fig/pca_plot_GSE22282"))
pca_plot_GSE11670 <- readRDS(here::here("R/IronDepleted/fig/pca_plot_GSE11670.rds"))



# Annotation ----------------------------------------------------------------------
# Check that the probe IDs are in the annotation package 
rownames(exprs) %in% keys(hgu133plus2.db) %>% summary  # all TRUE

# Retrieve gene id annotations
annot <- AnnotationDbi::select(
  x = hgu133plus2.db, 
  keys = rownames(exprs),
  columns = c("PROBEID", "ENSEMBL", "ENTREZID", "SYMBOL"),
  keytype = "PROBEID"
)
all(rownames(exprs) %in% annot$PROBEID)  # TRUE

# Check for duplicate probe ids (assigned to different genes)
dup.ids <- annot$PROBEID[duplicated(annot$PROBEID)] %>% 
  unique %>%
  sort
length(dup.ids)  # 4430

# Example probe
annot[annot$PROBEID == dup.ids[2], ]

# Collapse down the duplicated genes so each row is 
# a unique probe ID
collapser <- function(x){
  x %>% unique %>% sort %>% paste(collapse = "|")
}
annot <- annot %>% 
  group_by(PROBEID) %>%
  summarise_each(funs(collapser)) %>% 
  ungroup

# Filter out probes which don't have ensembl IDs
# annot %<>% dplyr::filter(!ENSEMBL == "")

head(annot)
dim(annot) # 54675 probes

# We also would want to restrict the expression data to these probes
# and get the ordering to match
exprs <- normData@assayData$exprs %>% as.data.frame %>%
  rownames_to_column("PROBEID") %>% 
  dplyr::filter(PROBEID %in% annot$PROBEID) %>%
  column_to_rownames("PROBEID") %>%
  magrittr::extract(annot$PROBEID, )

# Check the ordering of rownames is exactly the same
(rownames(exprs)==annot$PROBEID) %>% summary  # all TRUE


# LIMMA ANALYSIS -------------------------------------------------------------------------------------------------------
# Bg: Need the following objects from limma analysis to perform IRE enrichment analysis:
#     voomData, design, contrasts

design <- model.matrix(~0 + treatment, data = samples) %>%
  set_colnames(gsub(x = colnames(.), pattern = "treatment", replacement = ""))

contrasts <- makeContrasts(
  levels = colnames(design), 
  Hypoxia = hypoxia - normoxia
)

fit <- lmFit(log2(exprs), design) %>%
  contrasts.fit(contrasts)%>%
  eBayes()

hist(fit$Amean) # Many genes have low expression

results <- decideTests(fit, 
                       p.value = 0.05, 
                       adjust.method = "fdr")

results %>% summary()

volcanoplot(fit)

topgenes <- topTable(fit, coef = 1, n=1000) %>%
  as.data.frame %>%
  rownames_to_column("PROBEID") %>%
  left_join(annot)

# Gene level analysis -------------------------------------------------

# Expression matrix back to uncollapsed
exprs <- normData@assayData$exprs
annot <- AnnotationDbi::select(
  x = hgu133plus2.db, 
  keys = rownames(exprs),
  columns = c("PROBEID", "ENSEMBL", "ENTREZID", "SYMBOL"),
  keytype = "PROBEID"
)
all(rownames(exprs) %in% annot$PROBEID)  # TRUE

# Combine the expression data with annotations
datET <- exprs %>% 
  as.data.frame %>% 
  rownames_to_column("PROBEID") %>%
  left_join(annot, by = "PROBEID") %>% 
  dplyr::distinct(PROBEID, ENSEMBL, .keep_all = TRUE) 

# Collapse rows
collapsed <- WGCNA::collapseRows(datET[,2:7], 
                                 rowGroup = datET$ENSEMBL, 
                                 rowID = rownames(datET), 
                                 method ="MaxMean")

colExprs <- collapsed$datETcollapsed
 dim(colExprs)  # 23788  genes

colAnnot <- data.frame(ENSEMBL = rownames(colExprs)) %>%
  left_join(annot, by = "ENSEMBL") %>%
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

newFit <- lmFit(log2(colExprs), design) %>%
  contrasts.fit(contrasts) %>%
  eBayes()

results <- decideTests(newFit, 
                       p.value = 0.05, 
                       adjust.method = "fdr")

results %>% summary()  # still nothing


# IRE ENRICHMENT TESTING --------------------------------------------------------------

# Import in the enrichment function (fgsea; camera; mroast/fry) and the human gene sets
# with ENSEMBL ids. 
source(here::here("R/GSEA/combinedGSEA_ma.R"))
humanIreGenes <- readRDS(here::here("R/IREGenes/data/human_ireGenes.rds"))

# Build index specific to the ordering of the rownames in the expression matrix
humanIreIdx <- limma::ids2indices(humanIreGenes, 
                                  identifiers = rownames(colExprs))

ireEnr <- combinedGSEA_ma(log2(colExprs), 
                          fit = newFit, 
                          design = design, 
                          contrasts = contrasts, 
                          idx = humanIreIdx)
# Results
ireEnr

# Export
saveRDS(ireEnr, here::here("R/IronDepleted/data/GSE22282_ireEnr.rds"))


