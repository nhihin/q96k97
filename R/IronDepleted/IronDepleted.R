# library(hugene10sttranscriptcluster.db)
library(hgu133plus2.db)
library(affy)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(export)
theme_set(theme_bw())

# PROCESS RAW MICROARRAY DATA AND RETRIEVE PROBE-GENE ANNOTATIONS --------------------

# Import raw data downloaded from GEO
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE11670
affybatch <- ReadAffy(celfile.path = here::here("R/IronDepleted/data/GSE11670_RAW"))

# Convert from raw probe intensities to expression 
normData <- expresso(affybatch, 
                     bgcorrect.method = "rma", 
                     normalize.method = "constant", 
                     pmcorrect.method = "pmonly", 
                     summary.method = "avgdiff")

# Expression matrix
exprs <- normData@assayData$exprs

# Samples table from GEO
geo <- getGEO("GSE11670", GSEMatrix = TRUE)
samples <- geo$GSE11670_series_matrix.txt.gz%>% 
  phenoData() %>% 
  pData %>% 
  as_tibble %>%
  dplyr::select(title, 
                geo_accession,
                source_name_ch1) %>% 
  dplyr::rename(
    sample = title, 
    id = geo_accession, 
    treatment = source_name_ch1
  ) %>%
  dplyr::mutate(
    treatment = c("Control", "Control", 
                  "ICL670_10uM", "ICL670_10uM",
                  "ICL670_50uM", "ICL670_50uM")
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
pca_plot_GSE11670 <- pca$x %>% as.data.frame %>% 
  rownames_to_column("filename") %>% 
  left_join(samples%>%rownames_to_column("filename")) %>%
  ggplot(aes(x = PC1, y = PC2, colour = treatment)) + 
  geom_point(size=3) + 
  geom_text_repel(aes(label = filename)) +
  theme(aspect.ratio = 1) +
  labs(x = "Principal Component 1 (40.3%)", 
       y = "Principal Component 2 (20.6%)",
       colour = "Treatment")
pca_plot_GSE11670
export::graph2ppt(pca_plot_GSE11670, here::here("R/IronDepleted/fig/pca_plot_GSE11670"))
saveRDS(pca_plot_GSE11670, here::here("R/IronDepleted/fig/pca_plot_GSE11670.rds"))
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
annot %<>% dplyr::filter(!ENSEMBL == "")

head(annot)
dim(annot) # 43,285 probes

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
  IronChelator_ICL670_10uM = ICL670_10uM - Control,
  IronChelator_ICL670_50uM = ICL670_50uM - Control,
  DosageEffect = ICL670_50uM - ICL670_10uM
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
# dim(colExprs)  # 23788 genes

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
saveRDS(ireEnr, here::here("R/IronDepleted/data/GSE11670_ireEnr.rds"))
