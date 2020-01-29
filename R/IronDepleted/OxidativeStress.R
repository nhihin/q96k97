# Cells treated to oxidative stress via paraquat
# Load packages
library(oligo)
library(huex10sttranscriptcluster.db)
library(pd.huex.1.0.st.v2)
library(magrittr)
library(tibble)
library(dplyr)

# Import raw data from GEO ----------------------------------------------------------------
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21305
celFiles <- list.files(here::here("R/IronDepleted/data/GSE21305_RAW"), full.names=TRUE)
rawData <- read.celfiles(celFiles)
normData <- rma(rawData)

# Expression matrix
exprs <- normData@assayData$exprs
plotDensity(exprs)

# Check whether probe IDs are in the annotation package
rownames(exprs) %in% keys(huex10sttranscriptcluster.db) %>% summary  # all TRUE

# Retrieve gene id annotations using the probe ids
annot <- AnnotationDbi::select(
  x = huex10sttranscriptcluster.db, 
  keys = rownames(exprs),
  columns = c("PROBEID", "ENSEMBL", "ENTREZID", "SYMBOL"),
  keytype = "PROBEID"
)
all(rownames(exprs) %in% annot$PROBEID)  # TRUE

# Import sample data --------------------------------------------------------------------
geo <- getGEO("GSE21305", GSEMatrix = TRUE)
samples <- geo$GSE21305_series_matrix.txt.gz%>% 
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
  dplyr::arrange(treatment)  %>%
  dplyr::mutate(
    treatment = c(rep("Control", 5), rep("Paraquat", 5))
  ) 
samples <- normData@phenoData@data %>%
  rownames_to_column("filename") %>%
  mutate(id = gsub(x = filename, pattern = ".CEL.gz", replacement = "")) %>%
  dplyr::select(-index) %>%
  left_join(samples, by="id") %>%
  dplyr::arrange(treatment) %>%
  column_to_rownames("filename") 

# Principal Component Analysis ----------------------------------------------------------
pca <- exprs %>% log2 %>% t %>% prcomp()
summary(pca)
pca_GSE21305 <- pca$x %>% as.data.frame %>% 
  rownames_to_column("filename") %>% 
  left_join(samples%>%rownames_to_column("filename")) %>%
  ggplot(aes(x = PC1, y = PC2, colour = treatment)) + 
  geom_point(size = 3) + 
  geom_text_repel(aes(label = filename)) +
  theme(aspect.ratio = 1) +
  labs(x = "Principal Component 1 (29.9%)", 
       y = "Principal Component 2 (16.1%)", 
       colour = "Treatment")
pca_GSE21305
saveRDS(pca_GSE21305, here::here("R/IronDepleted/fig/pca_GSE21305.rds"))
export::graph2ppt(pca_GSE21305, here::here("R/IronDepleted/fig/pca_GSE21305"))

# DE Analysis ------------------------------------------------------------------------------
design <- model.matrix(~0 + treatment, data = samples) %>%
  set_colnames(gsub(x = colnames(.), pattern = "treatment", replacement = ""))

contrasts <- makeContrasts(
  levels = colnames(design), 
  OxidativeStress = Paraquat - Control
)

fit <- lmFit(exprs, design) %>%
  contrasts.fit(contrasts)%>%
  eBayes()

# Uniform p-value distribution
hist(fit$p.value)

results <- decideTests(fit, 
                       p.value = 0.05, 
                       adjust.method = "fdr")

# No significantly DE genes
results %>% summary()

# Volcano plot checks out, but log fold changes in general are v small
volcanoplot(fit)

# Gene level analysis ------------------------------------------------------------------------
# Combine the expression data with annotations
datET <- exprs %>% 
  as.data.frame %>% 
  rownames_to_column("PROBEID") %>%
  left_join(annot, by = "PROBEID") %>% 
  dplyr::filter(!is.na(ENSEMBL) | ENSEMBL != "") %>%
  dplyr::distinct(PROBEID, .keep_all = TRUE)%>%
  set_rownames(.$PROBEID)
 
# Collapse rows
collapsed <- WGCNA::collapseRows(datET[,2:11], 
                                 rowGroup = datET$ENSEMBL, 
                                 rowID = datET$PROBEID, 
                                 method ="MaxMean")

colExprs <- collapsed$datETcollapsed
# dim(colExprs)  # 17139 genes

colAnnot <- data.frame(ENSEMBL = rownames(colExprs)) %>%
  left_join(annot, by = "ENSEMBL") %>%
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

newFit <- lmFit(colExprs, design) %>%
  contrasts.fit(contrasts) %>%
  eBayes()

results <- decideTests(newFit, 
                       p.value = 0.05, 
                       adjust.method = "fdr")

results %>% summary()  # Still nothing at gene level


# IRE Enrichment Analysis -----------------------------------------------------------------

source(here::here("R/GSEA/combinedGSEA_ma.R"))
humanIreGenes<-readRDS(here::here("R/IREGenes/data/human_ireGenes.rds"))
humanIreIdx <- limma::ids2indices(humanIreGenes, identifiers = rownames(colExprs))

ireEnr <- combinedGSEA_ma(log2(colExprs), 
                          fit = newFit,
                          design = design,
                          contrasts = contrasts, 
                          idx = humanIreIdx)

# export results
saveRDS(ireEnr, here::here("R/IronDepleted/data/GSE21305_ireEnr.rds"))
