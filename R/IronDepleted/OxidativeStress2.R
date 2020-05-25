library(readr)

# This is for the human autism RNA-seq dataset

# Raw counts downloaded from GEO
# Import the counts 
counts <- read_tsv(here::here("R/IronDepleted/data/GSE143155_Rawcounts_ASD_study.txt"))
str(counts)  # 46478 genes
counts$Geneid  # Gene symbol

# Sample names
sampleNames <- counts[,-c(1)] %>% colnames
sampleNames

# Retrieve sample metadata from GEO
geo <- getGEO("GSE143155", GSEMatrix = TRUE)
samples <- geo$GSE143155_series_matrix.txt.gz@phenoData@data %>%
  dplyr::select(title, 
                geo_accession,
                source_name_ch1, 
                starts_with("characteristics"),
                -characteristics_ch1.3) %>%
  dplyr::rename(
    cell_type = source_name_ch1, 
    patient_id = characteristics_ch1,
    autism = characteristics_ch1.1,
    sex = characteristics_ch1.2, 
    treatment = characteristics_ch1.4
  ) %>%
  dplyr::mutate(
    patient_id = gsub(x=patient_id, pattern="^.* ", replacement = ""),
    autism = gsub(x=autism, pattern = "^.* ", replacement =""),
    sex = gsub(x = sex, pattern = ".* ", replacement = ""),
    treatment = gsub(x = treatment, pattern = ".*: ", replacement = "")
  ) %>%
  dplyr::filter(cell_type == "neuron") %>%
  dplyr::arrange(treatment) %>%
  dplyr::mutate(treatment = c(rep("hydrogen_peroxide", 12), rep("control", nrow(.)-12))) 
samples

# Reorder counts so the samples are in the same order. 
counts <- counts[, c("Geneid", as.character(samples$title))] %>% 
  as.data.frame %>%
  dplyr::distinct(Geneid, .keep_all = TRUE) %>% 
  column_to_rownames("Geneid") 

# Gene annotation using AnnotationHub

ah <- AnnotationHub()
ah %>%
  subset(grepl("sapiens", species)) %>%
  subset(rdataclass == "EnsDb")
ensDb <- ah[["AH64923"]]
genes <- genes(ensDb)
genes <- as.data.frame(genes)
head(genes)

annot <- data.frame(gene_name = rownames(counts)) %>%
  left_join(genes, by = "gene_name") %>%
  dplyr::distinct(gene_name, .keep_all = TRUE) %>%
  dplyr::filter(gene_id != "" | !is.na(gene_id)) %>%
  set_rownames(.$gene_id)

counts2 <- counts %>% as.data.frame %>% 
  rownames_to_column("gene_name") %>%
  left_join(annot[, c("gene_name", "gene_id")]) %>%
  dplyr::filter(gene_id != "" | !is.na(gene_id)) %>%
  dplyr::distinct(gene_id, .keep_all=TRUE) %>% 
  dplyr::select(-gene_name) %>%
  column_to_rownames("gene_id")


# DGEList object
dge <- DGEList(
  counts = counts2,
  samples = samples,
  genes = annot,
  remove.zeros = TRUE
) %>%
  calcNormFactors("TMM")


# Principal Component Analysis

pca <- dge %>% cpm(log=TRUE) %>% t %>% prcomp()

shape_factor <- paste0(samples$autism, "_", samples$treatment) %>% 
  factor(levels = unique(.))

pca_plot <- pca$x %>% magrittr::extract(, c("PC1", "PC2")) %>%
  set_colnames(c("PCa", "PCb")) %>%
  as.data.frame %>%
  rownames_to_column("samples") %>%
  left_join((dge$samples %>% rownames_to_column("samples")), by="samples") %>%
  ggplot(aes(x=PCa, y = PCb, shape = shape_factor, colour = as.factor(patient_id))) +
  geom_point(alpha = 0.7,size=4) + 
  scale_shape_manual(values = c(0,15,1,16),
                     labels = levels(shape_factor)) +
  labs(x = "Principal Component 1 (38.5%)", 
       y = "Principal Component 2 (18.9%)",
       colour = "Patient ID", 
       shape = "Autism") +
  theme(aspect.ratio = 1)
pca_plot

saveRDS(pca_plot, here::here("R/IronDepleted/fig/autism_pca_plot.rds"))
export::graph2ppt(pca_plot, here::here("R/IronDepleted/fig/autism_pca_plot"))

pca_plot2 <- pca$x %>% magrittr::extract(, c("PC2", "PC3")) %>%
  set_colnames(c("PCa", "PCb")) %>%
  as.data.frame %>%
  rownames_to_column("samples") %>%
  left_join((dge$samples %>% rownames_to_column("samples")), by="samples") %>%
  ggplot(aes(x=PCa, y = PCb, shape = shape_factor, colour = as.factor(patient_id))) +
  geom_point(alpha = 0.7,size=4) + 
  scale_shape_manual(values = c(0,15,1,16),
                     labels = levels(shape_factor)) +
  labs(x = "Principal Component 2 (18.9%)", 
       y = "Principal Component 3 (11.3%)",
       colour = "Patient ID", 
       shape = "Autism") +
  theme(aspect.ratio = 1)
pca_plot2
saveRDS(pca_plot2, here::here("R/IronDepleted/fig/autism_pca_plot2.rds"))
export::graph2ppt(pca_plot2, here::here("R/IronDepleted/fig/autism_pca_plot2"))

# PCA Notes
# Samples group by the patient they were taken from. 
# Some separation between autism and non autism across PC2. 
# Some separation between hydrogen peroxide and control possibly. 
# This might be clearer across PC 3 although it's not perfect. 

dge$samples %<>%
  dplyr::mutate(autism = gsub(x = autism, pattern = "\\+", replacement ="Autism"),
                autism=gsub(x=autism,pattern="\\-",replacement="NoAutism"))

dge$samples %<>% mutate(group = as.factor(paste0(treatment, "_", autism)))

design <- model.matrix(~0 + group, data = dge$samples) %>%
  set_colnames(gsub(x = colnames(.), pattern = "group", replacement = ""))

# Apply the voom method, which transforms discrete count data into continuous log-normal distribution. 
voomData <- voom(dge, design = design, plot = TRUE)

# dupcors <- duplicateCorrelation(voomData, design, block=dge$samples$patient_id)
dupcors <- 0.5833879

# Define the contrasts (comparisons) which we will test for differential expression.
contrasts <- makeContrasts(
  levels = colnames(design), 
  
  effectOfAutism = control_Autism-control_NoAutism,
  effectOfAutism_withH2O2 = hydrogen_peroxide_Autism-hydrogen_peroxide_NoAutism,
  
  effectOfH2O2 = hydrogen_peroxide_NoAutism-control_NoAutism,
  effectOfH2O2_inAutism = hydrogen_peroxide_Autism-control_Autism
)

# Do moderated t-test for each gene to test for differential expression. 
# The empirical Bayes step "borrows" information across genes to increase accuracy of variance estimation.
fit <- lmFit(voomData, design, block = dge$samples$patient_id, correlation = dupcors) %>%
  contrasts.fit(contrasts) %>%
  eBayes(robust = TRUE)

fit %>% saveRDS("R/IronDepleted/data/autism_fit.rds")
plotSA(fit, main = "Final model: Mean-variance trend")

# Apply fdr adjustment and adjusted p-value cutoff of 0.05 to define significance.
results <- decideTests(fit, p.value = 0.05, adjust.method = "fdr", method="global")
summary(results)

topTable(fit, coef = "effectOfH2O2", n=1000) %>% View


# Enrichment
humanIreGenes<-readRDS(here::here("R/IREGenes/data/human_ireGenes.rds"))
human_ire_idx <- ids2indices(humanIreGenes, rownames(voomData))

gseaResults_humanIREs <- combinedGSEA3(voomData, human_ire_idx, design, contrasts, fit)

gseaResults_humanIREs %>% saveRDS(here::here("R/IronDepleted/data/ireEnr_autism.rds"))
