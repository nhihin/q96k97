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


