# Running IRE enrichment for karissa

# import objects
contrasts <- readRDS("~/Downloads/contrasts.rds")
voom <- readRDS("~/Downloads/voom_RUVk1.rds")
design <- readRDS("~/Downloads/designRUV_k1.rds")

# function
source(here::here("R/GSEA/combinedGSEA.R"))

# import IRE gene sets
ireGenes <- readRDS(here::here("R/IREGenes/data/zebrafishIreGenes.rds"))

# create index
idx <- limma::ids2indices(ireGenes, rownames(voom))

# run enrichment
enrRes <- combinedGSEA(v = voom, idx = idx, design = design, contrasts = contrasts)

# export results
saveRDS(enrRes, "~/Downloads/enrRes_karissa.rds")

# results (combined p-values for each contrast)
enrRes$combTest