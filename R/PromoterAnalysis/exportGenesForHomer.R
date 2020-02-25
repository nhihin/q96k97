# Export out lists of IRE genes for HOMER

zebrafishIreGenes <- readRDS(here::here("R/IREGenes/data/zebrafishIreGenes.rds"))
data.frame(gene=zebrafishIreGenes$ire3_all) %>% write_tsv(here::here("R/IREGenes/data/zebrafishIreGenes_ire3all.txt"), col_names = FALSE)
data.frame(gene=zebrafishIreGenes$ire5_all) %>% write_tsv(here::here("R/IREGenes/data/zebrafishIreGenes_ire5all.txt"), col_names = FALSE)
data.frame(gene=zebrafishIreGenes$ire3_hq) %>% write_tsv(here::here("R/IREGenes/data/zebrafishIreGenes_ire3hq.txt"), col_names = FALSE)
data.frame(gene=zebrafishIreGenes$ire5_hq) %>% write_tsv(here::here("R/IREGenes/data/zebrafishIreGenes_ire5hq.txt"), col_names = FALSE)

humanIreGenes<-readRDS(here::here("R/IREGenes/data/human_ireGenes.rds"))
data.frame(gene=humanIreGenes$ire3_all) %>% write_tsv(here::here("R/IREGenes/data/humanIreGenes_ire3all.txt"), col_names = FALSE)
data.frame(gene=humanIreGenes$ire5_all) %>% write_tsv(here::here("R/IREGenes/data/humanIreGenes_ire5all.txt"), col_names = FALSE)
data.frame(gene=humanIreGenes$ire3_hq) %>% write_tsv(here::here("R/IREGenes/data/humanIreGenes_ire3hq.txt"), col_names = FALSE)
data.frame(gene=humanIreGenes$ire5_hq) %>% write_tsv(here::here("R/IREGenes/data/humanIreGenes_ire5hq.txt"), col_names = FALSE)

mouseIreGenes <- readRDS(here::here("R/Mouse/data/mouse_ireSets.rds"))
data.frame(gene=mouseIreGenes$ire3_all) %>% write_tsv(here::here("R/IREGenes/data/mouseIreGenes_ire3all.txt"), col_names = FALSE)
data.frame(gene=mouseIreGenes$ire5_all) %>% write_tsv(here::here("R/IREGenes/data/mouseIreGenes_ire5all.txt"), col_names = FALSE)
data.frame(gene=mouseIreGenes$ire3_hq) %>% write_tsv(here::here("R/IREGenes/data/mouseIreGenes_ire3hq.txt"), col_names = FALSE)
data.frame(gene=mouseIreGenes$ire5_hq) %>% write_tsv(here::here("R/IREGenes/data/mouseIreGenes_ire5hq.txt"), col_names = FALSE)