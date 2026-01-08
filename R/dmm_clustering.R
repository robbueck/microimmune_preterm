library(tidyverse)
library(ggplot2)
library(readxl)
library(phyloseq)
library(microbiome)
library(vegan)
library(metadeconfoundR)
library(DirichletMultinomial)
library(parallel)
library(gridExtra)
library(ggalluvial)
library(caret)
library(ggrepel)
library(patchwork)
library(mclust)
library(kml)
library(furrr)
library(mice)
library(igraph)
source("/fast/AG_Forslund/rob/PROSPER/immune_paper_data/R/functions.R")
setwd("/fast/AG_Forslund/rob/PROSPER/immune_paper_data/R/")

# switches ##################
dmm_step <- F
set.seed(123)


# load data:
ps_prosper_rf_filter <- readRDS(file = "/fast/AG_Forslund/rob/PROSPER/dada2/phyloseq_combined_rf_filter_genus.rds") 


immune_samples <- read_xlsx("/fast/AG_Forslund/rob/PROSPER/metadata/Hackathon_master_table_metadata_reduced_20250725_for RB.xlsx",
                            skip = 1,
                            sheet = "PBMC_chipcytometry_PRIMAL_IPOP")
immune_individuals <- immune_samples$MHH_PatID %>% unique
immune_data_abs <- read_xlsx("/fast/AG_Forslund/rob/PROSPER/immune_data/frequencies_of_mct.xlsx") %>%
  dplyr::select(-id, -timepoint) %>%
  column_to_rownames("chip") %>%
  as.matrix()

# get intersection of samples
keep_samples <- intersect(rownames(immune_data_abs), ps_prosper_rf_filter@sam_data$matching_chip_ID..chip_ID.)
kepp_individuals <- ps_prosper_rf_filter@sam_data %>% data.frame() %>% filter(matching_chip_ID..chip_ID. %in% keep_samples) %>%
  pull(MHH_PatID) %>%
  unique()
ps_prosper_rf_filter <- ps_prosper_rf_filter %>%
  subset_samples(MHH_PatID %in% kepp_individuals) %>%
  subset_samples(!Description %in% c("K301a", "K258a_2", "K295b_2")) # remove duplicated samples, keep the one closest to chip sample

# create time points
sample_data(ps_prosper_rf_filter) <- sample_data(ps_prosper_rf_filter) %>%
  data.frame() %>%
  mutate(time_point_immune = cut(day_of_life, breaks=c(-1, 16, 200, Inf),
                                 labels=c("1", "2", "3"))) %>%
  sample_data()


ps_prosper_rf_filter_raref <- ps_prosper_rf_filter %>%
  rarefy_even_depth(.,sample.size = 2000)

ps_prosper_rf_comp <- ps_prosper_rf_filter %>%
  prune_taxa(taxa_names(ps_prosper_rf_filter_raref),.) %>%
  microbiome::transform(transform = "compositional")

young_taxa <- ps_prosper_rf_filter_raref %>%
  subset_samples(day_of_life < 100) %>%
  prevalence(., detection = 10)
young_taxa <- names(young_taxa[young_taxa < 0.2])
old_taxa <- ps_prosper_rf_filter_raref %>%
  subset_samples(day_of_life > 100) %>%
  prevalence(., detection = 10)
old_taxa <- names(old_taxa[old_taxa < 0.2])
remove_taxa <- intersect(old_taxa, young_taxa) %>% unique

ps_prosper_rf_filter_raref <- ps_prosper_rf_filter_raref %>%
  merge_taxa2_fixed(.,taxa = remove_taxa, name = "Other")

ps_prosper_rf_comp <- ps_prosper_rf_comp %>%
  merge_taxa2_fixed(.,taxa = remove_taxa, name = "Other")

ps_prosper_rf_filter <- ps_prosper_rf_filter %>%
  prune_taxa(unique(c(taxa_names(ps_prosper_rf_filter_raref), remove_taxa)),.) %>%
  merge_taxa2_fixed(.,taxa = remove_taxa, name = "Other")


diversities <- estimate_richness(ps_prosper_rf_filter_raref, measures = c("Observed", "Shannon"))
sample_data(ps_prosper_rf_filter_raref)$Observed <- diversities$Observed
sample_data(ps_prosper_rf_filter_raref)$Shannon <- diversities$Shannon

mdat <- sample_data(ps_prosper_rf_filter_raref) %>% data.frame


# DMM clustering ###############################
# check dirichelet for the subset:
count_gen <- ps_prosper_rf_filter_raref %>%
  otu_table() %>%
  t() 
if(dmm_step){
  fit_gen <- mclapply(c(1:10), dmn, count=count_gen, verbose=TRUE, mc.cores = 10) 
  saveRDS(fit_gen, "/fast/AG_Forslund/rob/PROSPER/immune_paper_data/data/genus_dirichlet_model_immune_subset.rds")
} else {fit_gen <- readRDS("/fast/AG_Forslund/rob/PROSPER/immune_paper_data/data/genus_dirichlet_model_immune_subset.rds")}


lplc <- sapply(fit_gen, laplace)
aic <- sapply(fit_gen, AIC)
bic <- sapply(fit_gen, BIC)
plot(lplc, type="b", xlab="Number of Dirichlet Components",ylab="Model Fit")
lines(aic, type="b", lty = 2)
lines(bic, type="b", lty = 3)
(best <- fit_gen[[which.min(lplc)]])
# (best <- fit_gen[[3]])
mixturewt(best)

sample_data(ps_prosper_rf_filter_raref)$dmm_mixture_immune <- as.factor(mixture(best, assign = T))
mdat <- sample_data(ps_prosper_rf_filter_raref) %>% data.frame


## check in PCoA: ############
dist.bray_genus <- phyloseq::distance(ps_prosper_rf_filter_raref, method = "bray")
pcoa_bray_genus <- ape::pcoa(dist.bray_genus)
fit_adonis_genus <- adonis2(dist.bray_genus ~ day_of_life + dmm_mixture_immune + MHH_PatID,
                            data = mdat, by="margin", parallel = 10)

save(ps_prosper_rf_filter_raref, pcoa_bray_genus,
     file = "/fast/AG_Forslund/rob/PROSPER/immune_paper_data/data/for_paper_plots_dmm_pcoa.RData")
