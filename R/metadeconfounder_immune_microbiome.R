# includes analyses on all samples of individuals that have any kind of immune data
# for samples with paired immune data see script _small
library(tidyverse)
library(ggplot2)
library(readxl)
library(phyloseq)
library(microbiome)
library(vegan)
library(metadeconfoundR)
library(r2glmm)
library(furrr)
library(purrr)
library(lme4)
library(gridExtra)
library(gamlss)
library(fastDummies)
library(caret)
library(glue)
library(ggtext)
library(openxlsx)

source("/fast/AG_Forslund/rob/PROSPER/immune_paper_data/R/functions.R")
setwd("/fast/AG_Forslund/rob/PROSPER/immune_paper_data/R/")
set.seed(123)

# Switches ##############
metadec_step <- F
metadec_tp_1_step <- F
metadec_tp_2_step <- F
metadec_tp_3_step <- F
metadec_gamlss_step <- F
metadec_mediation_step <- F
# devtools::install_github("TillBirkner/metadeconfoundR@develop")


remove_cols <- c("cutoff", "biosample_LT_chronological", "MHH_unique_ID", "MHH_sampling_ID",
                 "sample_no", "X.SampleID", "SampleNr", "Description", "data_set",
                 "PRIMAL_sample_ID", "RESIST_sample_ID", "sample_name_stool",
                 "matching_mo.RNA.seq..MHH_sample_name.", "matching_pharyngeal_swab..MHH_sample_name.",
                 "matching_chip..MHH_sample_name.", "matching_chip_ID..chip_ID.", "birth_date",
                 "fecal_sample_ID_1", "PRIMAL_sample_ID_2", "PRIMAL_subject_ID", "group_median_prosper",
                 "group_sd_prosper", "group_median", "group_sd", "corrected_group_median", "corrected_group_sd",
                 "corrected_age_cat", "predicted_age",
                 "preterm_stat", "percentile_birth_weight", "MAZ_prosper", "MAZ_corrected", 
                 "corrected_age", "sampling_date", "time_point_overall",
                 "age_cat_2", "sample_sum", "Aktuelle_Menge_MM_ml_am_Probentag",
                 "Anteil_MM_Probeentnahmezeitpunkt", "MAZ", "Mittelwert_Anteil_MM_d1_30",
                 "MM_S100A8_A9_Konzentration", "Serum_S100A8_A9_Konzentration",
                 "Stuhl_S100A8_A9_Konzentration", "time_between_sampling_and_sepsis_lons", 
                 "storage_in_anything", "time_in_freezer", "time_spent_in_lockdown", "lockdown_binary",
                 "sepsis_eons_lons_binary", "metacluster_kilian", "group_lockdown_sum_LD_Post.expo",
                 "group_lockdown_sum_LD_Pre.expo", 
                 "group_lockdown_sum_LD_w.o.expo", 
                 "group_lockdown_sum_Non.LD_Pre.LD")

# load data:
ps_prosper_rf_filter <- readRDS(file = "/fast/AG_Forslund/rob/PROSPER/dada2/phyloseq_combined_rf_filter_genus.rds") %>%
  microbiome::aggregate_taxa(., level = "genus")


immune_data_rel <- read_xlsx("/fast/AG_Forslund/rob/PROSPER/immune_data/frequencies_of_pbmc.xlsx") %>%
  dplyr::select(-id, -timepoint) %>%
  column_to_rownames("chip") %>%
  as.matrix()

# data intersection:
# change rownames of immune data
immune_data_rel <- immune_data_rel[rownames(immune_data_rel) %in% ps_prosper_rf_filter@sam_data$matching_chip_ID..chip_ID.,]

id_dict <- ps_prosper_rf_filter@sam_data %>% data.frame() %>%
  select(matching_chip_ID..chip_ID.) %>%
  rownames_to_column() %>%
  filter(!duplicated(matching_chip_ID..chip_ID.)) %>%
  filter(matching_chip_ID..chip_ID. %in% rownames(immune_data_rel)) %>%
  column_to_rownames("matching_chip_ID..chip_ID.")
rownames(immune_data_rel) <- id_dict[rownames(immune_data_rel),]
immune_data_rel_df <- data.frame(immune_data_rel)

# get intersection of samples
keep_individuals <- ps_prosper_rf_filter@sam_data %>% data.frame() %>% 
  rownames_to_column() %>%
  filter(rowname %in% rownames(immune_data_rel)) %>%
  pull(MHH_PatID) %>%
  unique()

ps_prosper_rf_filter <- ps_prosper_rf_filter %>%
  subset_samples(MHH_PatID %in% keep_individuals)

# filter data
ps_prosper_rf_filter_raref <- ps_prosper_rf_filter %>%
  rarefy_even_depth(.,sample.size = 2000)
ps_prosper_rf_filter_raref_family <- ps_prosper_rf_filter %>%
  aggregate_taxa(., level = "family") %>%
  rarefy_even_depth(.,sample.size = 2000)


# age specific prevalence filter
young_taxa <- ps_prosper_rf_filter_raref %>%
  subset_samples(day_of_life < 100) %>%
  prevalence(., detection = 10)
young_taxa <- names(young_taxa[young_taxa < 0.1])
old_taxa <- ps_prosper_rf_filter_raref %>%
  subset_samples(day_of_life > 100) %>%
  prevalence(., detection = 10)
old_taxa <- names(old_taxa[old_taxa < 0.1])
remove_taxa <- intersect(old_taxa, young_taxa) %>% unique
# family level
young_families <- ps_prosper_rf_filter_raref_family %>%
  subset_samples(day_of_life < 100) %>%
  prevalence(., detection = 10)
young_families <- names(young_families[young_families < 0.2])
old_families <- ps_prosper_rf_filter_raref_family %>%
  subset_samples(day_of_life > 100) %>%
  prevalence(., detection = 10)
old_families <- names(old_families[old_families < 0.2])
remove_families <- intersect(old_families, young_families) %>% unique

ps_prosper_rf_filter_raref <- ps_prosper_rf_filter_raref %>%
  merge_taxa2_fixed(.,taxa = remove_taxa, name = "Other")
ps_prosper_rf_filter_raref_family <- ps_prosper_rf_filter_raref_family %>%
  merge_taxa2_fixed(.,taxa = remove_families, name = "Other")


# alpha diversity
diversities <- estimate_richness(ps_prosper_rf_filter_raref, measures = c("Observed", "Shannon"))
sample_data(ps_prosper_rf_filter_raref)$`Observed` <- diversities$Observed
sample_data(ps_prosper_rf_filter_raref)$`Shannon` <- diversities$Shannon

mdat_metadec <- sample_data(ps_prosper_rf_filter_raref) %>%
  data.frame()

# load nicer metadata names
display_names <- read_xlsx("/fast/AG_Forslund/rob/PROSPER/metadata/metadata_rename.xlsx")
display_names <- left_join(data.frame(old_name = colnames(mdat_metadec)), display_names) %>%
  mutate(new_name = ifelse(is.na(new_name), yes = old_name, no = new_name))
# split between metadata, immune data and microbiome data

# create mb-feature table with a-div and clusters:
# add other microbiome features, convert factors to dummy variables
mb_dat_in_metad <- c("Shannon", "Observed")

mb_special_data <- mdat_metadec %>% 
  rownames_to_column("rowname") %>%
  select(mb_dat_in_metad, rowname) %>%
  column_to_rownames("rowname")



mb_otu_matrix_df <- ps_prosper_rf_filter_raref %>% otu_table() %>% as.data.frame %>% t %>% as.data.frame %>%
  bind_cols(mb_special_data)
mb_otu_matrix_df_family <- ps_prosper_rf_filter_raref_family %>% otu_table() %>% as.data.frame %>% t %>% as.data.frame %>%
  bind_cols(mb_special_data)


mdat_metadec <- mdat_metadec %>% select(-c("Shannon", "Observed", 
                                           "two_tp_metaclusters",
                                           "two_tp_metaclusters_young"))


# metadeconfounder runs ###########################################
mdat_metadec <- mdat_metadec %>%
  dplyr::select(-remove_cols, -c("age_cat", "MHH_mother_ID")) %>%
  dplyr::select(!starts_with("dmm")) %>%
  select(-(nearZeroVar(., freqCut = nrow(.)/10, uniqueCut = 250/nrow(.)))) %>% # only keep things with at least two variables and more than 5 in each binary category
  dplyr::select(sex_male,  everything())

if(metadec_step){
  # mibi vs metadata
  all(rownames(mdat_metadec) == rownames(mb_otu_matrix_df))
  metad_mb <- metadeconfoundR::MetaDeconfound(featureMat = mb_otu_matrix_df,
                                              metaMat = mdat_metadec %>% select(-starts_with("ITJ")),
                                              adjustLevel = 2,
                                              returnLong = T,
                                              QCutoff = 0.05,
                                              DCutoff = 0.2,
                                              randomVar = c("cohort", "MHH_PatID"),
                                              nnodes = 10)
  metad_mb_ITJ <- metadeconfoundR::MetaDeconfound(featureMat = mb_otu_matrix_df,
                                              metaMat = mdat_metadec %>% select(starts_with("ITJ"), day_of_life),
                                              adjustLevel = 2,
                                              returnLong = T,
                                              QCutoff = 0.05,
                                              DCutoff = 0.2,
                                              randomVar = c("cohort", "MHH_PatID"),
                                              nnodes = 10)
  metad_mb_family <- metadeconfoundR::MetaDeconfound(featureMat = mb_otu_matrix_df_family,
                                                     metaMat = mdat_metadec,
                                                     # fixedVar = c("day_of_life"),
                                                     adjustLevel = 2, 
                                                     returnLong = T, 
                                                     QCutoff = 0.05,
                                                     DCutoff = 0.2, 
                                                     randomVar = c("cohort", "MHH_PatID"),
                                                     nnodes = 10)
  

    save(metad_mb, mdat_metadec, metad_mb_ITJ, metad_mb_family,
       file = "/fast/AG_Forslund/rob/PROSPER/immune_paper_data/data/metadec_mibi_immune_metadat.RData")
} else {
  load(file = "/fast/AG_Forslund/rob/PROSPER/immune_paper_data/data/metadec_mibi_immune_metadat.RData")
}


# time point specific runs #####################################################
## TP 1 ##############################
mdat_metadec_tp1 <- mdat_metadec %>% filter(day_of_life < 19) # 19-47, 415-490
mb_otu_matrix_df_tp1 <- mb_otu_matrix_df[rownames(mdat_metadec_tp1),]

if(metadec_tp_1_step){
  # mibi vs metadata
  metad_mb_tp1 <- metadeconfoundR::MetaDeconfound(featureMat = mb_otu_matrix_df_tp1,
                                              metaMat = mdat_metadec_tp1 %>% select(-starts_with("ITJ")),
                                              # fixedVar = c("day_of_life"),
                                              adjustLevel = 2,
                                              returnLong = T,
                                              QCutoff = 0.05,
                                              DCutoff = 0.2,
                                              randomVar = c("cohort", "MHH_PatID"),
                                              nnodes = 10)
  metad_mb_ITJ_tp1 <- metadeconfoundR::MetaDeconfound(featureMat = mb_otu_matrix_df_tp1,
                                                  metaMat = mdat_metadec_tp1 %>% select(starts_with("ITJ")),
                                                  # fixedVar = c("day_of_life"),
                                                  adjustLevel = 2,
                                                  returnLong = T,
                                                  QCutoff = 0.05,
                                                  DCutoff = 0.2,
                                                  randomVar = c("cohort", "MHH_PatID"),
                                                  nnodes = 10)
  
  save(metad_mb_tp1, metad_mb_ITJ_tp1,
       file = "/fast/AG_Forslund/rob/PROSPER/immune_paper_data/data/metadec_mibi_immune_metadat_tp1.RData")
} else {
  load(file = "/fast/AG_Forslund/rob/PROSPER/immune_paper_data/data/metadec_mibi_immune_metadat_tp1.RData")
}


## TP 2 ##############################
mdat_metadec_tp2 <- mdat_metadec %>% filter(day_of_life >= 19 & day_of_life <= 47) # 19-47, 415-490
mb_otu_matrix_df_tp2 <- mb_otu_matrix_df[rownames(mdat_metadec_tp2),]

if(metadec_tp_2_step){
  # mibi vs metadata
  metad_mb_tp2 <- metadeconfoundR::MetaDeconfound(featureMat = mb_otu_matrix_df_tp2,
                                                  metaMat = mdat_metadec_tp2 %>% select(-starts_with("ITJ")),
                                                  # fixedVar = c("day_of_life"),
                                                  adjustLevel = 2,
                                                  returnLong = T,
                                                  QCutoff = 0.05,
                                                  DCutoff = 0.2,
                                                  randomVar = c("cohort"),
                                                  nnodes = 10)
  metad_mb_ITJ_tp2 <- metadeconfoundR::MetaDeconfound(featureMat = mb_otu_matrix_df_tp2,
                                                  metaMat = mdat_metadec_tp2 %>% select(starts_with("ITJ")),
                                                  # fixedVar = c("day_of_life"),
                                                  adjustLevel = 2,
                                                  returnLong = T,
                                                  QCutoff = 0.05,
                                                  DCutoff = 0.2,
                                                  randomVar = c("cohort"),
                                                  nnodes = 10)
  
  save(metad_mb_tp2, metad_mb_ITJ_tp2,
       file = "/fast/AG_Forslund/rob/PROSPER/immune_paper_data/data/metadec_mibi_immune_metadat_tp2.RData")
} else {
  load(file = "/fast/AG_Forslund/rob/PROSPER/immune_paper_data/data/metadec_mibi_immune_metadat_tp2.RData")
}

## TP 3 ##############################
mdat_metadec_tp3 <- mdat_metadec %>% filter(day_of_life >= 255 & day_of_life <= 490) # 19-47, 415-490
mb_otu_matrix_df_tp3 <- mb_otu_matrix_df[rownames(mdat_metadec_tp3),]

if(metadec_tp_3_step){
  # mibi vs metadata
  metad_mb_tp3 <- metadeconfoundR::MetaDeconfound(featureMat = mb_otu_matrix_df_tp3 %>% select(-starts_with("dmm")),
                                                  metaMat = mdat_metadec_tp3 %>% select(-starts_with("ITJ")),
                                                  # fixedVar = c("day_of_life"),
                                                  adjustLevel = 2,
                                                  returnLong = T,
                                                  QCutoff = 0.05,
                                                  DCutoff = 0.2,
                                                  randomVar = c("cohort", "MHH_PatID"),
                                                  nnodes = 10)
  metad_mb_ITJ_tp3 <- metadeconfoundR::MetaDeconfound(featureMat = mb_otu_matrix_df_tp3 %>% select(-starts_with("dmm")),
                                                  metaMat = mdat_metadec_tp3 %>% select(starts_with("ITJ")),
                                                  # fixedVar = c("day_of_life"),
                                                  adjustLevel = 2,
                                                  returnLong = T,
                                                  QCutoff = 0.05,
                                                  DCutoff = 0.2,
                                                  randomVar = c("cohort", "MHH_PatID"),
                                                  nnodes = 10)
  
  save(metad_mb_tp3, metad_mb_ITJ_tp3,
       file = "/fast/AG_Forslund/rob/PROSPER/immune_paper_data/data/metadec_mibi_immune_metadat_tp3.RData")
} else {
  load(file = "/fast/AG_Forslund/rob/PROSPER/immune_paper_data/data/metadec_mibi_immune_metadat_tp3.RData")
}


# write excel files: #################
rename_variables <- function(df) {
  df %>% left_join(., display_names, by = c("metaVariable" = "old_name")) %>%
    mutate(metaVariable = as.character(metaVariable),
           metaVariable = ifelse(is.na(new_name), yes = metaVariable, no = new_name),
           metaVariable = factor(metaVariable),
           feature = gsub("unclassified_", "Unclassified ", feature),
           feature = gsub("Observed", "Observed taxa", feature),
           feature = gsub("Shannon", "Shannon diversity", feature),
           feature = gsub("dmm_mixture_immune_", "DMM cluster ", feature),
           Confounders = ifelse(grepl("C: ", status), yes = gsub("C: ", "", status), no = ""),
           status = ifelse(grepl("C: ", status), yes = "C", no = status),
           `Microbiome feature` = feature,
           Metavariable = metaVariable,
           `p-value` = Ps,
           `FDR-corrected p-value` = Qs,
           `Effect size` = Ds,
           `Confounding status` = status) %>%
    select(c("Microbiome feature", "Metavariable", "p-value", "FDR-corrected p-value", 
             "Effect size", "Confounding status", "Confounders"))}

# Create a new workbook and add a sheet
wb <- createWorkbook()
addWorksheet(wb, "Metadec_all_TPs_genus")
writeData(wb, "Metadec_all_TPs_genus", metad_mb %>% rename_variables())
addWorksheet(wb, "Metadec_all_TPs_genus_ITJs")
writeData(wb, "Metadec_all_TPs_genus_ITJs", metad_mb_ITJ %>% rename_variables())
addWorksheet(wb, "Metadec_all_TPs_family")
writeData(wb, "Metadec_all_TPs_family", metad_mb_family %>% rename_variables())
addWorksheet(wb, "Metadec_TP_1")
writeData(wb, "Metadec_TP_1", metad_mb_tp1 %>% rename_variables())
addWorksheet(wb, "Metadec_TP_1_ITJs")
writeData(wb, "Metadec_TP_1_ITJs", metad_mb_ITJ_tp1 %>% rename_variables())
addWorksheet(wb, "Metadec_TP_2")
writeData(wb, "Metadec_TP_2", metad_mb_tp2 %>% rename_variables())
addWorksheet(wb, "Metadec_TP_2_ITJs")
writeData(wb, "Metadec_TP_2_ITJs", metad_mb_ITJ_tp2 %>% rename_variables())
addWorksheet(wb, "Metadec_TP_3")
writeData(wb, "Metadec_TP_3", metad_mb_tp3 %>% rename_variables())
addWorksheet(wb, "Metadec_TP_3_ITJs")
writeData(wb, "Metadec_TP_3_ITJs", metad_mb_ITJ_tp3 %>% rename_variables())
saveWorkbook(wb, "/fast/AG_Forslund/rob/PROSPER/immune_paper_data/data/metadeconfounder_res.xlsx", overwrite = TRUE)


wb_mb <- createWorkbook()
addWorksheet(wb_mb, "Genus_level")
writeData(wb_mb, "Genus_level", mb_otu_matrix_df %>% rownames_to_column("sample_ID"))
addWorksheet(wb_mb, "Family_level")
writeData(wb_mb, "Family_level", mb_otu_matrix_df_family %>% rownames_to_column("sample_ID"))
addWorksheet(wb_mb, "Metadata")
writeData(wb_mb, "Metadata", mdat_metadec %>% select(day_of_life, GA_days, sex_male) %>% 
            rownames_to_column("sample_ID") %>%
            mutate(Age = day_of_life,
                   `Sex (m)` = sex_male, 
                   GA = GA_days,
                   time_point = case_when(Age < 19 ~ 1,
                                       Age >= 19 & Age <= 47 ~ 2,
                                       Age >= 255 ~ 3),
                   .keep = "unused"))
saveWorkbook(wb_mb, "/fast/AG_Forslund/rob/PROSPER/immune_paper_data/data/abundances.xlsx", overwrite = TRUE)