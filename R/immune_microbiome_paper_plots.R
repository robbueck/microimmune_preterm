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
library(ggprism)
library(ggtext)
source("/fast/AG_Forslund/rob/PROSPER/immune_paper_data/R/functions.R")
setwd("/fast/AG_Forslund/rob/PROSPER/immune_paper_data/R/")

load(file = "/fast/AG_Forslund/rob/PROSPER/immune_paper_data/data/metadec_mibi_immune_metadat.RData")
# metad_mb, metad_mb_bin, metad_immune, metad_mb_family,

# load clear names for variables
display_names <- read_xlsx("/fast/AG_Forslund/rob/PROSPER/metadata/metadata_rename.xlsx")
display_names <- left_join(data.frame(old_name = colnames(mdat_metadec)), display_names) %>%
  mutate(new_name = ifelse(is.na(new_name), yes = old_name, no = new_name))


# Fig 6: #########################
library(RColorBrewer)
library(circlize)
color_palette <- c("Microbiome" = "#30123BFF", 
                   "Clinical factors" = "#FB8022FF", 
                   "Immune trajectory" = "#28BBECFF")
metad_combined <- bind_rows(metad_mb %>% mutate(feature_type = "Microbiome", 
                                                feature_group = ifelse(feature %in% c("Shannon", "Observed"), 
                                                                       yes = "Microbiome_extra", no = "Microbiome"), 
                                                metavar_type = "Clinical factors"),
                            metad_mb_ITJ %>% mutate(feature_type = "Microbiome",
                                                    feature_group = "Microbiome", 
                                                    metavar_type = "Immune trajectory")) %>%
  left_join(., display_names, by = c("metaVariable" = "old_name")) %>%
  mutate(metaVariable = as.character(metaVariable),
         metaVariable = ifelse(is.na(new_name), yes = metaVariable, no = new_name),
         metaVariable = factor(metaVariable),
         feature = gsub("unclassified_", "Unclassified\n", feature),
         feature = gsub("Observed", "Observed taxa", feature),
         feature = gsub("Shannon", "Shannon diversity", feature),
         feature = gsub("dmm_mixture_immune_", "DMM cluster ", feature),
         # feature = if_else(grepl("\n| ", feature) |
         #                     feature == "Other",
         #                   true = glue("{feature}"),
         #                   false = glue("<i>{feature}</i>")),
         feature = factor(feature, levels = unique(c("Shannon diversity", "Observed taxa", feature))))%>%
  filter(#!grepl("dmm", feature),
         #!grepl("dmm", metaVariable),
         !grepl("group_lockdown", metaVariable),
         !metaVariable %in% c("storage_in_stabilisator", "immediate_freezing"))

df_sig <- metad_combined %>% 
  filter(Qs < 0.05, 
         abs(Ds) > 0.2, !is.na(Ds), 
         status != "NS",
         !is.infinite(Ds),
         metaVariable != "day_of_life")

# Plot positive associations (Ds > 0)
df_pos <- df_sig %>% filter(Ds > 0.2)
df_neg <- df_sig %>% filter(Ds < 0)


if(nrow(df_pos) > 0){
  pdf("/fast/AG_Forslund/rob/PROSPER/immune_paper_data/figs/paper_fig_6_a_mdc_res_immune_vs_mb_circ_pos.pdf", width = 7, height = 8.5)
  plot_circular_associations(df_subset = df_pos, title_text = NULL, type_colors = color_palette)
  dev.off()
}
# Plot negative associations (Ds < 0)
if(nrow(df_neg) > 0){
  pdf("/fast/AG_Forslund/rob/PROSPER/immune_paper_data/figs/paper_fig_6_b_mdc_res_immune_vs_mb_circ_neg.pdf", width = 7, height = 8.5)
  plot_circular_associations(df_neg, title_text = NULL, type_colors = color_palette)
  dev.off()
}

# Supplementary Figure S11 A: #####################
load("/fast/AG_Forslund/rob/PROSPER/immune_paper_data/data/for_paper_plots_dmm_pcoa.RData")

## Genera vs time points #################################
library(RColorBrewer)
taxa_names(ps_prosper_rf_filter_raref) <- gsub("unclassified_", "Unclassified ", taxa_names(ps_prosper_rf_filter_raref))
unique_taxa <- unique(taxa_names(ps_prosper_rf_filter_raref))
color_palette <- setNames(viridisLite::turbo(length(unique_taxa)), unique_taxa) 

# for time points: 
dmm_cluster_comps_immune_age <- ps_prosper_rf_filter_raref %>%
  subset_samples(., time_point_immune %in% c(1, 2, 3)) %>%
  microbiome::transform(transform = "compositional") %>%
  plot_composition(sample.sort = "day_of_life", otu.sort = "abundance", average_by = "time_point_immune") +
  scale_fill_manual("Genus", values = color_palette) +
  xlab("Time point") +
  ylab("Relative abundance") +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size= 15))
ggsave("/fast/AG_Forslund/rob/PROSPER/immune_paper_data/figs/paper_fig_s11_a_immune_time_point_composition.pdf",
       plot = dmm_cluster_comps_immune_age, scale = 1)
ggsave("/fast/AG_Forslund/rob/PROSPER/immune_paper_data/figs/paper_fig_s11_a_immune_time_point_composition.tiff",
       plot = dmm_cluster_comps_immune_age, scale = 1, dpi = 900)

## Generas vs DMMs ######################
dmm_cluster_comps_immune <- ps_prosper_rf_filter_raref %>%
  subset_samples(., time_point_immune %in% c(1, 2, 3)) %>%
  microbiome::transform(transform = "compositional") %>%
  plot_composition(sample.sort = "day_of_life", otu.sort = "abundance", average_by = "dmm_mixture_immune") +
  scale_fill_manual("Genus", values = color_palette) +
  xlab("DMM-cluster") +
  ylab("Relative abundance") +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size= 15))

ggsave("/fast/AG_Forslund/rob/PROSPER/immune_paper_data/figs/paper_fig_s11_c_immune_dmm_clusters_composition.pdf",
       plot = dmm_cluster_comps_immune, scale = 1)
ggsave("/fast/AG_Forslund/rob/PROSPER/immune_paper_data/figs/paper_fig_s11_c_immune_dmm_clusters_composition.tiff",
       plot = dmm_cluster_comps_immune, scale = 1, dpi = 900)


## PCoAs ###############################
(p1 <- plot_pcoa(ps_prosper_rf_filter_raref, pcoa_bray_genus, 
                 color="dmm_mixture_immune",
                 ellipses = T) +
   scale_color_manual("DMM-cluster", values = viridisLite::turbo(3),
                      guide = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1))) + 
   theme_bw() +
   theme(panel.border = element_blank(),
         axis.text.x = element_text(size = 13),
         axis.text.y = element_text(size = 13),
         axis.title.x = element_text(size = 15),
         axis.title.y = element_text(size= 15)))
ggsave("/fast/AG_Forslund/rob/PROSPER/immune_paper_data/figs/paper_fig_s11_e_dmm_clusters_pcoa.pdf",
       plot = p1, scale = 1, width = 8, height = 6.5)
ggsave("/fast/AG_Forslund/rob/PROSPER/immune_paper_data/figs/paper_fig_s11_e_dmm_clusters_pcoa.tiff",
       plot = p1, scale = 1, width = 8, height = 6.5, dpi = 900)


(p2 <- plot_pcoa(ps_prosper_rf_filter_raref, pcoa_bray_genus, 
                 color="time_point_immune",
                 ellipses = T) +
    scale_color_manual("Time point", values = viridisLite::turbo(3),
                       guide = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1))) + 
    theme_bw() +
    theme(panel.border = element_blank(),
          axis.text.x = element_text(size = 13),
          axis.text.y = element_text(size = 13),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size= 15)))
ggsave("/fast/AG_Forslund/rob/PROSPER/immune_paper_data/figs/paper_fig_s11_d_immune_age_pcoa.pdf",
       plot = p2, scale = 1, width = 8, height = 6.5)
ggsave("/fast/AG_Forslund/rob/PROSPER/immune_paper_data/figs/paper_fig_s11_d_immune_age_pcoa.tiff",
       plot = p2, scale = 1, width = 8, height = 6.5, dpi = 900)

## alpha diversity ##########################
df_p_val_a_div <- sample_data(ps_prosper_rf_filter_raref) %>% data.frame() %>%
  rstatix::wilcox_test(Shannon ~ time_point_immune, detailed = T) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj") %>%
  rstatix::add_xy_position()


(a_div_plot <- sample_data(ps_prosper_rf_filter_raref) %>% data.frame() %>%
  ggplot(aes(x = time_point_immune, y = Shannon)) +
  geom_boxplot()+
  xlab("Time point") +
  ylab("Shannon diversity") +
    add_pvalue(df_p_val_a_div, 
               label = "{p.adj.signif}",
               step.increase = 0.1,
               tip.length = 0.01) +
   theme_bw() +
   theme(panel.border = element_blank(),
         axis.text.x = element_text(size = 13),
         axis.text.y = element_text(size = 13),
         axis.title.x = element_text(size = 15),
         axis.title.y = element_text(size= 15)))
ggsave("/fast/AG_Forslund/rob/PROSPER/immune_paper_data/figs/paper_fig_s11_b_a_div_vs_age.pdf",
       plot = a_div_plot, width = 4, height = 3)
ggsave("/fast/AG_Forslund/rob/PROSPER/immune_paper_data/figs/paper_fig_s11_b_a_div_vs_age.tiff",
       plot = a_div_plot, width = 4, height = 3, dpi = 900)



# Supplementary Figure S12: #####################
# Microbiome vs clinical factors
## Overall ##############################
(mibi_vs_clinical_all <- BuildHeatmap_2(bind_rows(metad_mb) %>% 
                                filter(!metaVariable %in% c("ITJ1", "ITJ2", "ITJ3"),
                                       !metaVariable %in% c("storage_in_stabilisator", "immediate_freezing")) %>%
                                mutate(metaVariable = factor(metaVariable, 
                                                             levels = unique(c("day_of_life", "Vaginal_delivery", "antibiotics_neonatal",
                                                                               levels(metaVariable))))) %>%
                                mutate(#status = ifelse(grepl("C:", status), yes = "OK_sd", no = status),
                                       feature = gsub("unclassified_", "Unclassified ", feature),
                                       feature = gsub("Observed", "Observed taxa", feature),
                                       feature = gsub("Shannon", "Shannon diversity", feature),
                                       feature = gsub("dmm_mixture_immune_", "DMM cluster ", feature),
                                       feature = if_else(grepl("_| ", feature) |
                                                           feature == "Other",
                                                         true = glue("{feature}"),
                                                         false = glue("<i>{feature}</i>")),
                                       feature = factor(feature, levels = unique(c("Shannon diversity", "Observed taxa", feature)))), 
                              keepMeta = c("sepsis_eons_lons_binary", "Vaginal_delivery", "antibiotics_neonatal"),
                              metaVariableNames = display_names,
                              coloring = 1,
                              reOrder = "none",
                              q_cutoff = 0.05, d_cutoff = 0.2, starSize = 6) +#, starNudge_y = -0.5) +
  ggtitle("Clinical factor associations across the first year of life") +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_markdown(size = 14),
        axis.title = element_blank(),
        plot.subtitle = element_blank(),
        plot.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        strip.text.x = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)))
ggsave("/fast/AG_Forslund/rob/PROSPER/immune_paper_data/figs/paper_fig_s12_a_mdc_res_microbiome_vs_metadata.pdf", 
       width = 11, height = 7, plot = mibi_vs_clinical_all)
ggsave("/fast/AG_Forslund/rob/PROSPER/immune_paper_data/figs/paper_fig_s12_a_mdc_res_microbiome_vs_metadata.tiff", 
       width = 11, height = 7, plot = mibi_vs_clinical_all, dpi = 900)



## Per time point ##############################
load(file = "/fast/AG_Forslund/rob/PROSPER/immune_paper_data/data/metadec_mibi_immune_metadat_tp1.RData")
load(file = "/fast/AG_Forslund/rob/PROSPER/immune_paper_data/data/metadec_mibi_immune_metadat_tp2.RData")
load(file = "/fast/AG_Forslund/rob/PROSPER/immune_paper_data/data/metadec_mibi_immune_metadat_tp3.RData")

(mibi_vs_clinical_per_tp <- BuildHeatmap_2(bind_rows(bind_rows(metad_mb_tp1) %>% mutate(groupingVar = "TP1") %>% 
                                                       group_by(metaVariable) %>% filter(any(status != "NS")),
                                                     bind_rows(metad_mb_tp2) %>% mutate(groupingVar = "TP2") %>%
                                                       group_by(metaVariable) %>% filter(any(status != "NS")),
                                                     bind_rows(metad_mb_tp3) %>% mutate(groupingVar = "TP3") %>%
                                                       group_by(metaVariable) %>% filter(any(status != "NS"))) %>%
                                             filter(!metaVariable %in% c("ITJ1", "ITJ2", "ITJ3"),
                                                    !metaVariable %in% c("storage_in_stabilisator", "immediate_freezing")) %>%
                                             mutate(metaVariable = factor(metaVariable, 
                                                                          levels = unique(c("day_of_life", "Vaginal_delivery", "antibiotics_neonatal",
                                                                                            levels(metaVariable))))) %>%
                                             mutate(#status = ifelse(grepl("C:", status), yes = "OK_sd", no = status),
                                               feature = gsub("unclassified_", "Unclassified ", feature),
                                               feature = gsub("Observed", "Observed taxa", feature),
                                               feature = gsub("Shannon", "Shannon diversity", feature),
                                               feature = gsub("dmm_mixture_immune_", "DMM cluster ", feature),
                                               feature = if_else(grepl("_| ", feature) |
                                                                   feature == "Other",
                                                                 true = glue("{feature}"),
                                                                 false = glue("<i>{feature}</i>")),
                                               feature = factor(feature, levels = unique(c("Shannon diversity", "Observed taxa", feature)))), 
                                           keepMeta = c("sepsis_eons_lons_binary", "Vaginal_delivery", "antibiotics_neonatal"),
                                           metaVariableNames = display_names,
                                           coloring = 1,
                                           reOrder = "none",
                                           q_cutoff = 0.05, d_cutoff = 0.2, starSize = 4) +
   theme(axis.text.x = element_text(size = 14),
         axis.text.y = element_markdown(size = 14),
         axis.title = element_blank(),
         plot.subtitle = element_blank(),
         plot.title = element_blank(),
         axis.ticks = element_blank(),
         axis.line = element_blank(),
         strip.text.x = element_text(size = 14),
         legend.text = element_text(size = 12),
         legend.title = element_text(size = 12)))
ggsave("/fast/AG_Forslund/rob/PROSPER/immune_paper_data/figs/paper_fig_s12_b_mdc_res_microbiome_vs_metadata_per_tp.pdf", 
       width = 13.2, height = 6.3, plot = mibi_vs_clinical_per_tp)
ggsave("/fast/AG_Forslund/rob/PROSPER/immune_paper_data/figs/paper_fig_s12_b_mdc_res_microbiome_vs_metadata_per_tp.tiff", 
       width = 13.2, height = 6.3, plot = mibi_vs_clinical_per_tp, dpi = 900)

# Supplementary Figure S13: #####################
# Microbiome vs ITJs
## Overall ##############################
(mibi_vs_ITJs_all <- BuildHeatmap_2(bind_rows(metad_mb_ITJ) %>% 
                                      filter(metaVariable %in% c("ITJ1", "ITJ2", "ITJ3", "day_of_life")) %>%
                                      mutate(#status = ifelse(grepl("C:", status), yes = "OK_sd", no = status),
                                        feature = gsub("unclassified_", "Unclassified ", feature),
                                        feature = gsub("Observed", "Observed taxa", feature),
                                        feature = gsub("Shannon", "Shannon diversity", feature),
                                        feature = gsub("dmm_mixture_immune_", "DMM cluster ", feature),
                                        feature = if_else(grepl("_| ", feature) |
                                                            feature == "Other",
                                                          true = glue("{feature}"),
                                                          false = glue("<i>{feature}</i>")),
                                        feature = factor(feature, levels = unique(c("Shannon diversity", "Observed taxa", feature)))), 
                                    q_cutoff = 0.05, d_cutoff = 0.2, starSize = 5,# starNudge_y = -0.5,
                                    metaVariableNames = display_names,
                                    coloring = 1,
                                    reOrder = "none",
                                    keepMeta = c("ITJ1", "ITJ2", "ITJ3")) +
   theme(axis.text.x = element_text(size = 14),
         axis.text.y = element_markdown(size = 14),
         axis.title = element_blank(),
         axis.ticks = element_blank(),
         axis.line = element_blank(),
         plot.subtitle = element_blank(),
         plot.title = element_blank(),
         strip.text.x = element_blank(),
         legend.text = element_text(size = 12),
         legend.title = element_text(size = 12)))
ggsave("/fast/AG_Forslund/rob/PROSPER/immune_paper_data/figs/paper_fig_s13_a_mdc_res_microbiome_vs_ITJs.pdf", 
       width = 6, height = 5, plot = mibi_vs_ITJs_all)
ggsave("/fast/AG_Forslund/rob/PROSPER/immune_paper_data/figs/paper_fig_s13_a_mdc_res_microbiome_vs_ITJs.tiff", 
       width = 6, height = 5, plot = mibi_vs_ITJs_all, dpi = 900)





## Per time point #################################
(mibi_vs_ITJs_per_tp <- BuildHeatmap_2(bind_rows(bind_rows(metad_mb_ITJ_tp1) %>% mutate(groupingVar = "TP1"),
                                                 bind_rows(metad_mb_ITJ_tp2) %>% mutate(groupingVar = "TP2"),
                                                 bind_rows(metad_mb_ITJ_tp3) %>% mutate(groupingVar = "TP3")) %>%
                                         filter(metaVariable %in% c("ITJ1", "ITJ2", "ITJ3")) %>%
                                         mutate(#status = ifelse(grepl("C:", status), yes = "OK_sd", no = status),
                                           feature = gsub("unclassified_", "Unclassified ", feature),
                                           feature = gsub("Observed", "Observed taxa", feature),
                                           feature = gsub("Shannon", "Shannon diversity", feature),
                                           feature = gsub("dmm_mixture_immune_", "DMM cluster ", feature),
                                           feature = if_else(grepl("_| ", feature) |
                                                               feature == "Other",
                                                             true = glue("{feature}"),
                                                             false = glue("<i>{feature}</i>")),
                                           feature = factor(feature, levels = unique(c("Shannon diversity", "Observed taxa", feature)))), 
                                       q_cutoff = 0.05, d_cutoff = 0.2, starSize = 7,# starNudge_y = -0.5,
                                       metaVariableNames = display_names,
                                       coloring = 1,
                                       reOrder = "none",
                                       keepMeta = c("ITJ1", "ITJ2", "ITJ3")) +
   ggtitle("ITJ associations per time point") +
   theme(axis.text.x = element_text(size = 14),
         axis.text.y = element_markdown(size = 14),
         axis.title = element_blank(),
         plot.subtitle = element_blank(),
         plot.title = element_blank(),
         axis.ticks = element_blank(),
         axis.line = element_blank(),
         strip.text.x = element_text(size = 14),
         legend.text = element_text(size = 12),
         legend.title = element_text(size = 12)))
ggsave("/fast/AG_Forslund/rob/PROSPER/immune_paper_data/figs/paper_fig_s13_b_mdc_res_microbiome_vs_ITJs_per_tp.pdf", 
       width = 9, height = 3, plot = mibi_vs_ITJs_per_tp)
ggsave("/fast/AG_Forslund/rob/PROSPER/immune_paper_data/figs/paper_fig_s13_b_mdc_res_microbiome_vs_ITJs_per_tp.tiff", 
       width = 9, height = 3, plot = mibi_vs_ITJs_per_tp, dpi = 900)


# uploading data ######################################
# move samples used in this study:
used_samples <- c(list.files("/fast/AG_Forslund/rob/PROSPER/fastq_files/prosper/decontaminated/",
                             pattern = paste0("^(", paste(ps_prosper_rf_filter_raref@sam_data$X.SampleID, collapse = "|"), ")"),
                             full.names = F),
                  list.files("/fast/AG_Forslund/rob/PROSPER/fastq_files/primal/decontaminated/",
                             pattern = paste0("^(", paste(ps_prosper_rf_filter_raref@sam_data$PRIMAL_sample_ID, collapse = "|"), ")"),
                             full.names = F),
                  list.files("/fast/AG_Forslund/rob/PROSPER/fastq_files/RawData-Ext145/RawData/decontaminated/",
                             pattern = paste0("^(", paste(ps_prosper_rf_filter_raref@sam_data$X.SampleID, collapse = "|"), ")"),
                             full.names = F))
file.copy(used_samples, "/fast/AG_Forslund/rob/PROSPER/immune_paper_data/fastq_files/")

metadata_for_publishing <- ps_prosper_rf_filter_raref %>% sample_data() %>% data.frame() %>%
  select(X.SampleID, PRIMAL_sample_ID, Pseudonym, dmm_mixture, time_point_immune, GA_days, day_of_life) %>%
  mutate(sample_name = ifelse(is.na(X.SampleID), yes = PRIMAL_sample_ID, no = X.SampleID),
         library_ID = sample_name,
         title = "16S sequencing of human feces",
         library_strategy = "AMPLICON",
         library_source = "METAGENOMIC",
         library_selection = "PCR",
         library_layout = "paired",
         platform = "ILLUMINA",
         instrument_model = "Illumina MiSeq",
         design_description = "Libraries were prepared using primers targeting the V4 region of the 16S rRNA gene",
         filetype = "fastq",
         .keep = "unused")
metadata_for_publishing$filename <- map_chr(metadata_for_publishing$sample_name, 
                                            ~list.files("/fast/AG_Forslund/rob/PROSPER/immune_paper_data/fastq_files/",
                                               pattern = paste0("^", .x, ".*dec_1.fastq.gz")))
metadata_for_publishing$filename2 <- map_chr(metadata_for_publishing$sample_name,
                                             ~list.files("/fast/AG_Forslund/rob/PROSPER/immune_paper_data/fastq_files/",
                                                         pattern = paste0("^", .x, ".*dec_2.fastq.gz")))


biosample_attributes_for_publishing <- ps_prosper_rf_filter_raref %>% sample_data() %>% data.frame() %>%
  select(X.SampleID, PRIMAL_sample_ID) %>%
  mutate(sample_name = ifelse(is.na(X.SampleID), yes = PRIMAL_sample_ID, no = X.SampleID),
         source_material_id = sample_name,
         organism = "human gut metagenome",
         host = "Homo sapiens",
         isolation_source = "feces",
         collection_date = "2019/2022",
         geo_loc_name = "Germany:Wuerzburg",
         geo_loc_name_country = "Germany",
         geo_loc_name_country_continent = "Europe",
         lat_lon = "49.7913 N 9.9534 E",
         .keep = "unused")


write_tsv(metadata_for_publishing, "/fast/AG_Forslund/rob/PROSPER/immune_paper_data/data/metadata.tsv")
write_tsv(biosample_attributes_for_publishing, "/fast/AG_Forslund/rob/PROSPER/immune_paper_data/data/biosample_attributes.tsv")
