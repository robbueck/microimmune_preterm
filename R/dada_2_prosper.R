rm(list = ls())
maindir <- "/fast/AG_Forslund/rob/PROSPER/"
setwd(maindir)
library(dada2)
library("reshape2")
library("plyr")
library("dplyr")
library("ggplot2")
library("phangorn")
library("picante")
library("btools")
library("microbiome")
library("umap")
library("gridExtra")
library("Rtsne")
library("RColorBrewer")
library("optparse")
library(cowplot)


# load additional custom functions
source("/fast/AG_Forslund/rob/mm_index/R_scripts/dada_2_functions.R")

# switch parts of the pipeline on/off
filter_step <- F
error_and_asv_step <- F
remove_chimeras_step <- F
id_taxa_step <- F
tree_step <- F
phyloseq_step <- T
rarecurve_step <- F
diversity_step <- F
ordination_step <- F

n_cores <- 5

option_list = list(
  make_option(c("-t", "--threads"), type="numeric", default=NULL))
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$threads)){
  n_cores <- T
  n_cores <- 5
} else {
  n_cores <- opt$threads
}
# filter, learn errors and denoise seperately
read_counts_raw <- read.table("read_counts_dec_prosper.txt")
library_samples <- read_counts_raw$V1
print("Filter")
if(filter_step) {
  out <- filtering_fastq(sample_names = library_samples,
                         trimLeft = 10,
                         truncLen=c(220,190),
                         maxN=0,
                         maxEE=c(2,2),
                         truncQ=2,
                         verbose = T,
                         cores = n_cores)
  saveRDS(out, paste0(maindir, "dada2/prosper/out_table.rds")) # save for each run separately
} else {
  out <- readRDS(paste0(maindir, "dada2/prosper/out_table.rds"))
}
print("error estimation and asv step")
if(error_and_asv_step){
  total_bases <- c(sum(out[,2]) * 230, sum(out[,2]) * 190)
  # because error estimation is done for forward and reverse
  st.all <- error_and_asv(sample_names = library_samples, nbases = 0.25*total_bases)
  saveRDS(st.all, paste0(maindir, "dada2/prosper/libeqtab.rds")) # save for each run separately
} else {
  st.all <- readRDS(paste0(maindir, "dada2/prosper/libeqtab.rds"))
}

# Inspect distribution of sequence lengths
cdf_data_seqlen <- melt(nchar(getSequences(st.all)),
                        value.name = "seq_length")
ggplot(cdf_data_seqlen, aes(seq_length)) +
  stat_ecdf() +
  ylab("samples")
ggsave("dada2/prosper/cfd_seq_lengths.pdf", device = "pdf")

# remove chimeras
print("remove chimeras")
if(remove_chimeras_step){
  st.all.nochim <- removeBimeraDenovo(st.all, method="consensus", multithread=n_cores, verbose=TRUE)
  saveRDS(st.all.nochim, "dada2/prosper/st_all_nochim.rds")
} else {
  st.all.nochim <- readRDS("dada2/prosper/st_all_nochim.rds")
}
dim(st.all.nochim)
sum(st.all.nochim)/sum(st.all)


# where do most reads get lost?
# load read counts for decontamination
read_counts_raw <- read.table("read_counts_raw_prosper.txt", header = F,
                              col.names = c("sample_ID", "read_counts_raw"), row.names = 1)
rownames(out) <- gsub("_dec_1.fastq.gz", "", rownames(out))
loss_table <- merge(out, read_counts_raw, by = 0)
loss_table <- merge(loss_table,
                    cbind(denoised = rowSums(st.all), removed_chimeras = rowSums(st.all.nochim)),
                    by.x = "Row.names", by.y = 0)

track <- cbind(loss_ratio(loss_table$reads.in, loss_table$read_counts_raw),  # decontamination
               loss_ratio(loss_table$reads.out, loss_table$reads.in),        # filtering
               loss_ratio(loss_table$denoised, loss_table$reads.out),        # denoising
               loss_ratio(loss_table$removed_chimeras, loss_table$denoised), # remove chimeras
               loss_ratio(loss_table$removed_chimeras, loss_table$read_counts_raw))        # total
colnames(track) <- c("decontamination", "filtering", "denoising", "remove_chimeras", "total_loss")
rownames(track) <- loss_table$Row.names
head(track)

cdf_data_step_loss <- melt(track,
                 value.name = "lost_reads")
colnames(cdf_data_step_loss)[1:2] <- c("sample", "step")
# violin plot:
ggplot(cdf_data_step_loss, aes(x=step, y = lost_reads, color = step)) + 
  geom_violin() +
  ylim(0,1) +
  geom_boxplot(width=0.1) +
  coord_flip() +
  theme(legend.position="none")
ggsave("dada2/prosper/violin_filter_ratios_per_step.pdf", device = "pdf")

# cdf_data_step_loss <- ddply(cdf_data_step_loss, .(step), transform, ecd=ecdf(lost_reads)(lost_reads))
ggplot(cdf_data_step_loss, aes(x=lost_reads, col = step)) + 
  stat_ecdf() +
  ylab("samples")
cdf_data_step_loss[is.na(cdf_data_step_loss$lost_reads),]
ggsave("dada2/prosper/cfd_filter_ratios_per_step.pdf", device = "pdf")

# save loss table to be used for further stats
loss_table <- dplyr::transmute(loss_table, Row.names = Row.names, read_counts_raw = read_counts_raw,
                               decontaminated = reads.in, filtered = reads.out,
                               denoised = denoised, removed_chimeras = removed_chimeras)
write.csv(loss_table, "dada2/prosper/read_loss_per_step.csv")

# compare human with chicken decontamination:
read_counts_raw <- read.table("read_counts_raw_prosper.txt", header = F,
                              col.names = c("sample_ID", "read_counts_raw"), row.names = 1)
read_counts_decon_h <- read.table("read_counts_dec_prosper.txt", header = F,
                              col.names = c("sample_ID", "read_counts_dec"), row.names = 1)
read_counts_decon_c <- read.table("read_counts_decon_chicken.txt", header = F,
                              col.names = c("sample_ID", "read_counts_decon_chicken"), row.names = 1)
decon_table <- cbind(read_counts_raw, read_counts_decon_h, read_counts_decon_c)
decon_table <- decon_table %>% mutate(percent_chicken = (read_counts_raw - read_counts_decon_chicken) / read_counts_raw,
                                      percent_human = (read_counts_raw - read_counts_dec) / read_counts_raw)
decon_table_long <- melt(decon_table[,c("percent_chicken", "percent_human")], value.name = "percent", )
ggplot(decon_table_long, aes(x=percent, color = variable)) +
  geom_density()
ggsave("dada2/prosper/decontamination_step.pdf", device = "pdf")

# samples with too little reads:
metadata <- read.table("/fast/AG_Forslund/rob/PROSPER/metadata/final_metadata_combined.csv", sep = ";", header = T)
sample_list <- metadata[,c("X.SampleID", "SampleNr", "Description", 
                           "MHH_sampling_ID")]
write.csv(sample_list, "prosper_samples_robert.csv")


# assign taxonomy
library(DECIPHER)
print("id taxa")
dna <- DNAStringSet(getSequences(st.all.nochim)) # Create a DNAStringSet from the ASVs
start.time_silva <- Sys.time()
if(id_taxa_step) {
  load("/fast/AG_Forslund/rob/references/silva/SILVA_SSU_r138_2019.RData") 
  ids <- IdTaxa(dna, trainingSet, processors=n_cores, verbose=FALSE)
  saveRDS(ids, "dada2/prosper/taxa_ids.rds")
} else {
  ids <- readRDS("dada2/prosper/taxa_ids.rds")
}
end.time_silva <- Sys.time()
time_silva <- end.time_silva - start.time_silva

taxa <- taxa_to_matrix(ids)
# fill internal NAs
taxa[is.na(taxa[,5]),5] <- paste("unclassified", taxa[is.na(taxa[,5]),4], sep = "_")
taxa[is.na(taxa[,4]),4] <- paste("unclassified", taxa[is.na(taxa[,4]),3], sep = "_")

print("tree step")
if (tree_step) {
  tree <- build_tree(st.all.nochim)
  saveRDS(tree, "dada2/prosper/tree.rds")
} else {
  tree <- readRDS("dada2/prosper/tree.rds")
}

# handoff to phyloseq
library(phyloseq)
library(Biostrings)
library(ggplot2)
theme_set(theme_bw())
# load metadata
metadata <- readRDS("/fast/AG_Forslund/rob/PROSPER/metadata/final_metadata_combined.rds") %>%
  filter(!is.na(X.SampleID)) %>%
  as.data.frame()



if(phyloseq_step) {
  # add metadata
  rownames(metadata) <- metadata$X.SampleID
  all(rownames(metadata) %in% rownames(st.all.nochim))
  all(rownames(st.all.nochim) %in% rownames(metadata))
  ps <- phyloseq(otu_table(st.all.nochim, taxa_are_rows=FALSE), 
                 sample_data(metadata), 
                 tax_table(taxa), tree)
  dna <- Biostrings::DNAStringSet(taxa_names(ps))
  names(dna) <- taxa_names(ps)
  ps <- merge_phyloseq(ps, dna)
  taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
  # save object
  saveRDS(ps, file = "/fast/AG_Forslund/rob/PROSPER/dada2/prosper/phyloseq_prosper.rds")
} else {
  ps <- readRDS("/fast/AG_Forslund/rob/PROSPER/dada2/prosper/phyloseq_prosper.rds")
}

# remove controls and everything below 100 reads
sample_data(ps)$sample_sum <- sample_sums(ps)
ps <- ps %>% subset_samples(sample_sum > 100) %>%
  subset_samples(!is.na(day_of_life))



ps_filter <- ps %>% filter_taxa(function(x){sum(x > 0) > 12}, TRUE)  %>%# prevalence cutoff 1%
  filter_taxa(function(x) mean(x[x > 0]) > 10, TRUE)  # mean abundance cutoff, only counting abundances > 0

ps_genus <- ps_filter %>% aggregate_taxa(level = "genus")


# distribution of abundances
prev_ab_df <- data.frame(abundace = microbiome::transform(ps_filter, transform = "compositional") %>% 
                           otu_table %>%
                           apply(., 2, function(x) {mean(x[x>0])}), 
                         prevalence = otu_table(ps_filter) %>% 
                           apply(.,2, function(x){mean(x > 0)})) %>%
  melt(., value_name = value)

decon_table_long <- melt(decon_table[,c("percent_chicken", "percent_human")], value.name = "percent", )
ggplot(prev_ab_df, aes(x = value, color = variable)) +
  geom_density() +
  ylim(0, 100)
  
# rarefaction filtering
otu_matrix <- t(otu_table(ps_genus))
class(otu_matrix) <- "matrix"
# raredat <- rarecurve(otu_matrix, step = 10)

start_time <- Sys.time()
if(rarecurve_step) {
  raredat <- rarecurve(otu_matrix, step = 50, tidy = T)
  saveRDS(raredat, "/fast/AG_Forslund/rob/PROSPER/dada2/prosper/PROSPER_rarefaction.rds")
} else {
  raredat <- readRDS("/fast/AG_Forslund/rob/PROSPER/dada2/prosper/PROSPER_rarefaction.rds")
}

end_time <- Sys.time()
print("Total time for rarecurve step:")
print(end_time - start_time)
sample_data(ps) <- sample_data(ps) %>%
  data.frame %>%
  mutate(cutoff = 2000) %>%
  sample_data()
gc()
raredat <- left_join(raredat, sample_data(ps), by = c("Site" = "X.SampleID"))
p1 <- ggplot(raredat) +
  theme_bw() +
  # scale_fill_continuous() +
  # geom_density(data = n_reads, aes(x = n_reads)) +
  scale_color_discrete(guide = "none") +
  geom_line(aes(x = Sample, y = Species, group = Site, color = matching_ID)) +
  # geom_line(aes(x = Sample, y = dens)) +
  geom_vline(data = plyr::ddply(raredat, c("age_cat"),
                                summarize, min_sample_size = min(sample_sum)),
             aes(xintercept=min_sample_size)) +
  geom_vline(aes(xintercept = cutoff), linetype = "dashed", color = "red") +
  facet_grid(~age_cat) +
  xlim(0, 20000) +
  # ylim(0,500) +
  labs(title="Rarefaction curves") + xlab("Sequenced Reads") + ylab('ASVs Detected')


# apply rarefaction filter:
ps_rf_filter <- ps %>% subset_samples(sample_sums(ps) >= cutoff)
saveRDS(ps_rf_filter, file = "/fast/AG_Forslund/rob/PROSPER/dada2/prosper/phyloseq_prosper_rf_filter.rds")