library(TSRexploreR)
library(filesstrings)
library(tidyverse)
library(GenomicRanges)
library(readxl)
library(plyranges)
library(rtracklayer)
library(org.Dr.eg.db)
library(clusterProfiler)
library(ReactomePA)

#########################
### Data preparation ####
#########################

# Bash command to get CAGE bigWigs
# wget -r -np -nd -A plus.bw,minus.bw -I /downloads/danRer7/ http://zeprome.genereg.net/downloads/danRer7/

samples <- data.frame(sample_name = before_last_dot(before_last_dot(list.files("bigwigs", pattern = "plus"))), 
                      file_1 = list.files("bigwigs", full.names = TRUE, pattern = "plus"), 
                      file_2 = list.files("bigwigs", full.names = TRUE, pattern = "minus"),
                      condition = before_last_dot(before_last_dot(list.files("bigwigs", pattern = "plus"))))

# Read in annotation
annotation <- file.path("danRer7.refGene.gtf")

# Create TSRexploreR object
exp <- tsr_explorer(sample_sheet = samples, genome_annotation = annotation)

# Import TSS bigWigs
exp <- tss_import(exp)

# Change negative-strand TSS negative scores to positive
corrected_scores <- map(exp@experiment$TSSs, function(x) {
  as.data.frame(x) %>%
    mutate(score = ifelse(strand == "-", score * -1, score)) %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)
})

exp@experiment$TSSs <- corrected_scores

# Format TSSs
exp <- format_counts(exp, data_type = "tss")

# Thresholding
exp <- annotate_features(exp, data_type = "tss", feature_type = "transcript", 
                         upstream = 500, downstream = 500)

threshold_data <- explore_thresholds(exp, steps = 1, max_threshold = 20, 
                                     samples = c(
                                       "unfertilized.egg", "fertilized.egg", "64cells", 
                                       "512cells", "high", "oblong", "sphere.dome", "30p.dome",
                                       "shield", "somites", "prim6", "prim6.rep1",
                                       "prim6.rep2", "prim20"), 
                                     use_normalized = FALSE)

plot_threshold_exploration(threshold_data, ncol = 4, point_size = 1) +
  scale_color_viridis_c() +
  xlim(0,20)

# Cluster TSSs
exp <- tss_clustering(exp, threshold = 5, max_distance = 25, max_width = 250)

#########################
### Shifting analysis ###
#########################

### Unfertilized egg pairwise vs. all other time points

# Remove unfertilized egg from the list of samples to loop over in order to prevent self-comparison
the_shift_list <- filter(exp@meta_data$sample_sheet, sample_name != "unfertilized.egg") %>%
  dplyr::select(sample_name) %>%
  as_tibble

# Perform pairwise shift analyses
for(sample in the_shift_list$sample_name) {
  exp <- tss_shift(exp, 
                   sample_1 = c(TSS = "unfertilized.egg", TSR = "unfertilized.egg"),
                   sample_2 = c(TSS = sample, TSR = sample),
                   max_distance = 100, min_threshold = 10, n_resamples = 1000L,
                   comparison_name = str_c("unfertilized.egg.vs.", sample))
}

# Plot numbers of detected significant shifts for each comparison
map(exp@shifting$results, function(x) {
  x %>% 
  mutate(shift_status = case_when(
    shift_score < 0 ~ "upstream",
    shift_score > 0 ~ "downstream",
    TRUE ~ "unshifted"
  )) %>%
    group_by(shift_status) %>%
    count
}) %>% 
  as_tibble %>% 
  pivot_longer(everything(), names_to = "sample", values_to = "n") %>%
  mutate(sample = fct_relevel(sample, 
                              "unfertilized.egg.vs.fertilized.egg", "unfertilized.egg.vs.64cells", 
                              "unfertilized.egg.vs.512cells", "unfertilized.egg.vs.high", 
                              "unfertilized.egg.vs.oblong", "unfertilized.egg.vs.sphere.dome", 
                              "unfertilized.egg.vs.30p.dome", "unfertilized.egg.vs.shield", 
                              "unfertilized.egg.vs.somites", "unfertilized.egg.vs.prim6",
                              "unfertilized.egg.vs.prim6.rep1", "unfertilized.egg.vs.prim6.rep2",
                              "unfertilized.egg.vs.prim20")) %>%
  ggplot(aes(x = sample, fill = n$shift_status, y = n$n)) +
  geom_col(position = "stack") +
  scale_fill_viridis_d(end = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Generate a rank plots of shift scores
plot_shift_rank(exp) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  scale_fill_viridis_c()

samples <- c("unfertilized.egg.vs.fertilized.egg", "unfertilized.egg.vs.64cells", 
                  "unfertilized.egg.vs.512cells", "unfertilized.egg.vs.high", 
                  "unfertilized.egg.vs.oblong", "unfertilized.egg.vs.sphere", 
                  "unfertilized.egg.vs.30p_dome", "unfertilized.egg.vs.shield", 
                  "unfertilized.egg.vs.somites", "unfertilized.egg.vs.prim6",
                  "unfertilized.egg.vs.prim6.rep1", "unfertilized.egg.vs.prim6.rep2",
                  "unfertilized.egg.vs.prim20")

# Plot shift counts
plot_shift_count(exp, samples = samples) +
  scale_fill_viridis_d(end = 0.5)
          
# Annotate shifted TSRs
exp <- annotate_features(exp, data_type = "shift", feature_type = "transcript", 
                         upstream = 500, downstream = 500)

# Determine genomic distribution of shifts and plot
distribution <- genomic_distribution(exp, data_type = "shift")

plot_genomic_distribution(distribution) +
  scale_fill_viridis_d()

mbt_shift_sig <- filter(exp@shifting$results$pre_vs_post_mbt, shift_score != 0)

# Export merged shifts
export.bed(mbt_shift_sig, "tsrexplorer_shifts.bed")

#############################
### Functional annotation ###
#############################

# GO BP analysis
upstream_shift_genes <- exp@shifting$results$unfertilized.egg.vs.prim20 %>% 
  dplyr::filter(abs(distanceToTSS) <= 500 & shift_score < 0) %>% 
  dplyr::select(geneId)

enrichGO(upstream_shift_genes$geneId, OrgDb = org.Dr.eg.db, ont = "BP", keyType = "SYMBOL") %>%
  dotplot(showCategory = 10) +
  scale_color_viridis_c()

downstream_shift_genes <- exp@shifting$results$unfertilized.egg.vs.prim20 %>% 
  dplyr::filter(abs(distanceToTSS) <= 500 & shift_score > 0) %>% 
  dplyr::select(geneId) %>%
  as.data.frame

enrichGO(downstream_shift_genes$geneId, OrgDb = org.Dr.eg.db, ont = "BP", keyType = "SYMBOL") %>%
  dotplot(showCategory = 10) +
  scale_color_viridis_c()

# Reactome pathway analysis
upstream_shift_genes_entrez <- bitr(upstream_shift_genes$geneId, fromType = "SYMBOL",
                                    toType = "ENTREZID", OrgDb = org.Dr.eg.db) %>%
  dplyr::select(ENTREZID)

enrichPathway(upstream_shift_genes_entrez$ENTREZID, organism = "zebrafish") %>%
  dotplot(showCategory = 10) +
  scale_color_viridis_c()

downstream_shift_genes_entrez <- bitr(downstream_shift_genes$geneId, fromType = "SYMBOL",
                                      toType = "ENTREZID", OrgDb = org.Dr.eg.db) %>%
  dplyr::select(ENTREZID)

enrichPathway(downstream_shift_genes_entrez$ENTREZID, organism = "zebrafish") %>%
  dotplot(showCategory = 10) +
  scale_color_viridis_c()

#########################
### CAGEr comparisons ###
#########################

# CAGEr shifting promoters import
cager_shifts_gr <- read_xls("41586_2014_BFnature12974_MOESM327_ESM-2.xls", skip = 1) %>%
  dplyr::select(2:5) %>%
  makeGRangesFromDataFrame

export.bed(cager_shifts_gr, "~/Desktop/zf_cage/cager_shifts.bed")

overlap <- cager_shifts_gr %>%
  group_by_overlaps(makeGRangesFromDataFrame(exp@shifting$results$pre_vs_post_mbt, keep.extra.columns = TRUE))

export.bed(overlap, "~/Desktop/zf_cage/overlap_shifts.bed")

# Promoters to plot
# lipf/rnls
