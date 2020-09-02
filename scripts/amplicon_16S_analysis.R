#-------------------------------------------#
# Analysis of community composition data ####
#-------------------------------------------#

# recreates Figure 2 and Figure 3 of the manuscript
# creates Figure S2 and Figure S3

# clean workspace
rm(list = ls())

# load packages ####
library(DESeq2) # BiocManager::install("DESeq2")
library(phyloseq)
library(vegan)
library(patchwork) 
library(lme4)
library(flextable)
library(tidyverse)
library(adephylo)
library(ape)

# if not installed, install mctoolsr run remotes::install_github('leffj/mctoolsr')

# source extra functions
source('scripts/extra_functions.R')

#-------------------------------------#
# setup workspace and load in data ####
#-------------------------------------#

# set seed
set.seed(42)

# figure path
path_fig <- 'figures'

# load cleaned phyloseq object - removed any phyla which is assigned NA 
ps <- readRDS('data/phyloseq_16S_clean.rds')

# replace metadata with new metadata
meta_new <- read.csv('data/metadata.csv', stringsAsFactors = FALSE)
row.names(meta_new) <- meta_new$SampleID
sample_data(ps) <- sample_data(meta_new)

# show available ranks in the dataset
rank_names(ps)

# look at the number of reads per sample
sample_sums(ps)
min(sample_sums(ps)) # min of 28,000 Woof.

# initial subsetting - do not want nmc_t0 or wt_ancestor
to_keep <- filter(meta_new, ! treatment %in% c('nmc_t0', 'wt_ancestor'))
ps2 <- prune_samples(to_keep$SampleID, ps)
sample_sums(ps2)

# P. fluorescens SBW25 16S read
SBW25 = "ACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTCAAAACTGACTGACTAGAGTATGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTGATACTGACACTGAGGTGCGAAAGCGTGGGGAGCAA"

# remove SBW25 from the dataset
ps_SBW25 <- subset_taxa(ps2, rownames(tax_table(ps2)) %in% c(SBW25))
ps2 <- subset_taxa(ps2, ! rownames(tax_table(ps2)) %in% c(SBW25))

#--------------------------------------------------------------------------------#
# Analysis of the preadaptation with and without nmc on community composition ####
#--------------------------------------------------------------------------------#

# perfectly replicated for 1, 4 and 24 clones
# do these separately and then p-value correct

# 1. individual clones

# samples to keep - just individual clones
to_keep <- filter(meta_new, treatment %in% c('individual_clone'))
ps_sub <- prune_samples(to_keep$SampleID, ps2)

# transform counts to relative abundances for ordination
# weighted Unifrac - relative abundances
ps_sub_prop <- transform_sample_counts(ps_sub, function(x){x / sum(x)})

# make a data frame of the sample data
d_samp <- data.frame(sample_data(ps_sub_prop))
d_samp <- mutate(d_samp, nclones_fac = paste('C', n_clones, sep = '_'),
                 nclones_fac = as.factor(nclones_fac),
                 evol_fac = as.factor(evolution),
                 preadapt_pop = as.factor(preadapt_pop),
                 pop2 = preadapt_pop) %>%
  column_to_rownames(., 'SampleID')
levels(d_samp$pop2) <- letters[1:length(levels(d_samp$pop2))]

# calculate distance matrix
ps_wunifrac <- phyloseq::distance(ps_sub_prop, method = 'wunifrac')

# see whether individual clones within replicates are on average closer to each other than clones in other replicates

# calculate average distance of each clone to its sympatric populations and allopatric populations across reps
d_wunifrac <- dist_2_df(ps_wunifrac) %>%
  merge(., select(d_samp, X1_evol = evolution, X1_preadapt_pop = preadapt_pop) %>% rownames_to_column(., var = 'X1'), by = 'X1') %>%
  merge(., select(d_samp, X2_evol = evolution, X2_preadapt_pop = preadapt_pop) %>% rownames_to_column(., var = 'X2'), by = 'X2') %>%
  filter(., X1_evol == X2_evol) %>%
  mutate(same_rep = ifelse(X1_preadapt_pop == X2_preadapt_pop, 'Y', 'N')) %>%
  gather(., 'random', 'clone', c(X1, X2)) %>%
  group_by(clone, same_rep, X1_evol) %>%
  summarise(., mean = mean(dist)) %>%
  ungroup()

# check n - is correct
group_by(d_wunifrac, same_rep, X1_evol) %>%
  tally()

labels <- c(with_community = '(a) pre-adapted with nmc', without_community = '(b) pre-adapted without nmc')

ggplot(d_wunifrac, aes(same_rep, mean)) +
  MicrobioUoE::geom_pretty_boxplot(fill = 'black', col = 'black') +
  geom_point(fill = 'white', shape = 21, size = 3, position = position_jitter(width = 0.1)) +
  facet_wrap(~ X1_evol, labeller = labeller(X1_evol = labels)) +
  theme_bw(base_size = 16) +
  ylab('weighted Unifrac distance') +
  xlab('') +
  scale_x_discrete(labels = c('allopatric pair', 'sympatric pair')) +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0))

ggsave(file.path(path_fig, 'Figure_S3.png'), last_plot(), height = 5, width = 8)
ggsave(file.path(path_fig, 'Figure_S3.pdf'), last_plot(), height = 5, width = 8)

# See whether pre-adaptation treatment alters effect on community composition
mod_div_1 <- vegan::adonis(ps_wunifrac ~ evol_fac, data = d_samp, n_perm = 9999)
mod_betadisper_1 <- betadisper(ps_wunifrac, d_samp$evol_fac)

# 2. 4 related clones

# filter treatments to keep
to_keep <- filter(meta_new, treatment %in% c('4_related_clones'))
ps_sub <- prune_samples(to_keep$SampleID, ps2)

# transform counts to relative abundances for ordination
ps_sub_prop <- transform_sample_counts(ps_sub, function(x){x / sum(x)})

# make a data frame of the sample data
d_samp <- data.frame(sample_data(ps_sub_prop))
d_samp <- mutate(d_samp, nclones_fac = paste('C', n_clones, sep = '_'),
                 nclones_fac = as.factor(nclones_fac),
                 evol_fac = as.factor(evolution)) %>%
  column_to_rownames(., 'SampleID')

# calculate distance matrix
ps_wunifrac <- phyloseq::distance(ps_sub_prop, method = 'wunifrac')

# run permutational anova
mod_div_4 <- vegan::adonis(ps_wunifrac ~ evol_fac, data = d_samp, n_perm = 9999)
mod_betadisper_4 <- betadisper(ps_wunifrac, d_samp$evol_fac)

# 3. 24 related clones

# filter samples to keep
to_keep <- filter(meta_new, treatment %in% c('evolved_with_community', 'evolved_without_community'))
ps_sub <- prune_samples(to_keep$SampleID, ps2)

# transform counts to relative abundances for ordination
ps_sub_prop <- transform_sample_counts(ps_sub, function(x){x / sum(x)})

# prepare metadata
d_samp <- data.frame(sample_data(ps_sub_prop))
d_samp <- mutate(d_samp, nclones_fac = paste('C', n_clones, sep = '_'),
                 nclones_fac = as.factor(nclones_fac),
                 evol_fac = as.factor(evolution)) %>%
  column_to_rownames(., 'SampleID')

# calculate distance matrix
ps_wunifrac <- phyloseq::distance(ps_sub_prop, method = 'wunifrac')

# run an Adonis test
mod_div_24 <- vegan::adonis(ps_wunifrac ~ evol_fac, data = d_samp, n_perm = 9999)
mod_betadisper_24 <- betadisper(ps_wunifrac, d_samp$evol_fac)

# Make Figure 5
# Looking at effect of pre-adaptation history across different levels of diversity #

# make one big plot of all clones
to_keep <- filter(meta_new, ! treatment %in% c('4_unrelated_clones'))
ps_sub <- prune_samples(to_keep$SampleID, ps2)

# transform counts to relative abundances for ordination
ps_sub_prop <- transform_sample_counts(ps_sub, function(x){x / sum(x)})

# wrangle the metadata
d_samp <- data.frame(sample_data(ps_sub_prop))
d_samp <- mutate(d_samp, nclones_fac = paste('C', n_clones, sep = '_'),
                 nclones_fac = as.factor(nclones_fac),
                 evol_fac = as.factor(evolution)) %>%
  column_to_rownames(., 'SampleID') %>%
  unite(., 'id', nclones_fac, evol_fac, sep = ':')

# calculate distance matrix
ps_wunifrac <- phyloseq::distance(ps_sub_prop, method = 'wunifrac')

# run a betadisper
mod_betadisper <- betadisper(ps_wunifrac, d_samp$id)

# get correct eigenvalue contributions by applying a correction
correct_eigenvalues <- ape::pcoa(ps_wunifrac, correction = 'cailliez', d_samp$id) %>%
  .$values %>%
  pull(Rel_corr_eig)
correct_eigenvalues[1:2]

# grab centroids and other data
d_fig_preadapt <- get_betadisper_data(mod_betadisper)

# combine centroid and eigenvector dataframes for plotting
betadisper_lines <- merge(select(d_fig_preadapt$centroids, group, PCoA1, PCoA2), select(d_fig_preadapt$eigenvector, group, PCoA1, PCoA2), by = c('group'))

# add distances to eigenvector and lines data
betadisper_lines <- mutate(betadisper_lines, distances = dist_between_points(PCoA1.x, PCoA2.x, PCoA1.y, PCoA2.y))
d_fig_preadapt$eigenvector$distances <- d_fig_preadapt$distances$distances

# split up group into clones and evolution context
betadisper_lines <- separate(betadisper_lines, group, c('nclones', 'evol'), sep =':') %>%
  mutate(., evol = case_when(evol == 'NA' & nclones == 'C_1' ~ 'lacz_ancestor',
                             nclones == 'C_high' ~ 'negative_control',
                             TRUE ~ evol),
         nclones =forcats::fct_relevel(nclones,
                                   "C_1", "C_4", "C_24"))
d_fig_preadapt$centroids <- separate(d_fig_preadapt$centroids, group, c('nclones', 'evol'), sep =':') %>%
  mutate(., evol = case_when(evol == 'NA' & nclones == 'C_1' ~ 'lacz_ancestor',
                             nclones == 'C_high' ~ 'negative_control',
                             TRUE ~ evol),
         nclones =forcats::fct_relevel(nclones,
                                       "C_1", "C_4", "C_24"))
d_fig_preadapt$eigenvector <- separate(d_fig_preadapt$eigenvector, group, c('nclones', 'evol'), sep =':') %>%
  mutate(., evol = case_when(evol == 'NA' & nclones == 'C_1' ~ 'lacz_ancestor',
                             nclones == 'C_high' ~ 'negative_control',
                             TRUE ~ evol),
         nclones =forcats::fct_relevel(nclones,
                                       "C_1", "C_4", "C_24"))

d_fig_preadapt$eigenvalue <- mutate(d_fig_preadapt$eigenvalue, percent = eig/sum(eig)*100)

# clone labels
facet <- c(C_1 = '(a) single clone', C_4 = '(b) 4 clones', C_24 = '(c) 24 clones')

# plot PCoA
fig_preadapt <- ggplot() +
  # add negative control and lacz ancestor into background
  geom_point(aes(PCoA1, PCoA2, alpha = 1 - distances), select(filter(d_fig_preadapt$eigenvector, evol %in% c('lacz_ancestor')), -nclones), size = 0.75, col = '#e41a1c', alpha = 0.4) +
  geom_point(aes(PCoA1, PCoA2), select(filter(d_fig_preadapt$centroids, evol %in% c('lacz_ancestor')), -nclones), size = 5, col = '#e41a1c', alpha = 0.4) +
  # add points from preadapted clone treatments
  geom_point(aes(PCoA1, PCoA2, col = evol, alpha = 1 - distances), filter(d_fig_preadapt$eigenvector, ! evol %in% c('lacz_ancestor', 'negative_control')), size = 1) +
  geom_segment(aes(x = PCoA1.x, y = PCoA2.x, yend = PCoA2.y, xend = PCoA1.y, group = row.names(filter(betadisper_lines, ! evol %in% c('lacz_ancestor', 'negative_control'))), col = evol, alpha = 1 - distances), filter(betadisper_lines, ! evol %in% c('lacz_ancestor', 'negative_control'))) +
  geom_point(aes(PCoA1, PCoA2, col = evol), filter(d_fig_preadapt$centroids, ! evol %in% c('lacz_ancestor', 'negative_control')), size = 7) +
  theme_bw(base_size = 14, base_family = 'Helvetica') +
  ylab('PCoA Axis 2 (11.9 %)') +
  xlab('PCoA Axis 1 (26.5 %)') +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0)) +
  facet_wrap(~ nclones, labeller = labeller(nclones = facet)) +
  scale_alpha(range = c(0.0001, 1), guide = FALSE) +
  scale_color_manual('', values = c('dark grey', 'black'))

# save plot, other ways are available
ggsave(file.path(path_fig, 'Figure_5.png'), fig_preadapt, height = 5, width = 12)
ggsave(file.path(path_fig, 'Figure_5.pdf'), fig_preadapt, height = 5, width = 12)

# make Table S4: effect of pre-adaptation within levels of diversity
d_table <- tibble(clone = c(1, 4, 24), mod = list(mod_div_1, mod_div_4, mod_div_24))

tidy_adonis <- function(adonis_mod){
  stuff <- adonis_mod$aov.tab %>% data.frame(stringsAsFactors = FALSE) %>%
    tibble::rownames_to_column(var = 'factor') %>%
    janitor::clean_names() %>%
    rename(p_value = pr_f)
  return(stuff)
}

d_table <- mutate(d_table, output = purrr::map(mod, tidy_adonis)) %>%
  unnest(output) %>%
  filter(., factor != 'Total') %>%
  spread(factor, df) %>%
  fill(., Residuals, .direction = 'up') %>%
  filter(!is.na(p_value)) %>%
  unite(., 'd.f.', c(evol_fac, Residuals), sep = ', ') %>%
  select(clone, f_model, `d.f.`, r2, p_value) %>%
  mutate_at(., vars(f_model, r2), function(x) round(x, 2)) %>%
  mutate_all(., as.character)

# super_fp
super_fp <- officer::fp_text(vertical.align = "superscript", font.size = 16, font.family = 'Times')

# italics_fp
italic_fp <- officer::fp_text(italic = TRUE, font.size = 16, font.family = 'Times')

table <- flextable(d_table) %>%
  align(align = 'center', part = 'all') %>%
  set_header_labels(clone = 'Clonal diversity',
                    f_model = "F statistic",
                    p_value = "p value") %>%
  compose(., j = "r2", part = "header", 
          value = as_paragraph("R", as_chunk("2", props = super_fp))) %>%
  compose(., j = "d.f.", part = "header", 
          value = as_paragraph(as_chunk("d.f.", props = italic_fp))) %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 16, part = 'all') %>%
  autofit()

# save as png
save_as_image(table, 'Table_S4.png', webshot = 'webshot2', zoom = 3)

#---------------------------------------#
# Analysis of allopatry vs. sympatry ####
#---------------------------------------#

# analyse four clone mixes to see if there are differences between coevolved populations and allopatric populations

# filter treatments to keep
to_keep <- filter(meta_new, treatment %in% c('4_related_clones', '4_unrelated_clones'))
ps_sub <- prune_samples(to_keep$SampleID, ps2)

# transform counts to relative abundances for ordination
ps_sub_prop <- transform_sample_counts(ps_sub, function(x){x / sum(x)})

# make a data frame of the sample data
d_samp <- data.frame(sample_data(ps_sub_prop))
d_samp <- mutate(d_samp, nclones_fac = paste('C', n_clones, sep = '_'),
                 nclones_fac = as.factor(nclones_fac),
                 evol_fac = as.factor(evolution),
                 treatment_fac = as.factor(treatment)) %>%
  column_to_rownames(., 'SampleID')

# calculate distance matrix
ps_wunifrac <- phyloseq::distance(ps_sub_prop, method = 'wunifrac')

# run an Adonis test
mod_allopatry <- vegan::adonis(ps_wunifrac ~ treatment_fac, data = d_samp, n_perm = 9999)
mod_betadisper <- betadisper(ps_wunifrac, d_samp$treatment_fac)

plot(mod_betadisper)

# not significant but need to make the plot

#--------------------------------------------#
# Looking for overall effect of diversity ####
#--------------------------------------------#

# change some values of nclones
meta_new <- sample_data(ps2) %>% data.frame() %>%
  mutate(., n_clones = paste('C_', readr::parse_number(n_clones), sep = '')) %>%
  mutate(., n_clones = case_when(treatment == 'lacz_ancestor' ~ 'lacz_ancest',
                                 treatment == 'negative_control' ~ 'C_high',
                                 TRUE ~ as.character(n_clones)))
row.names(meta_new) <- meta_new$SampleID
sample_data(ps2) <- sample_data(meta_new)

# Look for impact of diversity within pre-adaptation history
# pre-adapted with nmc 1 vs 4 vs 24 clones
# pre-adapted without nmc 1 vs 4 vs 24 clones

# first preadapted with the microbial community

# filter samples
to_keep <- filter(meta_new, ! treatment %in% c('nmc_t0',  'wt_ancestor', 'negative_control')) %>%
  filter(evolution == 'with_community')
ps_sub <- prune_samples(to_keep$SampleID, ps2)

# make counts proportions
ps_prop <- transform_sample_counts(ps_sub, function(x){x / sum(x)})

# define distance metric
metric = 'wunifrac'

# get the distance matrix out of the data
ps_dist <- phyloseq::distance(ps_prop, method = metric)

# make a data frame of the sample data
d_samp <- data.frame(sample_data(ps_prop))
d_samp <- mutate(d_samp,
                 nclones_fac = as.factor(n_clones)) %>%
  column_to_rownames(., 'SampleID')

# run an Adonis test
mod_nclonefac <- vegan::adonis(ps_dist ~ nclones_fac, data = d_samp, n_perm = 9999)

# run a multiple comparison to see which treatments are different
mult_comp1 <- mctoolsr::calc_pairwise_permanovas(ps_dist, d_samp, 'nclones_fac', n_perm = 9999)

mult_comp1

# second preadapted without the microbial community

# filter samples
to_keep <- filter(meta_new, ! treatment %in% c('nmc_t0',  'wt_ancestor', 'negative_control')) %>%
  filter(evolution == 'without_community')
ps_sub <- prune_samples(to_keep$SampleID, ps2)

# make counts proportions
ps_prop <- transform_sample_counts(ps_sub, function(x){x / sum(x)})

# define distance metric
metric = 'wunifrac'

# get the distance matrix out of the data
ps_dist <- phyloseq::distance(ps_prop, method = metric)

# make a data frame of the sample data
d_samp <- data.frame(sample_data(ps_prop))
d_samp <- mutate(d_samp,
                 nclones_fac = as.factor(n_clones)) %>%
  column_to_rownames(., 'SampleID')

# run an Adonis test
mod_nclonefac <- vegan::adonis(ps_dist ~ nclones_fac, data = d_samp, n_perm = 9999)

# run a multiple comparison to see which treatments are different
mult_comp2 <- mctoolsr::calc_pairwise_permanovas(ps_dist, d_samp, 'nclones_fac', n_perm = 9999)

mult_comp2

# wrangle the multiple comparison datasets
mult_comp2 <- mutate(mult_comp2, evolution = 'without_community')
mult_comp1 <- mutate(mult_comp1, evolution = 'with_community')

mult_comp <- bind_rows(mult_comp1, mult_comp2)

d_num <- group_by(d_samp, nclones_fac, evolution) %>% tally() 
d_num <- ungroup(d_num) %>%
  bind_rows(., d_num) %>%
  mutate(evolution = rep(c('with_community', 'without_community'), each = 3))
d_num[3,'n'] <- 5

# create Table S3
mult_comp_to_save <- merge(mult_comp, rename(d_num, X1 = nclones_fac, n_1 = n), by = c('X1', 'evolution')) %>%
  merge(., rename(d_num, X2 = nclones_fac, n_2 = n), by = c('X2', 'evolution')) %>%
  unite(., 'n', c(n_1, n_2), sep = ', ') %>%
  mutate(., X1 = forcats::fct_recode(X1, `single clone` = 'C_1',
                                     `4 clones` = 'C_4',
                                     `24 clones` = 'C_24'),
         X2 = forcats::fct_recode(X2, `single clone` = 'C_1',
                                  `4 clones` = 'C_4',
                                  `24 clones` = 'C_24'),
         contrast = paste(X1, 'vs.', X2, sep = ' '),
         evolution = ifelse(evolution == 'with_community', 'pre-adapted with nmc', 'pre-adapted without nmc')) %>%
  select(., evolution, contrast, n, everything(),-c(X1, X2, pvalBon)) %>%
  mutate_at(., vars(4:ncol(.)), function(x) signif(x, 2)) %>%
  arrange(evolution)

# make table
super_fp <- officer::fp_text(vertical.align = "superscript", font.size = 16, font.family = 'Times')
sub_fp <- officer::fp_text(vertical.align = "subscript", font.size = 16, font.family = 'Times')
italic_fp <- officer::fp_text(italic = TRUE, font.size = 16, font.family = 'Times')
border_fp <- officer::fp_border(color="black", width = 1)

table <- flextable(select(mult_comp_to_save, evolution, contrast, n, R2, pval, pvalFDR)) %>%
  align(j = c('contrast'), align = 'left', part = 'all') %>%
  align(j = c('R2', 'pval', 'pvalFDR', 'n'), align = 'center', part = 'all') %>%
  add_footer_row(., colwidths = 6, values = '') %>% 
  bold(i = ~ pval < 0.05, j = ~ pval) %>%
  bold(i = ~ pvalFDR < 0.05, j = ~ pvalFDR) %>%
  set_header_labels(evolution = 'pre-adaptation treatment',
                    n = 'number of samples per treatment',
                    pval = "raw p value") %>%
  flextable::compose(., j = "R2", part = "header", 
                     value = as_paragraph("R", as_chunk("2", props = super_fp))) %>%
  flextable::compose(., j = "pvalFDR", part = "header", 
                     value = as_paragraph("p", as_chunk("adj", props = sub_fp))) %>%
  flextable::compose(., j = "contrast", part = "footer", 
                     value = as_paragraph("p", as_chunk("adj", props = sub_fp), 'was calculated using the false discovery rate method')) %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 16, part = 'all') %>%
  autofit() %>%
  padding(padding.top = 5, part = 'footer') %>%
  hline(i = 3, border = border_fp)

# get a png from the html file with webshot
save_as_image(table, 'Table_S3.png', zoom = 3, webshot = 'webshot2')

# now all levels of diversity pooled

# filter samples
to_keep <- filter(meta_new, ! treatment %in% c('nmc_t0',  'wt_ancestor', 'negative_control'))
ps_sub <- prune_samples(to_keep$SampleID, ps2)

# make counts proportions
ps_prop <- transform_sample_counts(ps_sub, function(x){x / sum(x)})

# define distance metric
metric = 'wunifrac'

# get the distance matrix out of the data
ps_dist <- phyloseq::distance(ps_prop, method = metric)

# make a data frame of the sample data
d_samp <- data.frame(sample_data(ps_prop))
d_samp <- mutate(d_samp,
                 nclones_fac = as.factor(n_clones)) %>%
  column_to_rownames(., 'SampleID')

# run an Adonis test
mod_nclonefac <- vegan::adonis(ps_dist ~ nclones_fac, data = d_samp, n_perm = 9999)

# run a multiple comparison to see which treatments are different
mult_comp <- mctoolsr::calc_pairwise_permanovas(ps_dist, d_samp, 'nclones_fac', n_perm = 9999)

d_num <- group_by(d_samp, nclones_fac) %>% tally()

# create multiple comparison table
mult_comp_to_save <- merge(mult_comp, rename(d_num, X1 = nclones_fac, n_1 = n), by = 'X1') %>%
  merge(., rename(d_num, X2 = nclones_fac, n_2 = n), by = 'X2') %>%
  unite(., 'n', c(n_1, n_2), sep = ', ') %>%
  mutate(., X1 = forcats::fct_recode(X1, `single clone` = 'C_1',
                                     `4 clones` = 'C_4',
                                     `24 clones` = 'C_24',
                                     `negative control` = 'C_high'),
         X2 = forcats::fct_recode(X2, `LacZ ancestor` = 'lacz_ancest',
                                  `4 clones` = 'C_4',
                                  `24 clones` = 'C_24',
                                  `negative control` = 'C_high'),
         contrast = paste(X1, 'vs.', X2, sep = ' ')) %>%
  select(., contrast, n, everything(),-c(X1, X2)) %>%
  mutate_at(., vars(3:ncol(.)), function(x) signif(x, 2)) %>%
  arrange(desc(contrast))

# make table
super_fp <- officer::fp_text(vertical.align = "superscript", font.size = 16, font.family = 'Times')
sub_fp <- officer::fp_text(vertical.align = "subscript", font.size = 16, font.family = 'Times')
italic_fp <- officer::fp_text(italic = TRUE, font.size = 16, font.family = 'Times')

table <- flextable(select(mult_comp_to_save, contrast, n, R2, pval, pvalFDR)) %>%
  align(j = c('contrast'), align = 'left', part = 'all') %>%
  align(j = c('R2', 'pval', 'pvalFDR', 'n'), align = 'center', part = 'all') %>%
  add_footer_row(., colwidths = 5, values = '') %>% 
  bold(i = ~ pval < 0.05, j = ~ pval) %>%
  bold(i = ~ pvalFDR < 0.05, j = ~ pvalFDR) %>%
  set_header_labels(n = 'number of samples per treatment',
                    pval = "raw p value") %>%
  flextable::compose(., j = "R2", part = "header", 
                     value = as_paragraph("R", as_chunk("2", props = super_fp))) %>%
  flextable::compose(., j = "pvalFDR", part = "header", 
                     value = as_paragraph("p", as_chunk("adj", props = sub_fp))) %>%
  flextable::compose(., j = "contrast", part = "footer", 
                     value = as_paragraph("p", as_chunk("adj", props = sub_fp), 'was calculated using the false discovery rate method')) %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 16, part = 'all') %>%
  autofit() %>%
  padding(padding.top = 5, part = 'footer')

# get a png from the html file with webshot
save_as_image(table, 'Table_S2.png', zoom = 3, webshot = 'webshot2')

# overwrite metadata to allow plotting of nclones_fac
sample_data(ps_prop) <- sample_data(d_samp)

# beta-diversity analysis - look at homogeneity of variances
mod1_dispers <- betadisper(ps_dist, d_samp$nclones_fac)

# get correct eigenvalue contributions by applying a correction
correct_eigenvalues <- ape::pcoa(ps_dist, correction = 'cailliez', d_samp$nclones_fac) %>%
  .$values %>%
  pull(Rel_corr_eig)
correct_eigenvalues[1:2]

# plot of model
plot(mod1_dispers)
boxplot(mod1_dispers)

# anova
anova(mod1_dispers)

# Permutation test for F
pmod <- permutest(mod1_dispers, pairwise = TRUE)

# Tukey's Honest Significant Differences
T_HSD <- TukeyHSD(mod1_dispers)

# create Figure 4

# get betadisper data
betadisper_dat <- get_betadisper_data(mod1_dispers)

# combine centroid and eigenvector dataframes for plotting
betadisper_lines <- merge(select(betadisper_dat$centroids, group, PCoA1, PCoA2), select(betadisper_dat$eigenvector, group, PCoA1, PCoA2), by = c('group'))

# add distances to eigenvector and lines data
betadisper_lines <- mutate(betadisper_lines, distances = dist_between_points(PCoA1.x, PCoA2.x, PCoA1.y, PCoA2.y))
betadisper_dat$eigenvector$distances <- betadisper_dat$distances$distances

# plot PCoA
p1 <- ggplot() +
  geom_point(aes(PCoA1, PCoA2, col = forcats::fct_relevel(group, c('lacz_ancest', 'C_1', 'C_4', 'C_24', 'C_high')), alpha = 0.5 - distances), betadisper_dat$eigenvector, size = 1.5) +
  geom_segment(aes(x = PCoA1.x, y = PCoA2.x, yend = PCoA2.y, xend = PCoA1.y, group = row.names(betadisper_lines), col = forcats::fct_relevel(group, c('lacz_ancest', 'C_1', 'C_4', 'C_24', 'C_high')), alpha = 0.5 - distances), betadisper_lines) +
  geom_point(aes(PCoA1, PCoA2, col = group), betadisper_dat$centroids, size = 7) +
  theme_bw(base_size = 16, base_family = 'Helvetica') +
  ylab('PCoA Axis 2 (11.2%)') +
  xlab('PCoA Axis 1 (23.7%)') +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 6) +
  ggtitle('(a)')

# plot axis 1
p2 <- ggplot(betadisper_dat$eigenvector, aes(forcats::fct_relevel(group, c('lacz_ancest', 'C_1', 'C_4', 'C_24', 'C_high')), PCoA1, col = forcats::fct_relevel(group, c('lacz_ancest', 'C_1', 'C_4', 'C_24', 'C_high')), fill = forcats::fct_relevel(group, c('lacz_ancest', 'C_1', 'C_4', 'C_24', 'C_high')))) +
  geom_hline(aes(yintercept = 0)) +
  MicrobioUoE::geom_pretty_boxplot() +
  geom_point(shape = 21, fill = 'white', size = 3, position = position_jitter(width = 0.1))  +
  theme_bw(base_size = 16, base_family = 'Helvetica') +
  ylab('PCoA Axis 1 (23.7%)') +
  xlab('') +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 6) +
  scale_fill_brewer(type = 'qual', palette = 6) +
  scale_x_discrete(labels = c('LacZ\nancestor', 'single\nclone', '4 clones', '24 clones', 'negative\ncontrol')) +
  ggtitle('(b)')

# calculate proportion of samples above 1
group_by(betadisper_dat$eigenvector, group) %>%
  summarise(above_1 = sum(PCoA1 > 0)/n())

p3 <- p1 + p2 + plot_layout(ncol = 2, widths = c(1, 1))

ggsave(file.path(path_fig, 'Figure_4.png'), p3, height = 6, width = 12)
ggsave(file.path(path_fig, 'Figure_4.pdf'), p3, height = 6, width = 12)

#-------------------------------------------#
# Look at differential abundance of OTUs ####
#-------------------------------------------#

# filter samples
to_keep <- filter(meta_new, ! treatment %in% c('nmc_t0',  'wt_ancestor', 'negative_control'))
ps_SBW25 <- prune_samples(to_keep$SampleID, ps)

# replace SBW25's taxonomic rank
tax_table1 <- data.frame(tax_table(ps_SBW25), stringsAsFactors = FALSE)
tax_table1[row.names(tax_table1) == SBW25, 1:7] <- 'focal species'
tax_table(ps_SBW25) <- as.matrix(tax_table1)

sort_taxa <- names(sort(taxa_sums(ps_SBW25), TRUE)) 
match(SBW25, sort_taxa)
# 39th most abundant taxa

# filter 100 most abundant taxa
n_otu <- 100

asv_abundant <- names(sort(taxa_sums(ps_SBW25), TRUE)[1:n_otu]) 
ps_sub2 <- prune_taxa(asv_abundant, ps_SBW25)

# assign otu number by total abundance
x <- tibble(otu = names(sort(taxa_sums(ps_sub2), TRUE)),
            abundance = paste('ASV', 1:n_otu, sep = ' ')) %>%
  merge(., tibble(otu = taxa_names(ps_sub2),
                  order = 1:n_otu), by = 'otu') %>%
  arrange(order)

taxa_names(ps_sub2) <- x$abundance

# add a column for difference between C_24 and everything else
sample_data(ps_sub2)$diversity <- ifelse(sample_data(ps_sub2)$n_clones == '24', 'C_high', 'C_low')

diagdds = phyloseq_to_deseq2(ps_sub2, ~diversity)

# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
dds_effectoftreatment <- DESeq(diagdds, test="LRT", reduced= ~diversity)
diagdds = DESeq(diagdds, fitType="local")
dds_effectoftreatment <- DESeq(diagdds, test="LRT", reduced=~ 1)

res = results(diagdds, pAdjustMethod = 'fdr')
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05

sigtab = cbind(as(res, "data.frame"), as(tax_table(ps_sub2)[rownames(res), ], "matrix")) %>%
  tibble::rownames_to_column(var = 'otu') %>%
  janitor::clean_names() %>%
  mutate(., just = ifelse(log2fold_change < 0, 2, -1),
         pos = ifelse(log2fold_change < 0, log2fold_change - lfc_se, log2fold_change + lfc_se),
         phylum2 = Hmisc::capitalize(gsub('_', ' ', phylum)))

head(sigtab)

sigtab2 <- filter(sigtab, padj < alpha)

# colours
cols <- c('#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6')

# add new place for labels

# plot
ggplot(sigtab2, aes(forcats::fct_reorder(otu, log2fold_change, .desc = TRUE), log2fold_change, col = phylum2)) +
  geom_point(size = 3) +
  geom_linerange(aes(x = otu, ymin = log2fold_change - lfc_se, ymax = log2fold_change+ lfc_se)) +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  geom_text(aes(y = pos, label = round(base_mean), vjust = just), col = 'black', size = MicrobioUoE::pts(12)) +
  theme_bw(base_size = 14) +
  labs(x = 'Phylum of each ASV',
       y = 'log2 fold change in abundance') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = arrange(sigtab, desc(log2fold_change)) %>% pull(phylum2)) +
  ylim(c(-2.5, 7)) +
  theme(legend.position = c(0.85, 0.75),
        legend.background = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.5, 2), 'cm')) +
  scale_color_manual('Phylum', values = cols)
  
ggsave(file.path(path_fig, 'Figure_6.pdf'), last_plot(), height = 8, width = 10)
ggsave(file.path(path_fig, 'Figure_6.png'), last_plot(), height = 8, width = 10)

# look at phylogenetic distance of the 100 most abundant taxa from SBW25

ps_tree <- prune_taxa(c(asv_abundant, SBW25), ps)

# calculate distance between all the nodes
d_otu_phylodist <- adephylo::distTips(phy_tree(ps_tree))
d_otu_phylodist <- dist_2_df(d_otu_phylodist)

# keep only distances from SBW25
d_otu_phylodist <- filter(d_otu_phylodist, X1 == SBW25 | X2 == SBW25) %>%
  mutate(., otu = ifelse(X1 == SBW25, X2, X1)) %>%
  select(otu, dist)

# have a look at changes in abundance using DESeq2
sigtab <- dplyr::rename(sigtab2, abundance = otu) %>%
  merge(., select(x, otu, abundance), by = 'abundance') %>%
  merge(., d_otu_phylodist, by = 'otu')

ggplot(sigtab, aes(dist, log2fold_change)) +
  geom_point() +
  theme_bw()

# look at correlation
cor.test(sigtab$log2fold_change, sigtab$dist, method = 'pearson')

