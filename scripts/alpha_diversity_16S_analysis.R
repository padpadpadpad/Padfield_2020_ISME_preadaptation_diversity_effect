# look at alpha diversity and pielous evenness

# load packages ####
library(phyloseq)
library(tidyverse)
library(patchwork)

# figure path
path_fig <- 'figures'

# load data - latest run which we are happy with ####
# these files need to be there
ps <- readRDS('data/phyloseq_16S_clean.rds')
# replace metadata with new metadata
# when wanting to add columns to metadata, it is better to edit metadata_creation and overwrite the metadata file as then it can be overwritten in all future files
meta_new <- read.csv('data/metadata.csv', stringsAsFactors = FALSE)
row.names(meta_new) <- meta_new$SampleID
sample_data(ps) <- sample_data(meta_new)

# show available ranks in the dataset
rank_names(ps)

# look at the number of reads per sample
sample_sums(ps)
min(sample_sums(ps)) # min of 28,000. Woof.

# alpha diversity estimates ####

# prune OTUs that are not present in any of the samples
ps_sub <- prune_taxa(taxa_sums(ps) > 0, ps)

# metadata ###
m <- sample_data(ps_sub) %>%
  data.frame() %>%
  select(., SampleID, sample_name, treatment, evolution, preadapt_pop, n_clones) %>%
  janitor::clean_names()

# calculate diversity measures of each sample ####
a_div <- estimate_richness(ps, measures = c('Shannon', 'Observed')) %>%
  mutate(., SampleID = row.names(.)) %>%
  mutate(., pielou = Shannon / log(Observed)) %>%
  janitor::clean_names() %>%
  merge(., m, by = 'sample_id') %>%
  mutate_at(., c('sample_id', 'treatment', 'sample_name', 'evolution'), as.character) 
a_div <- select(a_div, sample_id, treatment, evolution, preadapt_pop, observed, pielou, n_clones) %>%
  filter(! treatment %in% c('wt_ancestor', 'nmc_t0', 'negative_control')) %>%
  mutate(., n_clones = paste('C_', n_clones, sep = ''))

# look at individual clones only, does evolution with and without community impact diversity and evenness ####

# run a test on evenness and diversity across evolution lines
head(a_div)

a_div_preadapt <- filter(a_div, ! treatment %in% c('lacz_ancestor', '4_unrelated_clones')) %>%
  mutate(n_clones =forcats::fct_relevel(n_clones,
                                "C_1", "C_4", "C_24"))

facet <- c(C_1 = '(a) single clone', C_4 = '(b) 4 clones', C_24 = '(c) 24 clones')

# plot alpha diversity across pre-adaptation treatments
plot_div1 <- ggplot(a_div_preadapt, aes(evolution, observed, col = evolution, fill = evolution)) +
  MicrobioUoE::geom_pretty_boxplot() +
  geom_point(shape = 21, fill = 'white', position = position_jitter(width = 0.15), size = 3) +
  scale_x_discrete(labels = c('pre-adapted\nwith nmc', 'pre-adapted\nwithout nmc')) +
  facet_wrap(~n_clones, labeller = labeller(n_clones = facet)) +
  theme_bw(base_size = 14, base_family = 'Helvetica') +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        axis.title.x = element_blank()) +
  ylim(c(0, 1400)) +
  ylab('Alpha diversity\n(observed ASVs)') +
  xlab('') +
  scale_color_manual('', values = c('dark grey', 'black')) +
  scale_fill_manual('', values = c('dark grey', 'black'))

# model these differences
d_mods <- group_by(a_div_preadapt, n_clones) %>%
  nest() %>%
  mutate(mod_alpha = purrr::map(data, ~lm(observed ~ evolution, .x)),
         mod_evenness = purrr::map(data, ~lm(pielou ~ evolution, .x)))

d_mod_alpha <- d_mods %>%
  mutate(., tidy_mod = map(mod_alpha, broom::tidy)) %>%
  unnest(tidy_mod) %>%
  filter(term == 'evolutionwithout_community') %>%
  mutate(., model = 'diversity',
         p_adj = p.adjust(p.value, method = 'fdr'))
  
#--------------------------------------#
# analysis over levels of diversity ####
#--------------------------------------#

a_div <- mutate(a_div, n_clones = ifelse(treatment == 'lacz_ancestor', 'lacz_ancestor', n_clones))

# plot
plot_div2 <- ggplot(a_div, aes(forcats::fct_relevel(n_clones, c('lacz_ancestor', 'C_1', 'C_4', 'C_24')), observed, col = forcats::fct_relevel(n_clones, c('lacz_ancestor', 'C_1', 'C_4', 'C_24')), fill = forcats::fct_relevel(n_clones, c('lacz_ancestor', 'C_1', 'C_4', 'C_24')))) +
  MicrobioUoE::geom_pretty_boxplot() +
  geom_point(shape = 21, fill = 'white', position = position_jitter(width = 0.15), size = 3) +
  theme_bw(base_size = 14, base_family = 'Helvetica') +
  theme(legend.position = 'none',
        plot.title = element_text(size = 14)) +
  scale_x_discrete(labels = c('LacZ\nancestor', 'single\nclone', '4 clones', '24 clones')) +
  scale_color_brewer(type = 'qual', palette = 6) +
  scale_fill_brewer(type = 'qual', palette = 6) +
  ylim(c(0, 1400)) +
  ylab('Alpha diversity\n(observed ASVs)') +
  xlab(c('')) +
  ggtitle('(d)')

plot_div <- plot_div1 + plot_div2 + plot_layout(ncol = 1, heights = c(0.4, 0.6))

ggsave(file.path(path_fig, 'Figure_S5.png'), plot_div, height = 8, width = 9)
ggsave(file.path(path_fig, 'Figure_S5.pdf'), plot_div, height = 8, width = 9)

mod_div <- lm(observed ~ n_clones, a_div)
mod_div2 <- lm(observed ~ 1, a_div)
anova(mod_div, mod_div2)
emmeans::emmeans(mod_div, pairwise ~ n_clones)

