#-----------------------------------------------------#
# analysis of resource-use data from biolog plates ####
#-----------------------------------------------------#

# recreates Figure 2 and Figure 3 of the manuscript
# creates Figure S1 and Figure S4

# clear workspace
rm(list = ls())

# load packages ####
library(tidyverse)
library(lme4)
library(patchwork)
library(widyr)
library(corrr)
library(MicrobioUoE) # remotes::install_github('padpadpadpad/MicrobioUoE')

# figure path
path_fig <- 'figures'

# load in data ####
d <- read.csv('data/resource_use_data.csv', stringsAsFactors = FALSE) %>%
  mutate(time = parse_number(tp))

#-----------------------------------------#
# looking at what time point to choose ####
#-----------------------------------------#

# OD590 was only taken at time points 4, 5 & 6
# OD600 was taken at time points 1-6

# plot change in OD through time for OD600
ggplot(filter(d, od_wave == 600), aes(time, od_cor, group = interaction(sample, od_wave))) +
  geom_line(alpha = 0.1) +
  facet_wrap(~substrate, ncol = 4, labeller = labeller(substrate = number_facets)) +
  theme_bw(base_size = 12) +
  labs(
       x = 'Sampling time point',
       y = expression(OD[600])) +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size = 9))

# save this plot out
ggsave(file.path(path_fig, 'Figure_S1.png'), last_plot(), height = 11, width = 9)
ggsave(file.path(path_fig, 'Figure_S1.pdf'), last_plot(), height = 11, width = 9)

# look at variation between measurements at t590 and t600
d_var <- filter(d, time > 3) %>%
  group_by(., substrate, sample, time) %>%
  summarise(mean = mean(od_cor),
            diff = max(od_cor) - min(od_cor)) %>%
  ungroup()

# difference between wavelengths is almost always a magnitude smaller than the mean. They are very similar measurements. Use OD590 (time point 4) because that is what the Biolog team recommend.

#--------------------------------------------#
# do quality control on chosen time point ####
#--------------------------------------------#

# filter to keep t4 590. This was somewhat arbritary but see Figure S1 above
d_t4_590 <- filter(d, od_wave == 590 & tp == 'T4')

# 33 and 45 were identified as outliers that did not grow

# remove outliers
d_t4_590 <- filter(d_t4_590,! sample %in% c('33', '45'))

# add column for ranked mean OD_cor
d_t4_590 <- group_by(d_t4_590, substrate) %>%
  mutate(mean_od = mean(od_cor)) %>%
  ungroup() %>%
  mutate(., rank = dense_rank(desc(mean_od)))

# filter so that only wells where growth could occur were considered
d_t4_590 <- filter(d_t4_590, mean_od > 0.05)
# gives us 18 substrates

# plot Figure 2A
plot1 <- ggplot(d_t4_590) +
  geom_line(aes(forcats::fct_reorder(substrate, rank), od_cor, group = sample, col = evolved), alpha = 1) +
  theme_bw(base_size = 16) +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 25, hjust = 1, size = 10)) +
  ylab(expression(OD[590])) +
  xlab('substrate') +
  scale_color_manual(values = c('#e41a1c', 'dark grey', 'black')) +
  ggtitle('(a) Resource-use of pre-adapted clones')

# data wrangle
d_t4_590 <- mutate(d_t4_590, rep = parse_number(gsub('-', '', population)),
                   population2 = as.factor(case_when(evolved == 'without_community' ~ paste('pre-adapted without nmc', abs(rep), sep = ' '),
                                                    evolved == 'with_community' ~ paste('pre-adapted with nmc', abs(rep), sep = ' '),
                                                    evolved == 'ancestor' ~ 'ancestor')),
                   population2 = forcats::fct_shift(population2, 1))

# create Figure S4
ggplot(d_t4_590) +
  geom_line(aes(forcats::fct_reorder(as.character(rank), rank), od_cor, group = sample, col = evolved), alpha = 0.5) +
  stat_summary(aes(forcats::fct_reorder(as.character(rank), rank), od_cor, group = population, col = evolved), fun = mean, geom = "line") +
  theme_bw(base_size = 16) +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 25, hjust = 1, size = 9),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        panel.grid = element_blank()) +
  ylab(expression(OD[590])) +
  xlab('substrate rank') +
  facet_wrap(~ population2, labeller = labeller(population2 = MicrobioUoE::letter_facets), ncol = 3) +
  scale_color_manual(values = c('#e41a1c', 'dark grey', 'black'))

ggsave(file.path(path_fig, 'Figure_S4.pdf'), last_plot(), height = 11, width = 9)
ggsave(file.path(path_fig, 'Figure_S4.png'), last_plot(), height = 11, width = 9)

#----------------------------------#
# Calculate phenotypic variance ####
#----------------------------------#

# Phenotypic diversity, VP , is calculated as the average Euclidean distance across all pairs of clones in a population (Hall & Colegrave, 2006). For two genotypes from the same population, the Euclidean distance is the square root of the sum of the squared differences in OD between the two genotypes across all substrates. Biologically, this measures the differences in the metabolic profiles of two genotypes.

# Take all the distance matrices of each pairwise combination of clones in each 
# average Euclidean distance
pop_dists_df <- filter(d_t4_590, evolved %in% c('with_community', 'without_community', 'ancestor')) %>%
  group_by(., evolved, population) %>%
  pairwise_dist(sample, substrate, od_cor, upper = FALSE) %>%
  rename(clone_i = item1, clone_j = item2)

# create average phenotypic diversity per population
V_P <- group_by(pop_dists_df, evolved, population) %>%
  summarise(., V_P = mean(distance)) %>%
  data.frame()

# calculate phenotypic variance at the treatment level (pre-adaptation), not the population level
V_P_24clones <- filter(d_t4_590, evolved %in% c('with_community', 'without_community', 'ancestor')) %>%
  group_by(., evolved) %>%
  pairwise_dist(sample, substrate, od_cor, upper = FALSE) %>%
  rename(clone_i = item1, clone_j = item2) %>%
  ungroup() %>%
  group_by(evolved) %>%
  summarise(., V_P = mean(distance)) %>%
  data.frame()

#---------------------------------------#
# calculate V_G - genotypic variance ####
#---------------------------------------#

# Genotypic (clonal in the manuscript) variance VG is calculated as the average variance of clone performance across all substrates (Venail et al., 2008).

# for each substrate, we calculated the variance in performance among the clones of each population
# this was averaged across substrates to create the within-population mean genotypic variation
V_G <- group_by(d_t4_590, evolved, substrate, population) %>%
  summarise(V_G = var(od_cor)) %>%
  data.frame()
V_G_pop <- group_by(V_G, evolved, population) %>%
  summarise(V_G = mean(V_G)) %>%
  data.frame()

#-------------------------------------------#
# calculate V_E - environmental variance ####
#-------------------------------------------#

# Environmental variance, VE is calculated as the variance in the average clone performance across all substrates.

# calculate mean performance on each substrate per population
# calculate variance of mean performance across substrates
V_E <- group_by(d_t4_590, evolved, substrate, population) %>%
  summarise(V_E = mean(od_cor)) %>%
  data.frame()
V_E_pop <- group_by(V_E, evolved, population) %>%
  summarise(V_E = var(V_E)) %>%
  data.frame()

#----------------------------------------------------#
# Calculate G x E interaction for each population ####
#----------------------------------------------------#

# see Barrett et al. 2005 Am Nat and Venail et al. 2008 Nature

# 1. calculate responsiveness - indicates differences in the environmental variances and thus measures diversity of resource exploitation strategies (specialists and generalists)
# sum (sd_j - sd_i)^2/(2*n_genotypes(n_genotypes - 1))

# create dataframe for standard deviation per clone across environments
d_sd <- group_by(d_t4_590, evolved, population, sample) %>%
  summarise(., sd_E = sd(od_cor)) %>%
  data.frame(stringsAsFactors = FALSE)

# create 2 copies of this for merging later
sd_j_clone <- dplyr::rename(d_sd, clone_j = sample, sd_j = sd_E) 
sd_i_clone <- rename(d_sd, clone_i = sample, sd_i = sd_E)

# create every pairwise combination of 1:n (clones/genotypes) for each population
d_R <- group_by(d_sd, evolved, population) %>%
  do(data.frame(expand.grid(clone_j = .$sample, clone_i = .$sample, stringsAsFactors = FALSE))) %>%
  ungroup() %>%
  filter(., clone_j > clone_i) %>%
  merge(., sd_j_clone, by = c('clone_j', 'evolved', 'population')) %>%
  merge(., sd_i_clone, by = c('clone_i', 'evolved', 'population'))

# calculate R for each pairwise combination
d_R <- group_by(d_R, evolved, population) %>%
  mutate(., R_comb = (sd_j - sd_i)^2/(2*n())*(n()-1)) %>%
  ungroup()

# calculate responsiveness for each population
# sum of all the pairwise combinations
d_R_pop <- group_by(d_R, evolved, population) %>%
  summarise(., R_pop = sum(R_comb)) %>%
  data.frame()

# 2. calculate inconsistency: Inconsistency indicates non-correlations between genotypes over different environments and indicated resource specialisation.

# calculate correlations between each pairs of clones within a population
d_pearson <- group_by(d_t4_590, evolved, population) %>%
  nest() %>%
  mutate(., cor_col = purrr::map(data, pairwise_cor, item = sample, feature = substrate, value = od_cor, upper = FALSE)) %>%
  unnest_legacy(cor_col) %>%
  rename(clone_i = item1, clone_j = item2) 

# merge dataframe to responsiveness dataframe
d_inconsist <- merge(d_pearson, sd_i_clone, by = c('evolved', 'population', 'clone_i'), all.x = TRUE) %>%
  merge(., sd_j_clone, by = c('evolved', 'population', 'clone_j'), all.x = TRUE) %>%
  group_by(., evolved, population) %>%
  mutate(., i = (sd_j*sd_i*(1-correlation))/(n()*(n()-1))) %>%
  summarise(., I_pop = sum(i),
            pear_pop = mean(correlation)) %>%
  data.frame()

#----------------------------#
# analyse sequencing data ####
#----------------------------#

# read in data
d_snp <- read.csv('data/genetic_changes.csv', stringsAsFactors = FALSE, row.names = NULL)

# calculate genomic diversity
d_snp_diversity <- group_by(d_snp, preadaptation_population, evolved_with_without_community, pos) %>%
  summarise(prop = sum(af)/n()) %>%
  ungroup() %>%
  mutate(p = prop,
         q = 1-p,
         div = 1 - p^2 - q^2) %>%
  group_by(preadaptation_population, evolved_with_without_community) %>%
  summarise(., diversity = sum(div)) %>%
  ungroup()

# calculate genomic diversity for at the treatment level
d_snp_diversity_24clones <- group_by(d_snp, evolved_with_without_community, pos) %>%
  summarise(prop = sum(af)/n()) %>%
  ungroup() %>%
  mutate(p = prop,
         q = 1-p,
         div = 1 - p^2 - q^2) %>%
  group_by(evolved_with_without_community) %>%
  summarise(., diversity = sum(div)) %>%
  ungroup()

#-----------------#
# make figures ####
#-----------------#

# plot V_G
V_G_plot <- ggplot(V_G_pop, aes(evolved, V_G)) +
  geom_pretty_boxplot(aes(evolved, V_G, col = evolved, fill = evolved), filter(V_G_pop, evolved != 'ancestor')) +
  geom_point(aes(evolved, V_G, col = evolved), shape = 21, fill = 'white', size = 5, position = position_jitter(width = 0.1), filter(V_G_pop, evolved != 'ancestor')) +
  geom_point(aes(evolved, V_G), size = 7, col = '#e41a1c', filter(V_G_pop, evolved == 'ancestor')) +
  ylab('clonal variance') +
  xlab('') +
  theme_bw(base_size = 16) +
  theme(legend.position = 'none') +
  ggtitle(expression((b)~Clonal~variance~(V[C]))) +
  scale_x_discrete(labels = c('LacZ\nancestor', 'pre-adapted\nwith nmc', 'pre-adapted\nwithout nmc')) +
  scale_color_manual('', values = c('dark grey', 'black')) +
  scale_fill_manual('', values = c('dark grey', 'black'))

# plot V_E
V_E_plot <- ggplot(V_E_pop, aes(evolved, V_E)) +
  geom_pretty_boxplot(aes(evolved, V_E, col = evolved, fill = evolved), filter(V_E_pop, evolved != 'ancestor')) +
  geom_point(aes(evolved, V_E, col = evolved), shape = 21, fill = 'white', size = 5, position = position_jitter(width = 0.1), filter(V_E_pop, evolved != 'ancestor')) +
  geom_point(aes(evolved, V_E), size = 7, col = '#e41a1c', filter(V_E_pop, evolved == 'ancestor')) +
  ylab('environmental variance') +
  xlab('') +
  theme_bw(base_size = 16) +
  theme(legend.position = 'none',
        axis.text.x = element_text(size = 12, color = 'black')) +
  ggtitle(expression((b)~Environmental~variance~(V[E]))) +
  scale_x_discrete(labels = c('LacZ\nancestor', 'pre-adapted\nwith nmc', 'pre-adapted\nwithout nmc')) +
  scale_color_manual('', values = c('dark grey', 'black')) +
  scale_fill_manual('', values = c('dark grey', 'black'))

# Plot V_P
V_P_plot <- ggplot(V_P, aes(evolved, V_P)) +
  geom_pretty_boxplot(aes(evolved, V_P, col = evolved, fill = evolved), filter(V_P, evolved != 'ancestor')) +
  geom_point(aes(evolved, V_P, col = evolved), shape = 21, fill = 'white', size = 3, position = position_jitter(width = 0.1), filter(V_P, evolved != 'ancestor')) +
  geom_point(aes(evolved, V_P, col = evolved), shape = 21, fill = 'white', size = 7, filter(V_P_24clones, evolved != 'ancestor')) +
  geom_point(aes(evolved, V_P), col = '#e41a1c', size = 7, filter(V_P, evolved == 'ancestor')) +
  ylab('phenotypic variance') +
  xlab('') +
  theme_bw(base_size = 16) +
  theme(legend.position = 'none') +
  ggtitle(expression((a)~Phenotypic~variance~(V[P]))) +
  scale_x_discrete(labels = c('LacZ\nancestor', 'pre-adapted\nwith nmc', 'pre-adapted\nwithout nmc')) +
  scale_color_manual('', values = c('dark grey', 'black')) +
  scale_fill_manual('', values = c('dark grey', 'black')) +
  ylim(c(0,0.8))

# Plot responsiveness
r_plot <- ggplot(d_R_pop, aes(evolved, R_pop)) +
  geom_pretty_boxplot(aes(evolved, R_pop, col = evolved, fill = evolved), filter(d_R_pop, evolved != 'ancestor')) +
  geom_point(aes(evolved, R_pop, col = evolved), shape = 21, fill = 'white', size = 5, position = position_jitter(width = 0.1), filter(d_R_pop, evolved != 'ancestor')) +
  geom_point(aes(evolved, R_pop), col = '#e41a1c', size = 7, filter(d_R_pop, evolved == 'ancestor')) +
  ylab('responsiveness') +
  xlab('') +
  theme_bw(base_size = 16) +
  theme(legend.position = 'none') +
  ggtitle(expression((c)~Responsiveness)) +
  scale_x_discrete(labels = c('LacZ\nancestor', 'pre-adapted\nwith nmc', 'pre-adapted\nwithout nmc')) +
  scale_color_manual('', values = c('dark grey', 'black')) +
  scale_fill_manual('', values = c('dark grey', 'black'))

# plot inconsistency
I_plot <- ggplot(d_inconsist, aes(evolved, I_pop)) +
  geom_pretty_boxplot(aes(evolved, I_pop, col = evolved, fill = evolved), filter(d_inconsist, evolved != 'ancestor')) +
  geom_point(aes(evolved, I_pop, col = evolved), shape = 21, fill = 'white', size = 5, position = position_jitter(width = 0.1), filter(d_inconsist, evolved != 'ancestor')) +
  geom_point(aes(evolved, I_pop), col = '#e41a1c', size = 7, filter(d_inconsist, evolved == 'ancestor')) +
  ylab('Inconsistency') +
  xlab('') +
  theme_bw(base_size = 16) +
  theme(legend.position = 'none') +
  ggtitle(expression((d)~Inconsistency)) +
  scale_x_discrete(labels = c('LacZ\nancestor', 'pre-adapted\nwith nmc', 'pre-adapted\nwithout nmc')) +
  scale_color_manual('', values = c('dark grey', 'black')) +
  scale_fill_manual('', values = c('dark grey', 'black'))

figure_2 <- plot1 + {V_G_plot + r_plot + I_plot} + plot_layout(nrow = 2, heights = c(0.6, 0.4))

figure_2

# save plot, other ways are available
ggsave(file.path(path_fig, 'Figure_2.png'), figure_2, height = 12, width = 14)
ggsave(file.path(path_fig, 'Figure_2.pdf'), figure_2, height = 12, width = 14)

# plot genetic diversity
d_blank <- tibble(evolved_with_without_community = 'ancestor', diversity = 0)

snp_plot <- ggplot(d_snp_diversity, aes(evolved_with_without_community, diversity)) +
  geom_pretty_boxplot(aes(fill = evolved_with_without_community, col = evolved_with_without_community)) +
  geom_point(aes(col = evolved_with_without_community), fill = 'white', size = 3, position = position_jitter(width = 0.1), shape = 21) +
  geom_point(aes(col = evolved_with_without_community),data = d_snp_diversity_24clones, fill = 'white', shape = 21, size = 7) +
  geom_point(aes(evolved_with_without_community, diversity), col = '#e41a1c', size = 7, d_blank) +
  ylab('genomic diversity') +
  ggtitle('(b) Genomic diversity') +
  xlab('') +
  theme_bw(base_size = 16) +
  theme(legend.position = 'none') +
  scale_x_discrete(labels = c('LacZ\nancestor', 'pre-adapted\nwith nmc', 'pre-adapted\nwithout nmc')) +
  scale_color_manual('', values = c('dark grey', 'black')) +
  scale_fill_manual('', values = c('dark grey', 'black'))

figure_3 <- V_P_plot + snp_plot

ggsave(file.path(path_fig, 'Figure_3.png'), figure_3, height = 5, width = 10)
ggsave(file.path(path_fig, 'Figure_3.pdf'), figure_3, height = 5, width = 10)

#-------------------------------------------#
# analyse different metrics of diversity ####
#-------------------------------------------#

# 1. phenotypic diversity
mod_vp <- lm(V_P ~ evolved, filter(V_P, evolved != 'ancestor'))
mod_vp2 <- lm(V_P ~ 1, filter(V_P, evolved != 'ancestor'))
anova(mod_vp, mod_vp2)

# do both a t test and mann whitney u test for significance
vp_t <- t.test(filter(V_P, evolved != 'ancestor') %>% pull(V_P), 
               mu = filter(V_P, evolved == 'ancestor') %>% pull(V_P))
vp_mw <- wilcox.test(filter(V_P, evolved != 'ancestor') %>% pull(V_P), 
                mu = filter(V_P, evolved == 'ancestor') %>% pull(V_P))

# look at the percentage of samples greater than the ancestor
sum(filter(V_P, evolved != 'ancestor') %>% pull(V_P) > filter(V_P, evolved == 'ancestor') %>% pull(V_P))/length(filter(V_P, evolved != 'ancestor') %>% pull(V_P)) *100

# 2. genotypic/clonal diversity
mod_vg <- lm(V_G ~ evolved, filter(V_G_pop, evolved != 'ancestor'))
mod_vg2 <- lm(V_G ~ 1, filter(V_G_pop, evolved != 'ancestor'))
anova(mod_vg, mod_vg2)
vg_t <- t.test(filter(V_G_pop, evolved != 'ancestor') %>% pull(V_G), 
               mu = filter(V_G_pop, evolved == 'ancestor') %>% pull(V_G))
vg_mw <- wilcox.test(filter(V_G_pop, evolved != 'ancestor') %>% pull(V_G), 
                     mu = filter(V_G_pop, evolved == 'ancestor') %>% pull(V_G))

# look at the percentage of samples greater than the ancestor
sum(filter(V_G_pop, evolved != 'ancestor') %>% pull(V_G) > filter(V_G_pop, evolved == 'ancestor') %>% pull(V_G))/length(filter(V_G_pop, evolved != 'ancestor') %>% pull(V_G)) *100

# 3. environmental diversity
mod_ve <- lm(V_E ~ evolved, filter(V_E_pop, evolved != 'ancestor'))
mod_ve2 <- lm(V_E ~ 1, filter(V_E_pop, evolved != 'ancestor'))
anova(mod_ve, mod_ve2)
ve_t <- t.test(filter(V_E_pop, evolved != 'ancestor') %>% pull(V_E), 
               mu = filter(V_E_pop, evolved == 'ancestor') %>% pull(V_E))
ve_mw_1 <- wilcox.test(filter(V_E_pop, evolved == 'without_community') %>% pull(V_E), 
                     mu = filter(V_E_pop, evolved == 'ancestor') %>% pull(V_E))
ve_mw_2 <- wilcox.test(filter(V_E_pop, evolved == 'with_community') %>% pull(V_E), 
                       mu = filter(V_E_pop, evolved == 'ancestor') %>% pull(V_E))

# 4. responsiveness
mod_r <- lm(R_pop ~ evolved, filter(d_R_pop, evolved != 'ancestor'))
mod_r2 <- lm(R_pop ~ 1, filter(d_R_pop, evolved != 'ancestor'))
anova(mod_r, mod_r2)
r_t <- t.test(filter(d_R_pop, evolved != 'ancestor') %>% pull(R_pop), 
               mu = filter(d_R_pop, evolved == 'ancestor') %>% pull(R_pop))
r_mw <- wilcox.test(filter(d_R_pop, evolved != 'ancestor') %>% pull(R_pop), 
                     mu = filter(d_R_pop, evolved == 'ancestor') %>% pull(R_pop))

# look at the percentage of samples greater than the ancestor
sum(filter(d_R_pop, evolved != 'ancestor') %>% pull(R_pop) > filter(d_R_pop, evolved == 'ancestor') %>% pull(R_pop))/length(filter(d_R_pop, evolved != 'ancestor') %>% pull(R_pop)) *100

# 5. inconsistency
mod_i <- lm(I_pop ~ evolved, filter(d_inconsist, evolved != 'ancestor'))
mod_i2 <- lm(I_pop ~ 1, filter(d_inconsist, evolved != 'ancestor'))
anova(mod_i, mod_i2)
i_t <- t.test(filter(d_inconsist, evolved != 'ancestor') %>% pull(I_pop), 
              mu = filter(d_inconsist, evolved == 'ancestor') %>% pull(I_pop))
i_mw <- wilcox.test(filter(d_inconsist, evolved != 'ancestor') %>% pull(I_pop), 
                    mu = filter(d_inconsist, evolved == 'ancestor') %>% pull(I_pop))

# look at the percentage of samples greater than the ancestor
sum(filter(d_inconsist, evolved != 'ancestor') %>% pull(I_pop) > filter(d_inconsist, evolved == 'ancestor') %>% pull(I_pop))/length(filter(d_inconsist, evolved != 'ancestor') %>% pull(I_pop))*100

# 6. genomic diversity
mod_snp <- lm(diversity ~ evolved_with_without_community, d_snp_diversity)
mod_snp2 <- lm(diversity ~ 1, d_snp_diversity)
anova(mod_snp, mod_snp2)

snp_withnmc_mw <- wilcox.test(filter(d_snp_diversity, evolved_with_without_community == 'With_Community') %>% pull(diversity), 
                    mu = filter(d_snp_diversity_24clones, evolved_with_without_community == 'With_Community') %>% pull(diversity))
sum(filter(d_snp_diversity_24clones, evolved_with_without_community == 'With_Community') %>% pull(diversity) > filter(d_snp_diversity, evolved_with_without_community == 'With_Community') %>% pull(diversity))/length(filter(d_snp_diversity, evolved_with_without_community == 'With_Community') %>% pull(diversity))*100

snp_withoutnmc_mw <- wilcox.test(filter(d_snp_diversity, evolved_with_without_community == 'Without_Community') %>% pull(diversity), 
            mu = filter(d_snp_diversity_24clones, evolved_with_without_community == 'Without_Community') %>% pull(diversity))

# look at the percentage of samples greater than the ancestor
sum(filter(d_snp_diversity_24clones, evolved_with_without_community == 'Without_Community') %>% pull(diversity) > filter(d_snp_diversity, evolved_with_without_community == 'Without_Community') %>% pull(diversity))/length(filter(d_snp_diversity, evolved_with_without_community == 'Without_Community') %>% pull(diversity))*100

# p adjust p values if needed
d_pvals <- tibble(var = c('V_P', 'V_G', 'V_E1', 'V_E2', 'r', 'I', 'genomic'),
                  raw_pval = c(vp_mw$p.value, vg_mw$p.value, ve_mw_1$p.value, ve_mw_2$p.value, r_mw$p.value, i_mw$p.value, snp_withoutnmc_mw$p.value)) %>%
  mutate(., fdr = p.adjust(raw_pval, 'fdr'))

d_pvals
