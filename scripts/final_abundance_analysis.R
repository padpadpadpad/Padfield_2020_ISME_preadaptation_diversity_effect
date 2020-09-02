#------------------------------------------------------------------#
# look at changes in P. fluorescens abundance across treatments ####
#------------------------------------------------------------------#

# recreates Figure 7 of the manuscript

# load packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(flextable)
library(officer)
library(emmeans)

# figure path
path_fig <- 'figures'

# load data
d <- read.csv('data/metadata.csv', stringsAsFactors = FALSE) %>%
  janitor::clean_names()

# filter out treatments that are not used in the analyses
d_sub <- filter(d, !treatment %in% c('wt_ancestor', 'nmc_t0', 'negative_control')) %>%
  mutate(., log_density = log10(density_cfu_g)) %>%
  mutate_all(., function(x)ifelse(x == 0|is.infinite(x)|is.nan(x), NA, x))

# 1. Look at differences in abundance between pre-adaptation treatments WITHIN levels of diversity ####

# pairwise of pre_adaptation treatments
d_preadapt <- filter(d_sub, !treatment %in% c('4_unrelated_clones', 'lacz_ancestor')) %>%
  dplyr::select(log_density, evolution, n_clones) %>%
  filter(!is.na(log_density))

# quick plot of data
ggplot(d_preadapt, aes(evolution, log_density)) +
  geom_point(size = 3) +
  facet_wrap(~ n_clones) +
  theme_bw()

# run linear model across levels of diversity
mod_preadapt <- d_preadapt %>%
  nest(data = c(log_density, evolution)) %>%
  mutate(., model = purrr::map(data, ~lm(log_density ~ evolution, .x)),
         tidy_model = purrr::map(model, broom::tidy)) %>%
  dplyr::select(-c(data, model)) %>%
  unnest(tidy_model) %>%
  filter(term == 'evolutionwithout_community') %>%
  mutate(., p.value = round(p.value, 3))

# look at results
mod_preadapt
# all p values are > 0.05

# as there is no effect, can pool pre-adaptation treatments together to look at overall effect of diversity

#--------------------------------------------------------------#
# 2. look at effect of diversity and preadaptation as whole ####
#--------------------------------------------------------------#

# pool preadaptation treatment together to look at overall effect of diversity
d_nclones <- dplyr::select(d_sub, log_density, treatment, n_clones) %>%
  mutate(., n_clones = ifelse(treatment == 'lacz_ancestor', 'lacz_ancestor', n_clones)) %>%
  filter(!is.na(log_density))

# quick plot 
ggplot(d_nclones, aes(n_clones, log_density)) +
  MicrobioUoE::geom_pretty_boxplot(fill = 'black', col = 'black') +
  geom_point(fill = 'white', shape = 21, size = 3) +
  theme_bw()

# problem with doing unbalanced anova is usually homogeneity of variances assumption

# look at variance across samples and correlate that against the number of points
# would expect more points equals higher variance
d_var <- group_by(d_nclones, n_clones) %>%
  summarise(n = n(),
            var = var(log_density),
            sd = sd(log_density)) %>%
  ungroup()

cor.test(d_var$n, d_var$var, method = 'spearman')
# no correlation

# even though lacZ has lower n, has higher variance than individual clone where sample size is much bigger
# this is the justification for doing unbalanced anovas and linear models

d_nclones <- mutate(d_nclones, preadapt = ifelse(n_clones == 'lacz_ancestor', 'none', 'yes'))

mod_treat <- lm(log_density ~ n_clones, d_nclones)
mod_treat2 <- lm(log_density ~ 1, d_nclones)
anova(mod_treat, mod_treat2)
# significant difference 0.013

#--------------------#
# create Figure 7 ####
#--------------------#

ggplot(d_nclones, aes(forcats::fct_relevel(n_clones, c('lacz_ancestor', '1', '4', '24')), 10^log_density, col = forcats::fct_relevel(n_clones, c('lacz_ancestor', '1', '4', '24')), fill = forcats::fct_relevel(n_clones, c('lacz_ancestor', '1', '4', '24')))) +
  MicrobioUoE::geom_pretty_boxplot() +
  geom_point(shape = 21, fill = 'white', size = 3, position = position_jitter(width = 0.1))  +
  theme_bw(base_size = 16, base_family = 'Helvetica') +
  ylab(expression(italic(P.~fluorescens)~density~(CFU~g^-1~soil))) +
  xlab('') +
  theme(legend.position = 'none') +
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)), limits = c(10^3.5, 10^7.5)) +
  scale_color_brewer(type = 'qual', palette = 6) +
  scale_fill_brewer(type = 'qual', palette = 6) +
  scale_x_discrete(labels = c('LacZ\nancestor', 'single\nclone', '4 clones', '24 clones', 'negative\ncontrol'))

ggsave(file.path(path_fig, 'Figure_7.png'), last_plot(), height = 5, width = 6)
ggsave(file.path(path_fig, 'Figure_7.pdf'), last_plot(), height = 5, width = 6)

#------------------------------#
# create table of contrasts ####
#------------------------------#

# rename contrasts
levels <- c(`single clone vs. LacZ ancestor` = '1 - lacz_ancestor',
            `24 clones vs. LacZ ancestor` = '24 - lacz_ancestor',
            `4 clones vs. LacZ ancestor` = '4 - lacz_ancestor',
            `24 clones vs. 4 clones` = '24 - 4',
            `single clone vs. 24 clones` = '1 - 24',
            `single clone vs. 4 clones` = '1 - 4')

# create contrast table
contrasts <- emmeans::emmeans(mod_treat, pairwise ~ n_clones) %>%
  .$contrasts %>%
  data.frame() %>%
  mutate(contrast = forcats::fct_recode(contrast, !!!levels)) %>%
  mutate_if(is.numeric, function(x)round(x, 2)) %>%
  mutate_all(as.character)

# super_fp
super_fp <- fp_text(vertical.align = "superscript", font.size = 8, font.family = 'Times')

# italics_fp
italic_fp <- fp_text(italic = TRUE, font.size = 16, font.family = 'Times')

# create flextable
table <- flextable(contrasts) %>%
  align(align = 'center', part = 'all') %>%
  align(align = 'left', j = 'contrast', part = 'all') %>%
  set_header_labels(t.ratio = "t ratio",
                    p.value = "p value") %>%
  add_footer_row(., colwidths = 6, values = '') %>%
  compose(., j = "df", part = "header", 
          value = as_paragraph(as_chunk("d.f.", props = italic_fp))) %>%
  compose(., j = "contrast", part = "footer", 
          value = as_paragraph('P values were adjusted using the Tukey method for comparing a family of 4 estimates')) %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 16, part = 'all') %>%
  autofit() %>%
  bold(i = ~ p.value < 0.05, j = ~ p.value) %>%
  padding(padding.top = 5, part = 'footer')

# save out table as png
save_as_image(table, path = "Table_s5.png", webshot = 'webshot2')

