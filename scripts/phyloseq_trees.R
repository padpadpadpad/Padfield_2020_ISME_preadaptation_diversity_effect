# examining the phylogenetic trees of the phyloseq object

# this script plots the phylogenetic tree of different levels of filtering of the raw phyloseq object
# code is mainly adapted from Callahan et al. (2016) F1000 Research

# load packages
library(patchwork)
library(tidyverse)
library(phyloseq)

# figure path
path_fig <- 'figures'

# load data - latest run which we are happy with ####
ps <- readRDS('data/phyloseq_16S_raw.rds')

# look at object
ps
# SO MANY TREE NODES

# first...
rank_names(ps)

# look number of features for each phyla
table(tax_table(ps)[, "Phylum"], exclude = NULL)
# lots of NAs here. These are probably artifacts and can be removed

# raw tree
psraw_tree <- plot_tree(ps, method = "treeonly",
                        ladderize = "left",
                        title = "Raw tree")

# remove NA characterisation
ps0 <- subset_taxa(ps, !is.na(Phylum))
table(tax_table(ps0)[, "Phylum"], exclude = NULL)

# plot this tree
ps0_tree <- plot_tree(ps0, method = "treeonly",
          ladderize = "left",
          title = "Tree with NA phyla removed")

# plot side by side
plot_both <- grid.arrange(psraw_tree, ps0_tree, ncol = 2)

# explore prevalence of the dataframe
# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps0),
               MARGIN = ifelse(taxa_are_rows(ps0), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(prevalence = prevdf,
                    total_abundance = taxa_sums(ps0),
                    tax_table(ps0), stringsAsFactors = FALSE) %>%
  mutate(., feature = row.names(.))

# head(prevdf)

# Are there phyla that are comprised of mostly low-prevalence features? Compute the total and average prevalences of the features in each phylum.
prev_sum <- group_by(prevdf, Phylum) %>%
  summarise(., mean_prevalence = mean(prevalence)/nsamples(ps0)*100,
            total_prevalence = sum(prevalence)) %>%
  ungroup()

prev_sum

# get unique taxa for prevalence and plot
prevdf1 <- filter(prevdf, Phylum %in% get_taxa_unique(ps0, "Phylum"))

# plot prevalence of samples by 
ggplot(prevdf1, aes(total_abundance, prevalence / nsamples(ps0),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~ Phylum) + theme(legend.position="none")

# delete things under 5% prevalence
# define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(ps0)
prevalenceThreshold

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa <- filter(prevdf1, prevalence >=  prevalenceThreshold) %>%
  pull(feature)
ps_prevfilt <- prune_taxa(keepTaxa, ps0)

# plot tree again
ps_prevfilt_tree <- plot_tree(ps_prevfilt, method = "treeonly",
                              ladderize = "left",
                              title = "Tree after prevalence filtering")


# agglomerate taxa by Genus ####

# How many genera would be present after filtering?
length(get_taxa_unique(ps_prevfilt, taxonomic.rank = "Genus"))

# agglomerate taxa
ps_genus = tax_glom(ps_prevfilt, "Genus", NArm = TRUE)

# plot tree again
ps_genus_tree <- plot_tree(ps_genus, method = "treeonly",
                              ladderize = "left",
                              title = "Tree after pooling samples by genus")

ps_genus_tree

# plot side by side
plot_all <- psraw_tree + ps0_tree + ps_prevfilt_tree + ps_genus_tree + plot_layout(ncol = 2)

ggsave(file.path(path_fig, 'phyloseq_trees.pdf'), plot_all, height = 14, width = 12)
