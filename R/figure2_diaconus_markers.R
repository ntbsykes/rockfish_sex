# Setup ------------------------------------------------------------------------

#Install required packages
if(!require(tidyverse)){
  install.packages("tidyverse")
  library(tidyverse)
}
if(!require(cowplot)){
  install.packages("cowplot")
  library(cowplot)
}
if(!require(PNWColors)){
  install.packages("PNWColors")
  library(PNWColors)
}

# Set working directory
setwd("C:/Users/ntbsy/Desktop/rockfish/")

# Load colour palette
pal <- pnw_palette("Bay", 8) %>%
  append("grey90") %>%
  append("white")

# Loading the data -------------------------------------------------------------

# Load in markers and the sex popmap
S_diaconus_sex_markers <- read_tsv("S_diaconus/S_diaconus_markers_tidy.tsv")

S_diaconus_popmap <- read_tsv("S_diaconus/radsex/S_diaconus.popmap.tsv",
                 col_names = c("sample","sex"))

chr_map <- read_tsv("rockfish_chromosomes.tsv", 
                    col_names = c("chr", "chr_n"))

# Find mean allele frequency to see if the reference or alternate allele is the male allele.
mean_freq <- S_diaconus_sex_markers %>%  
  inner_join(S_diaconus_popmap) %>%
  group_by(chr, pos, sex) %>%
  summarize(mean_freq = mean(genotype)) 

male_allele_id <- mean_freq %>%
  ungroup() %>%
  pivot_wider(., names_from = sex, values_from = mean_freq) %>%
  mutate(male_allele = case_when(F > M ~ "ref",
                                 TRUE ~ "alt")) %>%
  select(pos,chr, male_allele)

# Plotting the heatmap (Figure 2A) ---------------------------------------------

# Chromosome 2 heatmap
pdf("rplots/diaconus_chr02.pdf", width = 8.25, height = 5.875)
plt <- S_diaconus_sex_markers %>%
  inner_join(S_diaconus_popmap) %>%
  inner_join(male_allele_id) %>%
  inner_join(chr_map) %>%
  filter(chr_n == "Chr02") %>%
  mutate(chr_n = recode(chr_n, Chr02 = "Chromosome 2")) %>%
  mutate(genotype = case_when(male_allele == "ref" ~ abs(2 - genotype),
                              TRUE ~ genotype)) %>%
  ggplot(.,aes(x = as.factor(pos), 
               y = fct_reorder(sample,sex), 
               fill = as.factor(genotype))) + 
  geom_tile() +
  facet_wrap(chr_n~sex, scales = "free") +
  scale_fill_manual(values = pal[c(1,5,8)]) +
  labs(fill = "Alleles") +
  theme_cowplot() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 4))
print(plt)
dev.off()

# Chromosome 12 heatmap
pdf("rplots/diaconus_chr12.pdf", width = 8.25, height = 5.875)
plt <- S_diaconus_sex_markers %>%
  inner_join(S_diaconus_popmap) %>%
  inner_join(male_allele_id) %>%
  inner_join(chr_map) %>%
  filter(chr_n == "Chr12") %>%
  mutate(chr_n = recode(chr_n, Chr12 = "Chromosome 12")) %>%
  mutate(genotype = case_when(male_allele == "ref" ~ abs(2 - genotype),
                              TRUE ~ genotype)) %>%
  ggplot(.,aes(x = as.factor(pos), 
               y = fct_reorder(sample,sex), 
               fill = as.factor(genotype))) + 
  geom_tile() +
  facet_wrap(chr_n~sex, scales = "free") +
  scale_fill_manual(values = pal[c(1,5,8)]) +
  labs(fill = "Alleles") +
  theme_cowplot() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 4))
print(plt)
dev.off()

# Plotting the bar graph (Figure 2B) -------------------------------------------

# Mean number of male-biased alleles per sample
pdf("rplots/diaconus_sexchr_alleles.pdf", width = 8.25, height = 5.875)
plt <- S_diaconus_sex_markers %>%
  inner_join(S_diaconus_popmap) %>%
  inner_join(male_allele_id) %>%
  inner_join(chr_map) %>%
  mutate(chr = recode(chr, 
                      Sebastes_aleutianus.PGA_scaffold_67__5_RagTag = "Chr12", 
                      Sebastes_aleutianus.PGA_scaffold_70__14_RagTag = "Chr02")) %>%
  mutate(genotype = case_when(male_allele == "ref" ~ abs(2 - genotype),
                              TRUE ~ genotype)) %>%
  group_by(sample, sex,chr) %>%
  summarize(mean_male_alleles = mean(genotype)) %>%
  ggplot(.,aes(x=fct_reorder(sample,mean_male_alleles),
               y=mean_male_alleles,
               fill=sex)) +
  geom_bar(stat="identity", width = 0.8) +
  facet_wrap(~chr, ncol = 1) +
  scale_fill_manual(values = pal[c(1,5)]) +
  theme_cowplot() +
  labs(y = "Mean number of male-biased alleles", fill = "Sex") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 5, angle = 90),
        axis.text.y = element_text(size = 6)) +
  facet_wrap(~chr, ncol = 2)
print(plt)
dev.off()
