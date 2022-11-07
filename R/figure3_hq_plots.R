# Setup ------------------------------------------------------------------------

# Install required packages
if(!require(tidyverse)){
  install.packages("tidyverse")
  library(tidyverse)
}
if(!require(zoo)){
  install.packages("zoo")
  library(zoo)
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
setwd("")

# Load colour palette
pal <- pnw_palette("Bay", 8) %>%
  append("grey90") %>%
  append("white")


# Loading the data -------------------------------------------------------------

# Set global p cutoff
pvalue_cutoff = 0.005

# Make a species list vector
Sebastes_list <- c("crocotulus", "miniatus")
chr_list <- c("Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", 
              "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", 
              "Chr15", "Chr16", "Chr17", "Chr18", "Chr19", "Chr20", "Chr21", 
              "Chr22", "Chr23", "Chr24")

# Load in chromosome map
chr_map_hq <- read_tsv("rockfish_chromosomes_hq.tsv", 
                       col_names = c("Contig", "chr_n"))

# Make a tibble for assembly of all the species information
Sebastes_cm <- tibble()

#This loop loads in all the species data one at a time and builds our dataframe
for (spp in Sebastes_list){
  #Load in allele frequency data and test significance
  spp_allelefreq <- read_table(paste0("S_",spp,"/variants/S_",spp,"_allelefreq_hq.tsv"))
  spp_allelefreq <- mutate(spp_allelefreq, 
                           p = pchisq(spp_allelefreq$chi_square, 1, lower.tail = F))
  spp_allelefreq <- mutate(spp_allelefreq, 
                           signif = case_when((p) < pvalue_cutoff ~ TRUE,
                                                              TRUE ~ FALSE))
  spp_allelefreq <- mutate(spp_allelefreq, stat = "allelefreq")
  
  #Load in the heterozygosity data and test significance
  spp_hetdif <- read_table(paste0("S_",spp,"/variants/S_",spp,"_hetdif_hq.tsv"))
  spp_hetdif <- mutate(spp_hetdif, p = pchisq(spp_hetdif$chi_square, 1, lower.tail = F))
  spp_hetdif <- mutate(spp_hetdif, signif = case_when((p) < pvalue_cutoff ~ TRUE,
                                                      TRUE ~ FALSE))
  spp_hetdif <- mutate(spp_hetdif, stat = "heterozygosity")
  
  #Load in missing sex data and test significance
  spp_missingsex <- read_table(paste0("S_",spp,"/variants/S_",spp,"_missingdata_hq.tsv"))
  spp_missingsex <- mutate(spp_missingsex, 
                           p = pchisq(spp_missingsex$chi_square, 1, lower.tail = F))
  spp_missingsex <- mutate(spp_missingsex, 
                           signif = case_when((p) < pvalue_cutoff & bias > 0 ~ TRUE,
                                                              TRUE ~ FALSE))
  spp_missingsex <- mutate(spp_missingsex, stat = "missingdata")
  
  #Now bind them all together into a single data frame & label with species name
  sp <- bind_rows(spp_allelefreq, spp_hetdif, spp_missingsex) %>%
    subset(., select = c(1,2,8,9,10))
  sp <- mutate(sp, species = paste0("S_",spp,""))
  Sebastes_cm <- rbind(Sebastes_cm, sp)
}

# Plotting the windowed genome (Figure 3A) -------------------------------------

# Load in chromosome lengths dataframe, species phylogeny
chr_lengths <- read_tsv("S_miniatus/S-miniatus_SEB-74.fasta.fai",
                        col_names = c("Contig", "length", "bits", 
                                      "spacer1", "spacer2")) %>%
  arrange(Contig)

spp_phylo <- read_tsv("spp_phylogeny.tsv", col_names = TRUE)

all_data <- Sebastes_cm %>%
  filter(p != "NA") %>%
  inner_join(spp_phylo)

# Select window size
window_size <- 250

# Make a new dataframe that includes cumulative position across the genome
chr_lengths %>%
  arrange(Contig) %>%
  select(Contig, length) %>%
  mutate(total = cumsum(length)-length) %>%
  select(-length) %>%
  left_join(all_data, ., by = c("Contig" = "Contig")) %>%
  arrange(Contig, Position) %>%
  mutate(cumulative_pos = Position + total) -> all_data_cumulative

chr_lengths <- inner_join(all_data_cumulative, chr_lengths) %>%
  subset(select = c(1,10,11,12,13)) %>%
  unique()

# Make a dataframe for the background colours
chr_breaks <- tibble(start = unique(all_data_cumulative$total),
                     end = unique(all_data_cumulative$total)+(chr_lengths$length)-1,
                     colours = rep(c(0,1),245))

# Make scaled dataframes to represent S. miniatus and crocotulus data
target_spp <- c("S_crocotulus", "S_miniatus")
target_stat <- c("missingdata","allelefreq", "heterozygosity")
filt_data <- all_data_cumulative %>%
  filter(species %in% target_spp, stat %in% target_stat) %>%
  subset(select = c(3,4,6,7,8,9))

# Plot figure 3A: whole genome rolling mean significance
pdf("rplots/rockfish_genome_hq.pdf", width = 8.25, height = 2.9175)
plt <- filt_data %>%
  group_by(species) %>%
  mutate(signif = as.numeric(signif)) %>%
  mutate(average_sig = zoo::rollmean(signif, k = window_size, fill=NA)) %>%
  ggplot() +
  geom_rect(data = chr_breaks, 
            aes(xmin = start,
              xmax = end,
              ymin = - Inf,
              ymax = Inf,
              fill = colours), 
            alpha = 0.1) +
  geom_line(aes(x = cumulative_pos, 
                y = average_sig, 
                color = as.factor(top)), 
            alpha = 1, 
            size = 0.2) +
  scale_color_manual(values = pnw_palette("Bay", 8)) +
  scale_x_continuous() +
  theme_cowplot() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        text = element_text(size = 7),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size = 0.2),
        strip.text = element_text(size = 6)) +
  labs(y = "Proportion of significant sites per window", color = "Species") +
  facet_wrap(~top, ncol = 1)
print(plt)
dev.off()


# Plotting marker data for chromosome 24 (Figure 3B) ---------------------------

# Load in markers, the sex popmap, and chromosome dataframe for S. miniatus
S_miniatus_chr24_markers <- read_tsv("S_miniatus/variants/S_miniatus_chr24_markers.tidy.tsv")
S_miniatus_popmap <- read_tsv("S_miniatus/radsex/S_miniatus.popmap.tsv",
                              col_names = c("sample", "sex"))
chr_map <- read_tsv("rockfish_chromosomes_hq.tsv", 
                    col_names = c("chr", "chr_n"))

# Find mean allele frequency to see if the reference or alternate allele is the male allele.
S_miniatus_mean_freq <- S_miniatus_chr24_markers %>% 
  inner_join(S_miniatus_popmap)

# Whip up a dataframe for the number of occurences of each sample
S_miniatus_min_count <- data.frame(table(S_miniatus_mean_freq$sample)) %>%
  rename("sample" = "Var1")

# Squeeze together a mean frequency dataframe
S_miniatus_mean_freq <- inner_join(S_miniatus_min_count, S_miniatus_mean_freq) %>%
  group_by(Freq, pos, sample, sex) %>%
  summarize(mean_freq = mean(alt_count)) 

# Now we're figuring out the alternate and reference alleles
S_miniatus_male_allele_id <- S_miniatus_mean_freq %>%
  ungroup() %>%
  pivot_wider(., names_from = sex, values_from = mean_freq) %>%
  mutate(male_allele = case_when(F > M ~ "ref",
                                 TRUE ~ "alt"))

# Get all the data together in one place here
S_miniatus_data <- S_miniatus_chr24_markers %>%
  inner_join(S_miniatus_popmap) %>%
  inner_join(S_miniatus_male_allele_id) %>%
  inner_join(chr_map) %>%
  filter(chr_n == "Chr24") %>%
  mutate(chr_n = recode(chr_n, Chr24 = "Chromosome 24")) %>%
  mutate(alt_count = case_when(male_allele == "ref" ~ abs(2 - alt_count),
                               TRUE ~ alt_count)) %>%
  subset(select = c(2,3,4,5,6))

# Expand the positions list to include all of them for each sample
S_miniatus_positions_expanded <- S_miniatus_data %>%
  expand(pos, sample)

# Merge the expanded dataframe with the rest of the data, keeping all positions
S_miniatus <- merge(S_miniatus_positions_expanded, S_miniatus_data, all = T) %>%
  subset(select = c(1:3)) %>%
  inner_join(S_miniatus_popmap) %>%
  inner_join(S_miniatus_min_count)

# Replace NA values in the alt_count column with "3"
S_miniatus[is.na(S_miniatus)] <- 3

# Plot chromosome 24 heatmap - S. miniatus
my_height = 0.7

pdf("rplots/miniatus_chr24.pdf", width = 8.25, height = 5.875)
plt <-  S_miniatus %>%
  ggplot(.,aes(x = as.factor(pos), y = reorder(-Freq, sex), fill = as.factor(alt_count))) + 
  geom_tile(height = my_height) +
  facet_wrap(~sex, scales = "free", ncol = 1) +
  scale_fill_manual(values = pal[c(4,4,4,9)]) +
  theme_cowplot() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())
print(plt)
dev.off()

# Load in markers, the sex popmap, and chromosome dataframe for S. crocotulus
S_crocotulus_chr24_markers <- read_tsv("S_crocotulus/variants/S_crocotulus_chr24_markers.tidy.tsv")
S_crocotulus_popmap <- read_tsv("S_crocotulus/radsex/S_crocotulus.popmap.tsv",
                                col_names = c("sample", "sex"))

# Find mean allele frequency to see if the reference or alternate allele is the male allele.
S_crocotulus_mean_freq <- S_crocotulus_chr24_markers %>% 
  inner_join(S_crocotulus_popmap)

# Whip up a dataframe for the number of occurences of each sample
S_crocotulus_min_count <- data.frame(table(S_crocotulus_mean_freq$sample)) %>%
  rename("sample" = "Var1")

# Squeeze together a mean frequency dataframe
S_crocotulus_mean_freq <- inner_join(S_crocotulus_min_count, S_crocotulus_mean_freq) %>%
  group_by(Freq, pos, sample, sex) %>%
  summarize(mean_freq = mean(alt_count)) 

# Now we're figuring out the alternate and reference alleles
S_crocotulus_male_allele_id <- S_crocotulus_mean_freq %>%
  ungroup() %>%
  pivot_wider(., names_from = sex, values_from = mean_freq) %>%
  mutate(male_allele = case_when(F > M ~ "ref",
                                 TRUE ~ "alt"))

# Get all the data together in one place here
S_crocotulus_data <- S_crocotulus_chr24_markers %>%
  inner_join(S_crocotulus_popmap) %>%
  inner_join(S_crocotulus_male_allele_id) %>%
  inner_join(chr_map) %>%
  filter(chr_n == "Chr24") %>%
  mutate(chr_n = recode(chr_n, Chr24 = "Chromosome 24")) %>%
  mutate(alt_count = case_when(male_allele == "ref" ~ abs(2 - alt_count),
                               TRUE ~ alt_count)) %>%
  subset(select = c(2,3,4,5,6))

# Expand the positions list to include all of them for each sample
S_crocotulus_positions_expanded <- S_crocotulus_data %>%
  expand(pos, sample)

# Merge the expanded dataframe with the rest of the data, keeping all positions
S_crocotulus <- merge(S_crocotulus_positions_expanded, S_crocotulus_data, all = T) %>%
  subset(select = c(1:3)) %>%
  inner_join(S_crocotulus_popmap) %>%
  inner_join(S_crocotulus_min_count)

# Replace NA values in the alt_count column with "3"
S_crocotulus[is.na(S_crocotulus)] <- 3

# Plot chromosome 24 heatmap - S. crocotulus
my_height = 0.7

pdf("rplots/crocotulus_chr24.pdf", width = 8.25, height = 5.875)
plt <-  S_crocotulus %>%
  ggplot(.,aes(x = as.factor(pos), y = reorder(-Freq, sex), fill = as.factor(alt_count))) + 
  geom_tile(height = my_height) +
  facet_wrap(~sex, scales = "free", ncol = 1) +
  scale_fill_manual(values = pal[c(5,5,5,9)]) +
  theme_cowplot() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())
print(plt)
dev.off()


# Plotting a comparison of contig 47 depth (Figure 3C) -------------------------

# Load in depth tables here
S_miniatus_depth_table <- read_tsv(paste0("S_miniatus/variants/S_miniatus_chr47_depth.tsv"))

S_crocotulus_depth_table <- read_tsv(paste0("S_crocotulus/variants/S_crocotulus_chr47_depth.tsv"))

# We're going to a little rearranging here and bring in the sexes
S_miniatus_m <- S_miniatus_depth_table %>%
  rename(Contig = CHROM) %>%
  pivot_longer(-c(Contig, POS), names_to = "sample", values_to = "depth") %>%
  inner_join(S_miniatus_popmap) %>%
  filter(sex == "M") %>%
  group_by(POS) %>%
  summarize(avg_depth = mean(depth)) %>%
  mutate(sex = "M")

S_miniatus_pos <- S_miniatus_m %>%
  subset(select = c(1))

S_miniatus_f <- S_miniatus_depth_table %>%
  rename(Contig = CHROM) %>%
  inner_join(S_miniatus_pos) %>%
  pivot_longer(-c(Contig, POS), names_to = "sample", values_to = "depth") %>%
  inner_join(S_miniatus_popmap) %>% 
  filter(sex == "F") %>%
  group_by(POS) %>%
  summarize(avg_depth = mean(depth)) %>%
  mutate(sex = "F")

S_miniatus <- bind_rows(S_miniatus_m, S_miniatus_f) %>%
  pivot_wider(names_from = sex, values_from = avg_depth) %>%
  mutate(., min_diff = M-F) %>%
  rename(min_M = M, min_F = F)

# Now for S. crocotulus
S_crocotulus_m <- S_crocotulus_depth_table %>%
  rename(Contig = CHROM) %>%
  pivot_longer(-c(Contig, POS), names_to = "sample", values_to = "depth") %>%
  inner_join(S_crocotulus_popmap) %>%
  filter(sex == "M") %>%
  group_by(POS) %>%
  summarize(avg_depth = mean(depth)) %>%
  mutate(sex = "M")

S_crocotulus_pos <- S_crocotulus_m %>%
  subset(select = c(1))

S_crocotulus_f <- S_crocotulus_depth_table %>%
  inner_join(S_crocotulus_pos) %>%
  rename(Contig = CHROM) %>%
  pivot_longer(-c(Contig, POS), names_to = "sample", values_to = "depth") %>%
  inner_join(S_crocotulus_popmap) %>%
  filter(sex == "F") %>%
  group_by(POS) %>%
  summarize(avg_depth = mean(depth)) %>%
  mutate(sex = "F")

S_crocotulus <- bind_rows(S_crocotulus_m, S_crocotulus_f) %>%
  pivot_wider(names_from = sex, values_from = avg_depth) %>%
  mutate(., croc_diff = M-F) %>%
  rename(croc_M = M, croc_F = F)

# Join the species frames together
croc_min <- inner_join(S_crocotulus, S_miniatus) %>%
  arrange(POS)

# Plotting figure 3C: rolling mean of depth plot

window_size = 250

# S. crocotulus
pdf("rplots/crocotulus_depth.pdf", width = 4.125, height = 1.45875)
plt <- croc_min %>%
  arrange(POS) %>%
  mutate(roll_min_diff = zoo::rollmean(min_diff, k = window_size, fill = NA)) %>%
  mutate(roll_croc_diff = zoo::rollmean(croc_diff, k = window_size, fill = NA)) %>%
  na.omit() %>%
  filter(POS >= 60000, POS <= 105000, croc_M > 0) %>%
  ggplot() +
  geom_point(aes(x = POS, 
                 y = roll_croc_diff), 
             alpha = 0.5, 
             size = 0.1,
             colour = "#EDC132") +
  theme_cowplot() +
  theme(legend.position = "none",
        panel.border = element_blank(),
        axis.line.y = element_blank(),
        text = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size = 0.2),
        strip.text = element_text(size = 6)) +
  labs(y = "Difference in average sample read depth (M-F)",
       x = "Position (bp)") +
  xlim(60000, 105000)
print(plt)
dev.off()

# S. miniatus
pdf("rplots/miniatus_depth.pdf", width = 4.125, height = 1.45875)
plt <- croc_min %>%
  arrange(POS) %>%
  mutate(roll_min_diff = zoo::rollmean(min_diff, k = window_size, fill = NA)) %>%
  mutate(roll_croc_diff = zoo::rollmean(croc_diff, k = window_size, fill = NA)) %>%
  na.omit() %>%
  filter(POS >= 60000, POS <= 105000, min_M > 0) %>%
  ggplot() +
  geom_point(aes(x = POS,
                 y = roll_min_diff),
             alpha = 0.5,
             size = 0.1,
             colour = "#ADBF5F") +
  theme_cowplot() +
  theme(legend.position = "none",
        panel.border = element_blank(),
        axis.line.y = element_blank(),
        text = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size = 0.2),
        strip.text = element_text(size = 6)) +
  labs(y = "Difference in average sample read depth (M-F)",
       x = "Position (bp)") +
  xlim(60000, 105000)
print(plt)
dev.off()