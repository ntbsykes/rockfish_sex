# Setup ------------------------------------------------------------------------

# Install required packages
if(!require(tidyverse)){
  install.packages("tidyverse")
  library(tidyverse)
}
if(!require(cowplot)){
  install.packages("cowplot")
  library(cowplot)
}
if(!require(zoo)){
  install.packages("zoo")
  library(zoo)
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

# Set global p cutoff
pvalue_cutoff = 0.005

# Make a species list and chromosome list vector
Sebastes_list <- c("carnatus", "chrysomelas", "crocotulus", "diaconus", 
                   "miniatus", "paucispinis", "pinniger", "ruberrimus")
chr_list <- c("Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", 
              "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", 
              "Chr15", "Chr16", "Chr17", "Chr18", "Chr19", "Chr20", "Chr21", 
              "Chr22", "Chr23", "Chr24")

# Load in contig/chromosome map
chr_map <- read_tsv("rockfish_chromosomes.tsv", col_names = c("Contig", "chr_n"))


# Load in the data -------------------------------------------------------------

# Make a tibble for assembly of all the species information
Sebastes <- tibble()

# This loop loads in all the species data one at a time and builds our dataframe
for (spp in Sebastes_list){
  # Load in allele frequency data and test significance
  spp_allelefreq <- read_table(paste0("S_",spp,"/variants/S_",spp,".allelefreq.tsv"))
  spp_allelefreq <- mutate(spp_allelefreq, p = pchisq(spp_allelefreq$chi_square, 1, 
                                                      lower.tail = F))
  spp_allelefreq <- mutate(spp_allelefreq, signif = case_when((p) < pvalue_cutoff ~ TRUE,
                                                              TRUE ~ FALSE))
  spp_allelefreq <- mutate(spp_allelefreq, stat = "allelefreq")
  
  # Load in the heterozygosity data and test significance
  spp_hetdif <- read_table(paste0("S_",spp,"/variants/S_",spp,".hetdif.tsv"))
  spp_hetdif <- mutate(spp_hetdif, p = pchisq(spp_hetdif$chi_square, 1, lower.tail = F))
  spp_hetdif <- mutate(spp_hetdif, signif = case_when((p) < pvalue_cutoff ~ TRUE,
                                                      TRUE ~ FALSE))
  spp_hetdif <- mutate(spp_hetdif, stat = "heterozygosity")
  
  # Load in missing sex data and test significance
  spp_missingsex <- read_table(paste0("S_",spp,"/variants/S_",spp,".missingsex.tsv"))
  spp_missingsex <- mutate(spp_missingsex, p = pchisq(spp_missingsex$chi_square, 1, 
                                                      lower.tail = F))
  spp_missingsex <- mutate(spp_missingsex, signif = case_when((p) < pvalue_cutoff & bias > 0 ~ TRUE,
                                                              TRUE ~ FALSE))
  spp_missingsex <- mutate(spp_missingsex, stat = "missingdata")
  
  # Load in RADSex data and edit columns
  spp_RADsex <- read_table(paste0("S_",spp,"/radsex/S_",spp,".radsex.tsv"), 
                           comment = "#", col_names = TRUE)
  spp_RADsex <- mutate(spp_RADsex, stat = "RADsex")
  spp_RADsex <- rename(spp_RADsex, p = P, signif = Signif, bias = Bias)
  
  # Now bind them all together into a single data frame & label with species name
  sp <- bind_rows(spp_allelefreq, spp_hetdif, spp_missingsex, spp_RADsex) %>%
    subset(., select = c(1,2,8,9,10))
  sp <- mutate(sp, species = paste0("S_",spp,""))
  Sebastes <- rbind(Sebastes, sp)
}

Sebastes <- inner_join(Sebastes, chr_map)


# Plotting the GWAS: Setup -----------------------------------------------------

# Load in chromosome lengths dataframe, and isolate desired data from Sebastes
chr_lengths <- read_tsv("S_diaconus/S_diaconus.genome.fasta.fai",
                        col_names = c("Contig", "length", "bits", 
                                      "spacer1", "spacer2")) %>%
  inner_join(chr_map) %>%
  filter(chr_n != "NA") %>%
  arrange(chr_n)

spp_phylo <- read_tsv("spp_phylogeny.tsv", col_names = TRUE)

all_data <- Sebastes %>%
  filter(p != "NA") %>%
  inner_join(spp_phylo)

# Make a new dataframe that includes cumulative position across the genome
chr_lengths %>%
  arrange(chr_n) %>%
  select(chr_n, length) %>%
  mutate(total = cumsum(length)-length) %>%
  select(-length) %>%
  left_join(all_data, ., by = c("chr_n" = "chr_n")) %>%
  arrange(chr_n, Position) %>%
  mutate(cumulative_pos = Position + total) -> all_data_cumulative

# Make a dataframe for the background colours
chr_breaks <- tibble(start = unique(all_data_cumulative$total),
                     end = unique(all_data_cumulative$total)+(chr_lengths$length)-1,
                     colours = c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1))

# Make a new dataframe with chromosome midpoint for labels
all_axis_df <- all_data_cumulative %>%
  group_by(chr_n) %>%
  summarize(centre = (max(cumulative_pos) + min(cumulative_pos))/2)


# Isolate the most telling data from each species ------------------------------

# S. carnatus & S.chrysomelas
target_spp <- c("S_carnatus", "S_chrysomelas")
target_stat <- c("allelefreq", "heterozygosity", "RADsex")
car_chr <- all_data_cumulative %>%
  filter(species %in% target_spp, stat %in% target_stat) %>%
  subset(select = c(3,4,6,7,8,9,10))

# S. crocotulus & S.miniatus
target_spp <- c("S_crocotulus", "S_miniatus")
target_stat <- c("missingdata","allelefreq", "heterozygosity")
cro_min <- all_data_cumulative %>%
  filter(species %in% target_spp, stat %in% target_stat) %>%
  subset(select = c(3,4,6,7,8,9,10))

# S. paucispinis, S. pinniger, S.ruberrimus
target_spp <- c("S_paucispinis", "S_pinniger", "S_ruberrimus")
target_stat <- c("RADsex")
pau_pin_rub <- all_data_cumulative %>%
  filter(species %in% target_spp, stat %in% target_stat) %>%
  subset(select = c(3,4,6,7,8,9,10))

# S. diaconus
dia <- all_data_cumulative %>%
  filter(species == "S_diaconus", stat == "RADsex") %>%
  subset(select = c(3,4,6,7,8,9,10))

filt_data <- rbind(car_chr, cro_min, pau_pin_rub, dia)
  

# Plot the whole genome rolling mean for each species (Figure 1) ---------------

# Select window size
window_size <- 250

# Plot it
pdf("rplots/line.pdf", width = 8.25, height = 11.75)
plt <- filt_data %>%
  group_by(species, chr_n) %>%
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
  scale_x_continuous(label = all_axis_df$chr_n, 
                     breaks = all_axis_df$centre) +
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
  labs(y = "Proportion of significant sites per window", 
       color = "Species") +
  facet_wrap(~top, ncol = 1)
print(plt)
dev.off()
