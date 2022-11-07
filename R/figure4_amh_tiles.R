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

# Set up colour palette
pal <- pnw_palette("Bay", 8) %>%
  append("grey90")

# Loading in the data ----------------------------------------------------------

amh_deletions <- read_tsv("amh_deletions.tsv")

amh_phylo <- read_tsv("amh_phylogeny.tsv")

amh <- inner_join(amh_deletions, amh_phylo)

amh <- amh %>%
  pivot_longer(!(Species | top), names_to = "del", values_to =  "pres")


# Plotting the deletions (Figure 4B) -------------------------------------------

pdf("rplots/amh_phylogeny_tiles_test.pdf", paper = "a4", 
    width = 8.25, height = 11.75)
plt <- amh %>%
  group_by(top) %>%
  mutate(del = recode(del, x131 = "131bp", x166 = "166bp")) %>%
  ggplot(aes(del, -top)) +
  geom_tile(aes(fill = as.factor(pres)), color = "white") +
  scale_fill_manual(values = pal[c(1,5)]) +
  theme_cowplot() +
  labs(x = "Deletions") +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_blank())
print(plt)
dev.off()
  
