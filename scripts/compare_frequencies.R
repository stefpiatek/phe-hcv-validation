library(readr)
library(here)
library(glue)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

# set up ggplot theme and palettes 
theme_custom <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      panel.border = element_blank(),
      legend.key = element_rect(fill = NA, colour = NA),
      panel.grid.minor.y = element_line(colour="white"),
      panel.grid.major = element_line(colour = "white"),
      axis.title.x=element_text(vjust=-1.5),
      plot.margin= unit(c(0.5, 0, 0.5, 0), "cm")# top right bottom left
    )
}
theme_set(theme_custom(base_size=12))

palette_5 <- c("#FFA347","#142952", "#B80000","#97F297", "#000000")



sample_root <- "170908_1_"





### Load data ----

consensus_fm <- read_delim(file = here("data", glue("{sample_root}quasi_frequency_matrix.txt")),
                           delim = "\t")
consensus_qb <- read_delim(file = here("data", glue("{sample_root}quasi.txt")),
                           delim = "\t") %>%
  select(-X22)  # extra column added in import

### Compare data where basecalls are discordant ----

discordant_basecalls <- consensus_qb %>%
  filter(RefN != Cons) %>%
  filter(!(Cons %in% c("R", "Y"))) 

ggplot(discordant_basecalls, aes(x=Pos, fill=Cons)) +
  geom_histogram(bins=50) +
  scale_x_continuous(limits=c(1, nrow(consensus_fm))) +
  labs(x = "Nucleotide position",
       y = "Count",
       title = "Discordant basecalls by position",
       fill = "Base called")
ggsave(glue(here("plots", "hist_discordant-calls_by-position.pdf"), device="pdf"))
# No T base miscalled, mostly just gaps and Ns

### Explore difference in base frequency ----

tmp_consensus_diff <- consensus_fm %>%
  full_join(consensus_qb, by="Pos", suffix=c("_consensus", "_quasibam")) %>%
  # keep only rows with 200 read depth -- is this too low?
  filter(Depth_quasibam >= 200) %>%
  # Keep only quasibam depth as depth
  mutate(depth = Depth_quasibam) %>%
  select(- starts_with("Depth_")) %>%
  # tidy data to have base, sample and value columns
  gather(key = base_sample, value = value, contains("_") ) %>%
  separate(col = base_sample, into = c("base", "sample"), sep = "_" ) %>%
  # get difference in percent
  group_by(Pos, base) %>%
  arrange(sample) %>%
  mutate(pc_diff = max(value) - min(value)) %>%
  # plot using log, so no negative. use factor to determine which value is greater
  mutate(consensus_greater = if_else(value < lag(value),
                                     "consensus value greater",
                                     "consensus value smaller",
                                     missing="quasibam_sample")) %>%
  mutate(consensus_greater = parse_factor(consensus_greater, 
                                          levels = c("consensus value greater",
                                                     "consensus value smaller",
                                                     "quasibam_sample"))) %>%
  ungroup()

# reshape quality data into long format and merge back into consensus_diff
quality_values <- tmp_consensus_diff %>%
  select(Pos, starts_with("q")) %>%
  gather(key = base, value = quality_value, starts_with("q")) %>%
  mutate(base = str_replace(base, "q", "")) %>%
  group_by(Pos, base) %>%
  slice(1) %>%
  ungroup()

consensus_diff <- tmp_consensus_diff %>%
  select(-starts_with("q")) %>%
  full_join(quality_values) %>%
  mutate(base = parse_factor(base, levels=c("A", "C", "G", "T", "Gap"))) %>%
  # for non-called bases, change base quality to NA
  mutate(quality_value = replace(quality_value,
                                 quality_value == 0 & value == 0,
                                 NA))

summary(consensus_diff)


### plot per base differences ----

consensus_diff %>%
  # remove samples with frequencies within 0.1% of eachother 
  filter(pc_diff >= 0.1) %>%
  filter(consensus_greater != "quasibam_sample") %>%
  ggplot(aes(x=pc_diff, fill=base)) +
  geom_histogram(bins=20) +
  scale_x_log10() + 
  facet_grid(consensus_greater ~ base) + 
  scale_fill_manual(values=palette_5) +
  labs(x = "Frequency difference",
       y = "Count",
       title = "Differences compared to consensus")
ggsave(glue(here("plots", "hist_frequency-difference.pdf"), device="pdf"))
# More gap reads in consensus, and less true base calls


consensus_diff %>%
  # remove samples with frequencies within 0.1% of eachother 
  filter(pc_diff >= 0.1) %>%
  filter(consensus_greater != "quasibam_sample") %>%
  ggplot() +
  geom_histogram(bins=20, aes(x=Pos, fill=base)) +
  facet_grid(consensus_greater ~ base) + 
  scale_fill_manual(values=palette_5) +
  labs(x = "Nucleotide position",
       y = "Count",
       title = "Differences by position")+
  theme(axis.text.x  = element_text(angle=45, vjust=0.5))
ggsave(glue(here("plots", "hist_position_frequency-difference.pdf"), device="pdf"))
# Similar positions as discodant bases, where gaps exist in consensus sequence
# As before, more gaps and fewer bases in conensus sequence


consensus_diff %>%
  filter(consensus_greater != "quasibam_sample") %>%
  filter(pc_diff >= 1) %>%
  ggplot() +
  geom_point(aes(x=Pos, y=pc_diff, colour=depth)) +
  scale_y_log10() + 
  scale_x_continuous(limits=c(1, nrow(consensus_fm))) + 
  facet_grid(consensus_greater ~ base) +
  scale_color_continuous(low = palette_5[1] , high = palette_5[2]) +
  labs(x = "Nucleotide position",
       y = "Frequency difference",
       title = "Percentage difference per base, over 1% with read depth") +
  theme(axis.text.x  = element_text(angle=45, vjust=0.5))
ggsave(glue(here("plots", "dot_frequency-difference_with-depth.pdf"), device="pdf"))
# Read depth of less than 1,500 seems to affect bases at ~ 1,000 nt, otherwise read depth at 2,000


consensus_diff %>%
  filter(consensus_greater != "quasibam_sample") %>%
  filter(pc_diff >= 1) %>%
  ggplot() +
  geom_point(aes(x=Pos, y=pc_diff, colour=quality_value)) +
  scale_y_log10() + 
  scale_x_continuous(limits=c(1, nrow(consensus_fm))) + 
  facet_grid(consensus_greater ~ base) +
  scale_color_continuous(low = palette_5[1] , high = palette_5[3]) +
  labs(x = "Nucleotide position",
       y = "Frequency difference",
       title = "Percentage difference per base, over 1% with mapping quality") +
  theme(axis.text.x  = element_text(angle=45, vjust=0.5))
ggsave(glue(here("plots", "dot_frequency-difference_with-quality.pdf"), device="pdf"))
# small proportion of consensus bases have no read depth (i.e. quality of NA) 
