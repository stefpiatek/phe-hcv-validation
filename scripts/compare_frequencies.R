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


### Compare consensus frequency matrix and quasibam output ----

consensus_fm <- read_delim(file = here("data", glue("{sample_root}quasi_frequency_matrix.txt")),
                           delim = "\t")
consensus_qb <- read_delim(file = here("data", glue("{sample_root}quasi.txt")),
                           delim = "\t") %>%
  select(-X22)  # extra column added in import

consensus_diff_temp <- consensus_fm %>%
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
  mutate(consensus_greater = if_else(value < lag(value),
                                     "consensus value greater",
                                     "consensus value smaller",
                                     missing="quasibam_sample")) %>%
  ungroup()

# reshape quality data into long format and merge back into consensus_diff
quality_values <- consensus_diff_temp %>%
  select(Pos, starts_with("q")) %>%
  gather(key = base, value = quality_value, starts_with("q")) %>%
  mutate(base = str_replace(base, "q", "")) %>%
  group_by(Pos, base) %>%
  slice(1)

consensus_diff <- consensus_diff_temp %>%
  select(-starts_with("q")) %>%
  full_join(quality_values) %>%
  mutate(base = parse_factor(base, levels=c("A", "C", "G", "T", "Gap")))


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


consensus_diff %>%
  filter(consensus_greater != "quasibam_sample") %>%
  filter(pc_diff >= 1) %>%
  ggplot() +
  geom_point(aes(x=Pos, y=pc_diff, colour=depth)) +
  scale_y_log10() + 
  scale_x_continuous(limits=c(1, nrow(consensus_fm))) + 
  facet_grid(consensus_greater ~ base) +
  scale_color_continuous(low = palette_5[1] , high = palette_5[3]) +
  labs(x = "Nucleotide position",
       y = "Frequency difference",
       title = "Percentage difference per base, over 1% with read depth") +
  theme(axis.text.x  = element_text(angle=45, vjust=0.5))


consensus_diff %>%
  filter(consensus_greater != "quasibam_sample") %>%
  filter(pc_diff >= 1) %>%
  ggplot() +
  geom_point(aes(x=Pos, y=pc_diff, colour=quality_value)) +
  scale_y_log10() + 
  scale_x_continuous(limits=c(1, nrow(consensus_fm))) + 
  facet_grid(consensus_greater ~ base) +
  scale_color_continuous(low = palette_5[1] , high = palette_5[2]) +
  labs(x = "Nucleotide position",
       y = "Frequency difference",
       title = "Percentage difference per base, over 1% with mapping quality") +
  theme(axis.text.x  = element_text(angle=45, vjust=0.5))


# Add consensus sequence to frequency matrix?
# 