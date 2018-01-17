## set up regions of interest ----
#+/- 30 bp, unless this would cause overlap

hcv_fasta <- read_delim(file = here("data", "reference",  "hcv1.fas"), 
                        delim = "\t", col_names = "hcv1") %>%
  filter(!startsWith(hcv1, ">")) 
hcv1_seq <- str_c(hcv_fasta$hcv1, collapse="")

padding = 30
roi <- tibble(region_name =  c("NS3", "NS5A",        "NS5B"),
              start = c(        3420,  6258,           7602 + padding),
              end =   c(        5312,  7601 - padding, 9374)
) %>%
  mutate(start = start - padding) %>%
  mutate(end = end + padding) %>%
  mutate(width = end - start + 1) %>%
  group_by(region_name) %>%
  mutate(sequence = str_sub(hcv1_seq, start = start, end = end)) 

## functions to filter to regions of interest ----
select_roi <- function(roi_alignment = NULL, unfiltered_file = NULL, name = ""){
  subject_alignment = subject(roi_alignment)
  unfiltered_file %>% 
    filter(between(Pos, start(subject_alignment), end(subject_alignment))) %>%
    mutate(region_name = name) %>%
    return()
}

filter_to_roi <- function(unfiltered_file = NULL){
  roi_alignment <- pairwiseAlignment(pattern = roi$sequence, 
                                     subject = str_c(unfiltered_file$RefN, collapse = ""),
                                     type = "subjectOverlap")
  select_roi(roi_alignment[1], unfiltered_file, name = "NS3") %>%
    bind_rows(select_roi(roi_alignment[2], unfiltered_file, name = "NS5A"))%>%
    bind_rows(select_roi(roi_alignment[3], unfiltered_file, name = "NS5B")) %>%
    group_by(region_name) %>%
    arrange(Pos) %>%
    mutate(Pos = 1:n()) %>%
    return()
}