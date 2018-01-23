library(Biostrings)  # bioconductor package
library(here)
library(dplyr)
library(stringr)
library(readr)

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
                                     type = "global-local")
  select_roi(roi_alignment[1], unfiltered_file, name = "NS3") %>%
    bind_rows(select_roi(roi_alignment[2], unfiltered_file, name = "NS5A"))%>%
    bind_rows(select_roi(roi_alignment[3], unfiltered_file, name = "NS5B")) %>%
    group_by(region_name) %>%
    arrange(Pos) %>%
    mutate(Pos = 1:n()) %>%
    return()
}

## Deal with indels: currently only insertions in pipeline quasibams are dealt with ----

# create indel_tbl to keep track of all indels
indel_tbl <- tibble(
  sample_name = col_character(),
  Pos = col_integer(),
  RefN = col_character(),
  A = col_integer(),
  C = col_integer(),
  T = col_integer(),
  G = col_integer(),
  Gap = col_integer(),
  Depth = col_integer(),
  region_name = col_character(),
  indel_type = col_character())



### Compensate for an individual region indels ----
compensate_indel <- function(indel_table = NULL, ref_table = NULL, region = NULL){
  # filter both tables to region of interest only
  roi_indel_table <- indel_table %>%
    filter(region_name == region)
  
  tmp_ref <- ref_table %>%
    filter(region_name == region)
  
  ref_sequence <- str_c(tmp_ref$RefN, collapse="")
  subject_sequence <- str_c(roi_indel_table$RefN, collapse = "")
  
  # align sequences
  roi_alignment <- pairwiseAlignment(pattern = ref_sequence, 
                                     subject = subject_sequence,
                                     type = "global-local")
  pattern_alignment = pattern(roi_alignment)
  # note: del here is a deletion compared to reference, i.e. to delete from indel_table
  del_start <- start(deletion(roi_alignment)[[1]])
  del_end <- end(deletion(roi_alignment)[[1]])
  
  # If there are bases to delete from indel_table
  if(length(del_end > 0)){
    for(del_i in 1:length(del_start)){
      # Add information to indel table (global variable)
      indel_tbl <<- roi_indel_table %>%
        filter(between(Pos, del_start[del_i], del_end[del_i])) %>%
        select(Pos:Gap, sample_name:region_name, Cons) %>%
        mutate(indel_type = "insertion") %>%
        bind_rows(indel_tbl)
      
      roi_indel_table <- roi_indel_table %>%
        mutate(Pos = if_else(between(Pos, del_start[del_i], del_end[del_i]),
                             -1,
                             as.double(Pos)))
    }
  }
  
  return(roi_indel_table)  
}


### Run all compensations for indels -----
compensate_for_indels <- function(indel_table = NULL, sample = "170908_18", ref_table = NULL){
  # filter reference table to only the sample of interest
  ref_table <- ref_table %>%
    filter(sample_name == sample)
  
  output_table <- compensate_indel(indel_table, ref_table, "NS3") %>%
    bind_rows(compensate_indel(indel_table, ref_table, "NS5A")) %>%
    bind_rows(compensate_indel(indel_table, ref_table, "NS5B"))
  
  return(output_table)
}

## printing summary information for filtered differences ----

analyse_threshold <- function(quality_threshold = 1, source_name = "consensus"){
  pass_filter <<- consensus_diff %>%
    filter(pc_diff >= quality_threshold) %>%
    filter(quality_value > 1) %>%
    filter(source == source_name) 
  
  print(glue("Positions with >= {quality_threshold}% difference"))
  all_positions <- pass_filter %>%
    group_by(Pos, region_name) %>%
    filter(Pos == max(Pos)) %>%
    group_by(region_name) %>%
    summarise(positions = n()) 
  
  position_diff_export <<- pass_filter %>%
    filter(!sample_name %in% c('171009_10', '171009_11', '171009_12')) %>%
    group_by(Pos, region_name) %>%
    filter(Pos == max(Pos)) %>%
    group_by(region_name) %>%
    summarise(filtered_positions = n()) %>%
    full_join(all_positions, by = "region_name") 
  
  
  bases_pass <- pass_filter %>%
    filter(base != "Gap") %>%
    group_by(region_name) %>%
    summarise(bases = n())
  
  conf_int <- tidy(prop.test(sum(bases_pass$bases), 
                             sum(total_bases$max_length))) %>%
    select(conf.low:conf.high) %>%
    mutate(filtering = "All bases") %>%
    select(filtering, conf.low:conf.high)
  
  bases_pass_filtered <- pass_filter %>%
    filter(!sample_name %in% c('171009_10', '171009_11', '171009_12')) %>%
    filter(base != "Gap") %>%
    group_by(region_name) %>%
    summarise(bases = n())
  
  conf_int_export <<- tidy(prop.test(sum(bases_pass_filtered$bases), 
                                     sum(total_bases_filtered$max_length))) %>%
    select(conf.low:conf.high) %>%
    mutate(filtering = "remove 171009_10:12") %>%
    select(filtering, conf.low:conf.high) %>%
    bind_rows(conf_int) 
}

