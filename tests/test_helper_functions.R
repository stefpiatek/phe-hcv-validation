library(testthat)
library(here)
library(dplyr)
library(stringr)
source(here("scripts", "helper_functions.R"))



## test select_roi() ----

input <- tibble(RefN =  c("t", "h", "i", "s", "i", "s", "_", "a","n","_", "e","x", "z"),
              Pos = c(1:13),
              Depth = c(10:22)) 

test_pattern = "an_ex"
test_region = "test_region"

test_alignment <- pairwiseAlignment(pattern = test_pattern, 
                                   subject = str_c(input$RefN, collapse = ""),
                                   type = "global-local")

output <- select_roi(test_alignment, input, test_region)


expect_equal(str_c(output$RefN, collapse = ""), test_pattern)
expect_equal(output$Pos, 8:12)
expect_equal(output$Depth, 17:21)
expect_equal(unique(output$region_name), test_region)

## test filter_to_roi() ----

# overwrite roi for test
roi <- tibble(sequence =  c("thi", "is_", "an_ex")) 

output <- filter_to_roi(input)

ns3 <- tibble(
  RefN = c("t", "h", "i"),
  Pos = 1:3,
  Depth = 10:12,
  region_name = "NS3"
)


output %>%
  ungroup() %>%
  filter(region_name == "NS3") %>%
expect_equal(ns3)

ns5a <- tibble(
  RefN = c("i", "s", "_"),
  Pos = 1:3,
  Depth = 14:16,
  region_name = "NS5A"
)

output %>%
  ungroup() %>%
  filter(region_name == "NS5A") %>%
  expect_equal(ns5a)

ns5b <- tibble(
  RefN = c("a", "n", "_", "e", "x"),
  Pos = 1:5,
  Depth = 17:21,
  region_name = "NS5B"
)

output %>%
  ungroup() %>%
  filter(region_name == "NS5B") %>%
  expect_equal(ns5b)
  

## test compensate_indel() ----

# simplified indel tbl for tests
indel_tbl <- tibble(
  sample_name = col_character(),
  Pos = col_integer(),
  RefN = col_character(),
  Gap = col_integer(),
  Depth = col_integer(),
  region_name = col_character(),
  indel_type = col_character(),
  Cons = col_character()
  )


indel_input <- tibble(
  Pos = c(1:4, 1:7),
  RefN =  c("t", "h", "i", "s", 
            "a", "c", "t", "g", "n", "a", "t"),
  Depth = c(10:20), 
  Gap = 1:11,
  sample_name = c(rep("sample", 11)),
  region_name = c(rep("one_region", 4), 
                  rep("two_region", 7)),
  Cons = c(rep("A", 11))
  )

ref_tbl <- tibble(RefN =  c("t", "h", "i", "s", 
                            "a", "c", "t","g", "a", "t"),
                  region_name = c(rep("one_region", 4), 
                                  rep("two_region", 6)),
                  Pos = c(1:10)) 
# should set n to position -1
output_indel <- compensate_indel(indel_input, ref_tbl, "two_region")


expect_equal(min(output_indel$Pos), -1)
expect_equal(nrow(output_indel), 7)

# no differences, should just be subsetted on region
output_match <- compensate_indel(indel_input, ref_tbl, "one_region")

indel_input %>%
  filter(region_name == "one_region") %>%
expect_equal(output_match)


# indel table to have one entry
expect_equal(nrow(indel_tbl), 1)

