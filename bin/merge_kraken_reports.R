#!/usr/bin/env Rscript

###############################################################################
# AUTHOR : OLABIYI ADEREMI OBAYOMI
# DESCRIPTION: A script to merge several kraken reports in a directory into one species table.
# E-mail: obadbotanist@yahoo.com
# Created: September 2025
# example: Rscript merged_kraken_reports.R --reports-dir '/path/to/kraken_reports'
###############################################################################

library(optparse)

######## -------- Get input variables from the command line ----##############

version <- 1.0 

# Input options 
option_list <- list(
  
  make_option(c("-r", "--reports-dir"), type="character", default="./", 
              help="A directory containing kraken2 report files",
              metavar="path"),
  
  make_option(c("-o", "--output-prefix"), type="character", default="", 
              help="Unique name to tag onto output files. Default: empty string.",
              metavar=""),
  
  make_option(c("-a", "--assay-suffix"), type="character", default="_GLmetagenomics", 
              help="Genelab assay suffix.", metavar="_GLmetagenomics"),
  
  make_option(c("--version"), action = "store_true", type="logical", 
              default=FALSE,
              help="Print out version number and exit.", metavar = "boolean")
)




opt_parser <- OptionParser(
  option_list=option_list,
  usage = "Rscript %prog \\
                  --reports-dir 'data/read_taxonomy/kraken2/kraken2/'",
  description = paste("Author: Olabiyi Aderemi Obayomi",
                      "\nEmail: olabiyi.a.obayomi@nasa.gov",
                      "\n A script to merge multiple kraken reports into one species table.",
                      "\nIt outputs a merge kraken species table ",
                      sep="")
)


opt <- parse_args(opt_parser)

# print(opt)
# stop()


if (opt$version) {
  cat("merge_kraken_reports.R version: ", version, "\n")
  options_tmp <- options(show.error.messages=FALSE)
  on.exit(options(options_tmp))
  stop()
}



library(tidyverse)
library(pavian)
library(glue)

reports_dir <- opt[["reports-dir"]]

reports <- read_reports(reports_dir)

samples <- names(reports) %>% str_split("-") %>% map_chr(function(x) pluck(x, 1))
merged_reports  <- merge_reports2(reports, col_names = samples)
taxonReads <- merged_reports$taxonReads
cladeReads <- merged_reports$cladeReads
tax_data <- merged_reports[["tax_data"]]

species_table <- tax_data %>% 
  bind_cols(cladeReads) %>%
  filter(taxRank %in% c("U","S")) %>% 
  select(-contains("tax")) %>%
  zero_if_na() %>% 
  filter(name != 0) %>%  # drop unknown taxonomies
  group_by(name) %>% 
  summarise(across(everything(), sum)) %>% 
  ungroup() %>% 
  as.data.frame() %>% 
  rename(species=name)

species_names <- species_table[,"species"]
rownames(species_table) <- species_names

output_prefix <- opt[["output-prefix"]]
assay_suffix <- opt[["assay-suffix"]]
write_tsv(x = species_table, 
          file = glue("{output_prefix}merged_kraken_table{assay_suffix}.tsv"))
