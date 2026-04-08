#!/usr/bin/env Rscript


###############################################################################
# AUTHOR : OLABIYI ADEREMI OBAYOMI
# DESCRIPTION: A script to process humann function tables generated from read based processing.
# E-mail: obadbotanist@yahoo.com
# Created: March 2026
# example: Rscript process_humann_table.R \
#                  --table 'Gene-families-cpm.tsv' \
#                  --type 'uniref'

###############################################################################

library(optparse)



######## -------- Get input variables from the command line ----##############

version <- 1.0 


option_list <- list(
  
  make_option(c("-t", "--table"), type="character",
              default=NULL, 
              help="path to a tab separated contig or gene level taxonomy or functions(KO) table. 
              i.e Gene-families-cpm.tsv or Gene-families-KO-cpm.tsv or Pathway-abundances-cpm.tsv",
              metavar="path"),
  make_option(c("-T", "--type"), type="character", default="taxonomy", 
              help="The type of assembly annotation performed in the input table.
              either  uniref, KO and pathway. Default: pathway.",
              metavar="type"),
  make_option(c("-o", "--output-prefix"), type="character", default="", 
              help="Unique name to tag onto output files. Default: empty string.",
              metavar=""),
  make_option(c("-a", "--assay-suffix"), type="character", default="_GLMetagenomics", 
              help="Genelab assay suffix.", metavar="GLMetagenomics"),
  make_option(c("--version"), action = "store_true", type="logical", 
              default=FALSE,
              help="Print out version number and exit.", metavar = "boolean")
)


library(tibble)
library(tidyr)
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(magrittr)
library(glue)

opt_parser <- OptionParser(
  option_list=option_list,
  usage = "Rscript %prog \\
                  --table 'Pathway-abundances-cpm.tsv'
                  --type 'pathway'" ,
  description = paste("Author: Olabiyi Aderemi Obayomi",
                      "\nEmail: olabiyi.a.obayomi@nasa.gov",
                      "\n  A script to reformat humann functions table.",
                      "\nIt outputs a reformatted uniref, KO or pathway table.",
                      sep="")
)


opt <- parse_args(opt_parser)

# print(opt)
# stop()


if (opt$version) {
  cat("process_humann_table.R version: ", version, "\n")
  options_tmp <- options(show.error.messages=FALSE)
  on.exit(options(options_tmp))
  stop()
}



if(is.null(opt[["table"]])) {
  stop("Path to a tab separated humann table file must be set.")
}

if(is.null(opt[["type"]])) {
  stop("You must provide a valid input table type. One of 'KO', 'uniref', and 'pathway'.")
}



humann_table <- opt[["table"]]
type <- opt[["type"]]
prefix <- opt[["output-prefix"]]
suffix <- opt[["assay-suffix"]] 



feature_table <- read_delim(file = humann_table, delim = "\t")


if(type == "KO"){
  
  table2write <- feature_table %>% 
    rename(KO=`# Gene Family`) %>% 
    set_names(colnames(.) %>% str_replace_all("_Abundance-CPM", "")) %>%
    as.data.frame()
  
  
  write_tsv(x = table2write, file = glue("{prefix}Gene-families-KO{suffix}.tsv"))
  
}else if(type == "uniref"){
  
  table2write <-  feature_table  %>% 
    rename(Uniref90=`# Gene Family`) %>%
    mutate(Uniref90=str_replace_all(Uniref90, "UniRef90_", "")) %>%
    set_names(colnames(.) %>% str_replace_all("_Abundance-CPM", "")) %>%
    as.data.frame()
  
   write_tsv(x = table2write, file = glue("{prefix}Gene-families-uniref{suffix}.tsv"))
  
}else if(type == "pathway"){
  
  # Pathway
  table2write <- feature_table  %>% 
    rename(Pathway=`# Pathway`) %>%
    set_names(colnames(.) %>% str_replace_all("_Abundance-CPM", "")) %>%
    as.data.frame()
  
  write_tsv(x = table2write, file = glue("{prefix}Pathway-abundances{suffix}.tsv"))
}else{
  
  stop("You must provide a valid input table type. One of 'KO', 'uniref', and 'pathway'.")
  
}

