#!/usr/bin/env Rscript


###############################################################################
# AUTHOR : OLABIYI ADEREMI OBAYOMI
# DESCRIPTION: A script to convert a raw kaiju table into a conventional species table.
# E-mail: obadbotanist@yahoo.com
# Created: November 2025
# example: Rscript process_kaiju_table.R --merged-table 'merged_kaiju_table.tsv'
###############################################################################

library(optparse)



######## -------- Get input variables from the command line ----##############

version <- 1.0 


option_list <- list(
  
  make_option(c("-m", "--merged-table"), type="character",
              default='merged_kaiju_table.tsv', 
              help="path to a merged kaiju table (merged_kaiju_table.tsv) 
              i.e the output of running kaiju2table.",
              metavar="path"),
  make_option(c("-o", "--output-prefix"), type="character", default="", 
              help="Unique name to tag onto output files. Default: empty string.",
              metavar=""),
  make_option(c("-a", "--assay-suffix"), type="character", default="_GLMetagenomics", 
              help="Genelab assay suffix.", metavar="GLMetagenomics"),
  make_option(c("--version"), action = "store_true", type="logical", 
              default=FALSE,
              help="Print out version number and exit.", metavar = "boolean")
)




library(tidyverse)
library(glue)




opt_parser <- OptionParser(
  option_list=option_list,
  usage = "Rscript %prog \\
                  --merged-table 'merged_kaiju_table.tsv' " ,
  description = paste("Author: Olabiyi Aderemi Obayomi",
                      "\nEmail: olabiyi.a.obayomi@nasa.gov",
                      "\n  A script to convert a raw kaiju table into a conventional species table.",
                      "\nIt outputs a conventional tab separated species table.",
                      sep="")
)


opt <- parse_args(opt_parser)

# print(opt)
# stop()


if (opt$version) {
  cat("process_kaiju_table.R version: ", version, "\n")
  options_tmp <- options(show.error.messages=FALSE)
  on.exit(options(options_tmp))
  stop()
}



if(is.null(opt[["merged-table"]])) {
  stop("Path to a merged kaiju table file must be set.")
}


# ------- Functions ------------------- #
get_last_assignment <- function(taxonomy_string, split_by =';', remove_prefix = NULL){
  # A function to get the last taxonomy assignment from a taxonomy string 
  split_names <- strsplit(x =  taxonomy_string , split = split_by) %>% 
    unlist()
  
  level_name <- split_names[[length(split_names)]]
  
  if(level_name == "_"){
    return(taxonomy_string)
  }
  
  if(!is.null(remove_prefix)){
    level_name <- gsub(pattern = remove_prefix, replacement = '', x = level_name)
  }
  
  return(level_name)
}

mutate_taxonomy <- function(df, taxonomy_column="taxonomy"){
  
  # make sure that the taxonomy column is always named taxonomy
  col_index <- which(colnames(df) == taxonomy_column)
  colnames(df)[col_index] <- 'taxonomy'
  df <- df %>% dplyr::mutate(across( where(is.numeric), \(x) tidyr::replace_na(x,0)  ) )%>% 
    dplyr::mutate(taxonomy=map_chr(taxonomy,.f = function(taxon_name=.x){
      last_assignment <- get_last_assignment(taxon_name) 
      last_assignment  <- gsub(pattern = "\\[|\\]|'", replacement = '',x = last_assignment)
      trimws(last_assignment, which = "both")
    })) %>% 
    as.data.frame(check.names=FALSE, StringAsFactor=FASLE)
  # Ensure the taxonomy names are unique by aggregating duplicates
  df <- aggregate(.~taxonomy,data = df, FUN = sum)
  return(df)
}


process_kaiju_table <- function(file_path, taxon_col="taxon_name") {
  
  kaiju_table <- read_delim(file = file_path,
                            delim = "\t",
                            col_names = TRUE) 
  
  # Create  a sample colname if the file column wasn't pre-edited
  if(colnames(kaiju_table)[1] ==  "file" ){
    
    kaiju_table <-  kaiju_table %>%
      filter(!str_detect(file, "dmp")) %>%
      mutate(file=str_replace_all(file, ".+/(.+)_kaiju.out", "\\1")) %>%
      rename(sample=file)
  }
  
  abs_abun_df <- kaiju_table %>% # read input table
    select(sample, reads, taxonomy=!!sym(taxon_col)) %>%
    pivot_wider(names_from = "sample", values_from = "reads", 
                names_sort = TRUE) %>% # convert long dataframe to wide dataframe
    mutate_taxonomy # mutate the taxonomy column such that it contains only lowest taxonomy assignment
  
  # Set the taxon names as row names, drop the taxonomy column and convert to a matrix
  rownames(abs_abun_df) <- abs_abun_df[,"taxonomy"]
  abs_abun_df <- abs_abun_df[,-(which(colnames(abs_abun_df) == "taxonomy"))]
  # Drop unwanted columns if present
  unwanted_columns <- grep("\\.dmp$", colnames(abs_abun_df))
  if(length(unwanted_columns) > 0) { abs_abun_df <- abs_abun_df[,-unwanted_columns] }
  abs_abun_matrix <- as.matrix(abs_abun_df)
  
  return(abs_abun_matrix)
}


table_path <- opt[["merged-table"]]
prefix <- opt[["output-prefix"]]
suffix <- opt[["assay-suffix"]]

feature_table <- process_kaiju_table(file_path=table_path)
table2write <- feature_table  %>%
  as.data.frame() %>%
  rownames_to_column("Species")

write_tsv(x = table2write, file = glue("{prefix}kaiju_species_table{suffix}.tsv"))



