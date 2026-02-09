#!/usr/bin/env Rscript


###############################################################################
# AUTHOR : OLABIYI ADEREMI OBAYOMI
# DESCRIPTION: A script to perform feature decontamination with decontam.
# E-mail: obadbotanist@yahoo.com
# Created: November 2025
# example: Rscript process_metaphlan.R \
#                  --metaphlan-table 'Metaphlan-taxonomy_GLmetagenomics.tsv' \
#                  --read-count 'reads_per_sample.txt' \
#                  --assay-suffix '_GLlbnMetag'
###############################################################################

library(optparse)



######## -------- Get input variables from the command line ----##############

version <- 1.0 



option_list <- list(
  
  make_option(c("-m", "--metaphlan-table"), type="character",
              default='Metaphlan-taxonomy_GLmetagenomics.tsv', 
              help="path to a metaphlan taxonomy table 
              i.e the output of running metaphlan on the commandline.",
              metavar="path"),
  make_option(c("-r", "--read-count"), type="character",
              default='reads_per_sample.txt', 
              help="path to a 2-column tab delimited file with 
              sample names and read counts as columns 1 and 2, respectively.",
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
                  --merged-table 'process_metaphalan_table.tsv' " ,
  description = paste("Author: Olabiyi Aderemi Obayomi",
                      "\nEmail: olabiyi.a.obayomi@nasa.gov",
                      "\n  A script to convert a raw metaphlan table into a conventional species table.",
                      "\nIt outputs a conventional species table with the file column as species names \
                      and other columns as per sample species count.",
                      sep="")
)


opt <- parse_args(opt_parser)

# print(opt)
# stop()


if (opt$version) {
  cat("process_metaphlan_table.R version: ", version, "\n")
  options_tmp <- options(show.error.messages=FALSE)
  on.exit(options(options_tmp))
  stop()
}



if(is.null(opt[["metaphlan-table"]])) {
  stop("Path to a metaphlan table file must be set.")
}

if(is.null(opt[["read-count"]])) {
  stop("Path to a reads count per sample file must be set.")
}

# Function to read an input table into a dataframe
read_input_table <- function(file_name){
  
  df <- read_delim(file = file_name, comment = "#")
  colnames(df)[1] <- "taxonomy"
  
  df <- df %>%
    filter(str_detect(taxonomy, "UNCLASSIFIED|s__") & 
             str_detect(taxonomy, "t__", negate = TRUE))
  return(df)
}


# Read read-based species table
read_species_table <- function(df) {
  
  taxon_levels <- c("Kingdom", "Phylum", "Class", "Order",
                    "Family", "Genus", "Species")
  
  
  df <-  df %>%
    mutate(Species=str_replace_all(taxonomy, '\\w__', "")) %>% 
    separate(Species, into=taxon_levels, sep="\\|") %>%
    mutate(across(where(is.character),
                  function(x) replace_na(x, "UNCLASSIFIED")
    )
    ) %>% 
    mutate(Species=str_replace_all(Species, "_", " ")) %>% 
    select(-taxonomy, -Kingdom, -Phylum, -Class, -Order, -Family, -Genus) %>% 
    select(Species, everything())
  
  
  return(df)
}

species_file <- opt[["metaphlan-table"]]
read_count <- opt[["read-count"]]
prefix <- opt[["output-prefix"]]
suffix <- opt[["assay-suffix"]]

# Read table
abund_species_table <-  read_input_table(species_file) %>%
  read_species_table() %>% as.data.frame()
rownames(abund_species_table) <- abund_species_table$Species
abund_species_table <- abund_species_table[,-match("Species", colnames(abund_species_table))]

# Make max abundance equal to 1
tab2 <- (abund_species_table %>% t) / 100


# Read raw read count by sample
counts <- read_delim(read_count, col_names = c("Sample_ID", "Reads"), skip = 1) %>%
         	mutate(Sample_ID=map_chr(Sample_ID, function(val) { 
						  str_replace_all(val,  "(.+)_R[12].+", "\\1" ) 
                                                   } ) ) %>%
                distinct()

max_val <- counts$Reads %>% max
# If counts are expressed as decimals
# multiply by a million
if(max_val < 10) { counts <- counts %>%  mutate(Reads= Reads * 1000000)}
counts <- counts %>% as.data.frame()


# Set rownames as sample names
rownames(counts) <- counts$Sample_ID
# Drop the Sample_ID column
counts <- counts[, -1, drop=FALSE]

samples <- intersect(rownames(counts), rownames(tab2))
tab2 <- tab2[samples,]
counts <- counts[samples, ,drop=FALSE]

# Convert relative abundance to raw count
species_table <- map2(tab2 %>% as.data.frame, 
                      colnames(tab2), function(col, specie) {
                        df <- col * counts
                        colnames(df) <- specie
                        return(df) 
                      }) %>% list_cbind() %>% t


table2write <- species_table  %>%
  as.data.frame() %>%
  rownames_to_column("Species")

write_tsv(x = table2write, file = glue("{prefix}metaphlan_species_table{suffix}.tsv"))
