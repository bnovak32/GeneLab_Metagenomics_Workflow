#!/usr/bin/env Rscript


###############################################################################
# AUTHOR : OLABIYI ADEREMI OBAYOMI
# DESCRIPTION: A script to perform feature decontamination with decontam.
# E-mail: obadbotanist@yahoo.com
# Created: January 2026
# example: Rscript run_decontam.R \
#                  --feature-table 'kaiju_species_table_GLlbnMetag.csv' \
#                  --feature-column 'Species' \
#                  --metadata-table 'mapping/metadata.csv' \
#                  --samples-column 'Sample_ID' \
#                  --prevalence-column 'NTC' \
#                  --frequency-column 'concentration' \
#                  --ntc-name  'true' \
#                  --threshold 0.1 \
#                  --classification-method 'kaiju' \
#                  --assay-suffix '_GLlbnMetag' 
###############################################################################

library(optparse)



######## -------- Get input variables from the command line ----##############

version <- 1.0 

# Input options 
option_list <- list(


  make_option(c("-i", "--feature-table"), type="character", default=NULL, 
              help="path to a tab or comma separated samples feature table 
              i.e. species/functions table with species/functions as the first column 
              and samples as other columns.",
              metavar="path"),
    
  make_option(c("-f", "--feature-column"), type="character", default=NULL, 
              help="Feature column name in feature table ['Species', 'species', 'KO_ID'].
              Default: empty string",
              metavar="Feature_Column"),
 
  make_option(c("-m", "--metadata-table"), type="character", default=NULL, 
              help="path to a comma separated samples metadata file with the 
               --prevalence-column and/or --frequency-column to be analyzed.",
              metavar="path"),
  
  make_option(c("-s", "--samples-column"), type="character", default="Sample Name", 
              help="Column in metadata containing the sample names in the feature table. \
                    Deafault: 'Sample Name' ",
              metavar="Sample_Name"),
  
  make_option(c("-p", "--prevalence-column"), type="character", default=NULL, 
              help="Column in metadata to be use for prevalence based analysis. 
              Default: 'Sample_or_Control' ",
              metavar="prevalence_column"),
  
  make_option(c("-n", "--ntc-name"), type="character", default="Control_Sample", 
              help="Name of NTC in prevalence column. \
                    Deafault: 'Control_Sample' ",
              metavar="NTC_Name"),
  
  make_option(c("-F", "--frequency-column"), type="character", default=NULL, 
              help="Column in metadata to be use for frequency based analysis. \
                    Deafault: 'concentration' ",
              metavar="frequency_column"),
  
  make_option(c("-t", "--threshold"), type="numeric", default=0.5, 
              help="Dencontam's threshold (0-1) for both prevalence and frequency based analyses. \
                    Deafault: 0.1 ",
              metavar="threshold"),

  make_option(c("-M", "--classification-method"), type="character", default="", 
              help="Taxonomy or functional method used to generate the input 
              feature table. The supplied string will be added to output file names
              ['kaiju', 'kraken2', 'metaphlan', 'contig-taxonomy',
               'gene-taxonomy', 'gene-function']. Default: empty string.",
              metavar=""),
  
  make_option(c("-o", "--output-prefix"), type="character", default="", 
              help="Unique name to tag onto output files. Default: empty string.",
              metavar=""),

  make_option(c("-a", "--assay-suffix"), type="character", default="_GLMetagenomics", 
              help="Genelab assay suffix.", metavar="GLMetagenomics"),
  
  make_option(c("--version"), action = "store_true", type="logical", 
              default=FALSE,
              help="Print out version number and exit.", metavar = "boolean")
)




opt_parser <- OptionParser(
  option_list=option_list,
  usage = "Rscript %prog \\
                  --metadata-table 'mapping/metadata.csv' \\
                  --feature-table 'kaiju_species_table_GLlbnMetag' \\
                  --feature-column 'Species' \\ # ['Species', 'species', 'KO']
                  --prevalence-column 'NTC' \\
                  --frequency-column 'concentration' \\
                  --classification-method 'taxonomy' \\
                  --threshold 0.1 " ,
  description = paste("Author: Olabiyi Aderemi Obayomi",
                      "\nEmail: olabiyi.a.obayomi@nasa.gov",
                      "\n  A script to perform feature decontamination with decontam.",
                      "\nIt outputs decontam's primary output.",
                      "along with a decontaminated feature table",
                      sep="")
)


opt <- parse_args(opt_parser)


if (opt$version) {
  cat("run_decontam.R version: ", version, "\n")
  options_tmp <- options(show.error.messages=FALSE)
  on.exit(options(options_tmp))
  stop()
}



if(is.null(opt[["metadata-table"]])) {
  stop("Path to a metadata file must be set.")
}

if(is.null(opt[["feature-table"]])) {
  stop("Path to a feature table e.g. species/functions table file must be set.")
}

if(opt[["samples-column"]] == "Sample Name") {
  message("I will assume that the sample names are in a column named 'Sample Name' \n")
}


library(decontam)
library(phyloseq)
library(glue)
library(tidyverse)

# Feature table decontamination with decontam
run_decontam <- function(feature_table, metadata, contam_threshold=0.1,
                         prev_col=NULL, freq_col=NULL, ntc_name = "true") {
  
  # Retain metadata for only the samples present in the inpute feature table
  sub_metadata <- metadata[colnames(feature_table),]
  # Modify NTC concentration
  # Often times the user may set the NTC concentration to zero because they think nothing 
  # should be in the negative control but decontam fails if the value is set to zero.
  # To prevent decontam from failing, we replace zero with a very small concentration value
  # 0.0000001
  if (!is.null(freq_col)) {
    
    sub_metadata <- sub_metadata %>% 
      mutate(!!freq_col := map_dbl(!!sym(freq_col), .f= function(conc) { 
        if(conc == 0) return(0.0000001) else return(conc) 
      } 
      )
      )
    sub_metadata[, freq_col] <- as.numeric(sub_metadata[,freq_col])
    
  }
  
  # Create phyloseq object
  ps <- phyloseq(otu_table(feature_table, taxa_are_rows = TRUE), sample_data(sub_metadata))
  
  # In our phyloseq object, `prev_col` is the sample variable that holds the negative 
  # control information. We'll summarize the data as a logical variable, with TRUE for control 
  # samples, as that is the form required by isContaminant.
  sd <- as.data.frame(sample_data(ps)) # Extract sample metadata
  sd[,"is.neg"] <- 0 # Initialize
  sd[,"is.neg"] <- sample_data(ps)[[prev_col]] == ntc_name # Assign boolean value
  sample_data(ps) <- sd
  
  # Run Decontam 
  if (!is.null(freq_col) && !is.null(prev_col)) {   
    
    # Run decontam in both prevalence and frequency modes
    contamdf <- isContaminant(ps, neg="is.neg", conc=freq_col, threshold=contam_threshold) 
    
  } else if(!is.null(freq_col)) {
    
    # Run decontam in frequency mode
    contamdf <- isContaminant(ps, conc=freq_col, threshold=contam_threshold) 
    
  } else if(!is.null(prev_col)){
    
    # Run decontam in prevalence mode
    contamdf <- isContaminant(ps, neg="is.neg", threshold=contam_threshold)
    
  } else {
    
    cat("Both freq_col and prev_col cannot be set to NULL.\n")
    cat("Please supply either one or both column names in your metadata")
    cat("for frequency and prevalence based analysis, respectively\n")
    stop()
    
  }
  
  return(contamdf)
}


# Set input variables
feature_table_file <- opt[["feature-table"]] # 'kaiju_species_table_GLlbnMetag.csv'
metdata_file <- opt[["metadata-table"]] # "metadata.csv"
samples_column <-  opt[["samples-column"]] # 'Sample_ID'
freq_col <- opt[["frequency-column"]] # "input_conc_ng"
prev_col <- opt[["prevalence-column"]] # "NTC"
threshold <- opt[["threshold"]] # 0.5
ntc_name <- opt[["ntc_name"]] # "true"
# "kaiju", "kraken2", "metaphlan", "contig-taxonomy", "gene-taxonomy", "gene-function"
method <- opt[["classification-method"]] # 'kaiju'
feature_column <- opt[["feature-column"]] # 'Species'
prefix <- opt[["output-prefix"]]
suffix <- opt[["assay-suffix"]] # GLlbnMetag



# Prepare feature table
feature_table <- read_delim(feature_table_file) %>%  as.data.frame()
rownames(feature_table) <- feature_table[[1]]
feature_table <- feature_table[,-1]  %>% as.matrix()


# Prepare metadata
metadata <- read_delim(metdata_file) %>% as.data.frame()
row.names(metadata) <- metadata[,samples_column]

# Subset metadata and feature table  to contain the samples
samples <- intersect(colnames(feature_table), rownames(metadata))
metadata <- metadata[samples,]
feature_table <- feature_table[,samples]

# Run decontam
contamdf <- run_decontam(feature_table, metadata, threshold, prev_col, freq_col, ntc_name) 

contamdf <- as.data.frame(contamdf) %>% rownames_to_column(feature_column)

type <- "species"
if(method == "gene-function")  { type <- "KO"}

# Write decontaminated feature table and decontam's primary results
outfile <- glue("{prefix}{method}_decontam_{type}_results{suffix}.tsv")
write_tsv(x = contamdf, file = outfile)



# Get the list of contaminants identified by decontam
contaminants <- contamdf %>%
                   filter(contaminant == TRUE) %>%
                   pull(!!sym(feature_column))

# Drop contaminants(s) if detected
if(length(contaminants) > 0){
  
# Drop contaminant features identified by decontam
decontaminated_table <- feature_table %>% 
  as.data.frame()  %>% 
  rownames_to_column(feature_column) %>% 
  filter(str_detect(!!sym(feature_column), 
                    pattern = str_c(contaminants,
                                    collapse = "|"),
                    negate = TRUE))

rownames(decontaminated_table) <- decontaminated_table[[feature_column]]
decontaminated_table <- decontaminated_table[,-1] %>% as.matrix

outfile <- glue("{prefix}{method}_decontam_{type}_table{suffix}.tsv")
write_tsv(x = decontaminated_table, file = outfile)

}else{
  
  message("No contaminant was detected by Decontam")
}
