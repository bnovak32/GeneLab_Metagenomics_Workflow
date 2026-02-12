#!/usr/bin/env Rscript


###############################################################################
# AUTHOR : OLABIYI ADEREMI OBAYOMI
# DESCRIPTION: A script to filter a feature table based on user defined criteria.
# E-mail: obadbotanist@yahoo.com
# Created: January 2026

# max_value - relative abundance filtering approach to group rare species
# example: Rscript filter_feature_table.R \\
#                  --feature-table 'kaiju_species_table_GLlbnMetag.tsv' \\
#                  --mode  'max_value' \\
#                  --threshold 0.05 \\
#                  --output-file  'kaiju_filtered_species_table_GLlbnMetag.tsv'

# across_samples - relative read-based filtering accross samples
# example: Rscript filter_feature_table.R \
#                  --feature-table 'kraken_species_table_GLlbnMetag.tsv' \
#                  --mode  'across_samples' \
#                  --threshold 0.5 \\
#                  --output-file  'kraken_filtered_species_table_GLlbnMetag.tsv'

# values_sum - assembly based CPM filtering
# example: Rscript filter_feature_table.R \\
#                  --feature-table 'Gene_function_table_GLMetagenomics.tsv' \\
#                  --mode  'values_sum' \\
#                 --threshold 1000 \\ # filter out features (KO_ID) with with less than 1000CPM across samples  
#                 --output-file  'Gene_function_table_filtered_GLlbnMetag.tsv'
###############################################################################

library(optparse)



######## -------- Get input variables from the command line ----##############

version <- 1.0 



# Input options 
option_list <- list(
  
  make_option(c("-f", "--feature-table"), type="character", default=NULL, 
              help="path to a tab or comma separated samples feature table 
              i.e. species/functions table with species/functions as the first column 
              and samples as other columns.",
              metavar="path"),
  make_option(c("-o", "--output-file"), type="character", default='filtered_table.csv', 
              help="path to a tab separated samples feature table 
              i.e. species/functions table with species/functions as the first column 
              and samples as other columns.",
              metavar="path"),
  make_option(c("-m", "--mode"), type="character", default="", 
              help="Feature(s) filtering mode to use to use. 
              Options are ['max_value', 'across_samples', 'values_sum'].
              max_value - group rare species [standard rare filtering approach doesn't consider the abundance across all samples]
              across_samples - relative abundance of features per sample (read_based filtering considering the abundance of the feature/taxa across all samples)
              values_sum  - sum of feature values (CPM) [assembly based]
              Default: across_samples .",
              metavar=""),
  make_option(c("-d", "--features-to-drop"), type="character", default=NULL, 
              help="Comma separated list of feature names to drop. Default: NULL.",
              metavar=""),
  make_option(c("-t", "--threshold"), type="numeric", default=NULL, 
              help="A threshold for filtering out rare features (taxa and funtions). When \
              mode == 'max_value' it will be a  > 0 and < 100 \
              mode == 'accros_samples' it will be a  > 0 and < 100 \
              mode == 'values_sum' it will be a value greater than 0. default is 1000 for 1000 CPM",
              metavar="threshold"),
  make_option(c("-p", "--output-prefix"), type="character", default="", 
              help="Unique name to tag onto output files. Default: empty string.",
              metavar=""),
  
  make_option(c("-a", "--assay-suffix"), type="character", default="_GLMetagenomics", 
              help="Genelab assay suffix.", metavar="GLMetagenomics"),
  
  make_option(c("--version"), action = "store_true", type="logical", 
              default=FALSE,
              help="Print out version number and exit.", metavar = "boolean")
)


# ----------   Based of sum of feature values ------------#
# values_sum [assembly based data]
get_abundant_features <- function(mat, cpm_threshold=1000){
  
  # mat - matrix with features as rows and samples as columns
  # cpm_threshold - threshold to filter abundant features
  features <- rowSums(mat) %>% sort()
  
  abund_features <- features[features > cpm_threshold] %>% names
  
  abund_features.m <- mat[abund_features,]
  
  return(abund_features.m)
  
}



#-------------------  Based on relative abundance -----------------------#
# across_samples [read based data]
# Filter out rare and non_microbial assignments 
filter_rare <- function(species_table, non_microbial, threshold=1) {
  
  # species_table - species table with rows as species and columns as samples
  # non_microbial - a regular expression of names to be dropped from the species table before filtering by relative abundance
  # threshold - relative abundance filtering threshold in percentage where 1 stands for 1%
  
  # Drop species listed in 'non_microbial' regex
  clean_tab_count  <-  species_table %>% 
    as.data.frame %>% 
    rownames_to_column("Species") %>% 
    filter(str_detect(Species, non_microbial, negate = TRUE))  
  # Calculate species relative abundance
  clean_tab <- clean_tab_count %>% 
    mutate( across( where(is.numeric), function(x) (x/sum(x, na.rm = TRUE))*100 ) )
  # Set rownames as species name and drop species column  
  rownames(clean_tab) <- clean_tab$Species
  clean_tab  <- clean_tab[,-1] 
  
  
  # Get species with relative abundance less than `threshold` in all samples
  rare_species <- map(clean_tab, .f = function(col) rownames(clean_tab)[col < threshold])
  rare <- Reduce(intersect, rare_species)
  
  # Set rownames as species name and drop species column  
  rownames(clean_tab_count) <- clean_tab_count$Species
  clean_tab_count  <- clean_tab_count[,-1] 
  # Drop rare species
  abund_table <- clean_tab_count[!(rownames(clean_tab_count) %in% rare), ]
  
  return(abund_table)
}


# Function to group rare taxa or return a table with the rare taxa
group_low_abund_taxa <- function(abund_table, threshold=0.05,
                                 rare_taxa=FALSE) {
  # abund_table is a relative abundance matrix with taxa as columns and  samples as rows
  #rare_taxa is a boolean specifying if only rare taxa should be returned
  #If set to TRUE then a table with only the rare taxa will be returned 
  #intialize an empty vector that will contain the indices for the
  #low abundance columns/ taxa to group
  taxa_to_group <- c()
  #intialize the index variable of species with low abundance (taxa/columns)
  index <- 1
  
  #loop over every column or taxa check to see if the max abundance is less than the set threshold
  #if true save the index in the taxa_to_group vector variable
  for (column in ncol(abund_table)){
    if(max(abund_table[,column], na.rm = TRUE) < threshold ){
      #print(column)
      taxa_to_group[index] <- column
      index = index + 1
    }
  }
  
  if(is.null(taxa_to_group)){
    
    
    message(glue::glue("Rare taxa were not grouped. please provide a higher 
                       threshold than {threshold} for grouping rare taxa, 
                       only numbers are allowed."))
    return(abund_table)
  }
  
  
  
  if(rare_taxa){
    abund_table <- abund_table[,taxa_to_group,drop=FALSE]
  }else{
    #remove the low abundant taxa or columns
    abundant_taxa <-abund_table[,-(taxa_to_group), drop=FALSE]
    #get the rare taxa
    # rare_taxa <-abund_table[,taxa_to_group]
    rare_taxa <- subset(x = abund_table, select = taxa_to_group)
    #get the proportion of each sample that makes up the rare taxa
    rare <- rowSums(rare_taxa)
    #bind the abundant taxa to the rae taxa
    abund_table <- cbind(abundant_taxa,rare)
    #rename the columns i.e the taxa
    colnames(abund_table) <- c(colnames(abundant_taxa),"Rare")
  }
  
  return(abund_table)
}


library(tidyverse)
library(glue)

opt_parser <- OptionParser(
  option_list=option_list,
  usage = "Rscript %prog \\
                  --feature-table 'kaiju_species_table_GLlbnMetag.csv' \\
                  --mode  'across_value' \\
                  --threshold 0.5 \\
                  --output-file  'kaiju_filtered_species_table_GLlbnMetag.csv' " ,
  description = paste("Author: Olabiyi Aderemi Obayomi",
                      "\nEmail: olabiyi.a.obayomi@nasa.gov",
                      "\n  A script to filter a feature table based on user defined criteria..",
                      "\nIt outputs a filtered feature (species or functions) table with the first column as feature name and other columns as sample names",
                      sep="")
)



opt <- parse_args(opt_parser)



if (opt$version) {
  cat("filter_feature_table.R version: ", version, "\n")
  options_tmp <- options(show.error.messages=FALSE)
  on.exit(options(options_tmp))
  stop()
}



feature_table_file <- opt[['feature-table']]
threshold <- opt[['threshold']]


if(!is.null(opt[['features-to-drop']])){
  
non_microbial <-  opt[['features-to-drop']]
  
}else{
  
  non_microbial <- "UNCLASSIFIED|Unclassifed|unclassified|Homo sapien|cannot|uncultured|unidentified"
}

feature_table <- read_delim(feature_table_file) %>% as.data.frame()
feature_name <- colnames(feature_table)[1]
rownames(feature_table) <- feature_table[,1]
feature_table <- feature_table[, -1]

if( opt[['mode']] == "values_sum"){
  
  # Assembly CPM tables  
  table2write <- get_abundant_features(feature_table, cpm_threshold=threshold) %>%
                      as.data.frame() %>% 
                      rownames_to_column(feature_name)
  
  
}else if(opt[['mode']] == "across_samples"){

  # Read- based cont tables
  table2write <- filter_rare(feature_table, non_microbial, threshold=threshold) %>%
                     as.data.frame() %>% 
                     rownames_to_column(feature_name)
  
}else if(opt[['mode']] == "max_value"){
  # Relative abundance table for which rare taxa should be grouped
  
  # convert count table to a relative abundance matrix
  abund_table <- feature_table %>% rownames_to_column(feature_name) %>%  
                 mutate( across( where(is.numeric), function(x) (x/sum(x, na.rm = TRUE))*100 ) ) %>%
                    as.data.frame()
  
  rownames(abund_table) <- abund_table[,1]
  abund_table  <- abund_table[,-1] %>% t 
  
  table2write <-  group_low_abund_taxa(abund_table, threshold=threshold)  %>%
                        t %>% as.data.frame() %>%
                       rownames_to_column(feature_name)
  
}else{
  
  message("Please provide a valid option for filtering features. One of ['max_value', 'across_samples', 'values_sum']")
  stop()
}

#GLDS-XXX_metagenomics-lowbiomass-longread_kaiju_filtered_taxon_counts_GLlbnMetag.tsv

write_tsv(x = table2write, file = opt[['output-file']])
