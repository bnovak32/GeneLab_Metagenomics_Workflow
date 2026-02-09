#!/usr/bin/env Rscript


###############################################################################
# AUTHOR : OLABIYI ADEREMI OBAYOMI
# DESCRIPTION: A script to create static and interactive barplots.
# E-mail: obadbotanist@yahoo.com
# Created: January 2025
# example: Rscript make_barplot.R \
#                  --metadata-table 'metadata.txt' \
#                  --feature-table 'kaiju_species_table_GLlbnMetag.tsv \
#                  --group-column 'Description' \
#                  --feature-column 'Species' \
#                  --samples-column 'Sample Name'  \
#                  --output-prefix  'filtered-kaiju_species'

###############################################################################

library(optparse)



######## -------- Get input variables from the command line ----##############

version <- 1.0 

# Input options 
option_list <- list(
  
  make_option(c("-m", "--metadata-table"), type="character", default=NULL, 
              help="path to a tab or comma separated samples metadata file with the 
               --prevalence-column and/or --frequency-column to be analyzed.",
              metavar="path"),
  
  make_option(c("-i", "--feature-table"), type="character", default=NULL, 
              help="path to a tab or comma separated samples feature table 
              i.e. species/functions table with species/functions as the first column 
              and samples as other columns.",
              metavar="path"),
  
  make_option(c("-f", "--feature-column"), type="character", default=NULL, 
              help="Feature column name in feature table ['Species', 'species', 'KO_ID'].
              Default: empty string",
              metavar="Feature_Column"),
  
  make_option(c("-s", "--samples-column"), type="character", default="Sample Name", 
              help="Column in metadata containing the sample names in the feature table. \
                    Deafault: 'Sample Name' ",
              metavar="Sample_Name"),
  
  make_option(c("-g", "--group-column"), type="character", default=NULL, 
              help="Column in metadata to be use to facet/group plot. 
              Default: NULL ",
              metavar="group_column"),
  
  
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
                  --metadata-table 'mapping/GLDS-487_amplicon_v1_runsheet.csv' \\
                  --feature-table '' \\
                  --feature-column 'Species' \\ # ['Species', 'species', 'KO']
                  --group-column 'Description' \\
                  --frequency-column 'concentration' \\
                  --method 'taxonomy' \\
                  --threshold 0.1 " ,
  description = paste("Author: Olabiyi Aderemi Obayomi",
                      "\nEmail: olabiyi.a.obayomi@nasa.gov",
                      "\n  A script to make static and interactive stacked barplots.",
                      "\nIt outputs a static and inetractive barplot.",
                      sep="")
)


opt <- parse_args(opt_parser)


if (opt$version) {
  cat("make_barplot.R version: ", version, "\n")
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


library(glue)
library(plotly)
library(htmlwidgets)
library(tidyverse)



# Convert species count matrix to relative abundance matrix
count_to_rel_abundance <- function(species_table) {
  
  abund_table <- species_table %>% 
    as.data.frame %>% 
    mutate( across(everything(), function(x) (x/sum(x, na.rm = TRUE))*100 ) )  %>% # calculate species relative abundance per sample
    select(
      where( ~all(!is.na(.)) )
    )  %>% # drop columns where none of the reads were classified or were non-microbial (NA)
    rownames_to_column("Species") 
  
  # Set rownames as species name and drop species column  
  rownames(abund_table) <- abund_table$Species
  abund_table <- abund_table[, -match(x = "Species", colnames(abund_table))] %>% t
  
  return(abund_table)
}

# Make bar plot
make_plot <- function(abund_table, metadata, colors2use, publication_format,
                      samples_column="Sample_ID", prefix_to_remove="barcode"){
  
  abund_table_wide <- abund_table %>% 
    as.data.frame() %>% 
    rownames_to_column(samples_column) %>% 
    inner_join(metadata) %>% 
    select(!!!colnames(metadata), everything()) %>% 
    mutate(!!samples_column := !!sym(samples_column) %>% str_remove(prefix_to_remove))
  
  
  abund_table_long <- abund_table_wide  %>%
    pivot_longer(-colnames(metadata), 
                 names_to = "Species",
                 values_to = "relative_abundance")
  
  p <- ggplot(abund_table_long, mapping = aes(x=!!sym(samples_column), 
                                              y=relative_abundance, fill=Species)) +
    geom_col() +
    scale_fill_manual(values = colors2use) + 
    labs(x=NULL, y="Relative Abundance (%)") + 
    publication_format
  
  return(p)
}


publication_format <- theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.ticks.length=unit(-0.15, "cm"),
        axis.text.x=element_text(margin=ggplot2::margin(t=0.5,r=0,b=0,l=0,unit ="cm")),
        axis.text.y=element_text(margin=ggplot2::margin(t=0,r=0.5,b=0,l=0,unit ="cm")), 
        axis.title = element_text(size = 18,face ='bold.italic', color = 'black'), 
        axis.text = element_text(size = 16,face ='bold', color = 'black'),
        legend.position = 'right', legend.title = element_text(size = 15,face ='bold', color = 'black'),
        legend.text = element_text(size = 14,face ='bold', color = 'black'),
        strip.text =  element_text(size = 14,face ='bold', color = 'black'))


custom_palette <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F", "#FF7F00",
                    "#CAB2D6","#6A3D9A","#FF00FFFF","#B15928","#000000","#FFC0CBFF","#8B864EFF","#F0027F",
                    "#666666","#1B9E77", "#E6AB02","#A6761D","#FFFF00FF","#FFFF99","#00FFFFFF",
                    "#B2182B","#FDDBC7","#D1E5F0","#CC0033","#FF00CC","#330033",
                    "#999933","#FF9933","#FFFAFAFF",colors()) 
# remove white colors
colors2use <- custom_palette[-c(21:23,
                                    grep(pattern = "white|snow|azure|gray|#FFFAFAFF|aliceblue",
                                         x = custom_palette, ignore.case = TRUE))]


feature_table_file <- opt[["feature-table"]] # 'kaiju_species_table_GLlbnMetag.csv'
metdata_file <- opt[["metadata-table"]] # "metadata.csv"
samples_column <- opt[["samples-column"]] # 'Sample_ID'
feature_column <- opt[["feature-column"]] # 'Species'
prefix <-  opt[["output-prefix"]] # "filtered-kaiju_species"
suffix <- opt[["assay-suffix"]]
facet_by <- reformulate(opt[["group-column"]]) # 'Description'


# Prepare feature table
feature_table <- read_delim(feature_table_file) %>%  as.data.frame()
rownames(feature_table) <- feature_table[[1]]
feature_table <- feature_table[,-1]

number_of_species <- nrow(feature_table)

if( number_of_species > length(colors2use)){

    N <- number_of_species / length(colors2use)
    colors2use  <- rep(colors2use, times=N*2)
}

# Prepare metadata
metadata <- read_delim(metdata_file) %>% as.data.frame()
row.names(metadata) <- metadata[,samples_column]

abund_table <- count_to_rel_abundance(feature_table)
p <- make_plot(abund_table , metadata, colors2use, publication_format, samples_column) +
     facet_wrap(facet_by, nrow=1, scales = "free_x") + 
     theme(axis.text.x = element_text(angle = 90))

static_plot <- p
number_of_species <- p$data$Species %>% unique() %>% length()

# Don't save legend if the number of species to plot is greater than 30
if(number_of_species > 30) {
  
  static_plot <- static_plot + theme(legend.position = "none")
  
}

width <- 2 * nrow(metadata) # 3.6 * number_of_samples
if(width < 14) { width = 14 } # set minimum width to 14 inches
if(width > 50) { width = 50 } # Cap plot with at 50 inches
# Save Static
ggsave(filename = glue("{prefix}_barplot{suffix}.png"), 
       plot = static_plot,
       device = 'png', width =width,
       height = 10, units = 'in', dpi = 300 , limitsize = FALSE)


# Save interactive
htmlwidgets::saveWidget(ggplotly(p), glue("{prefix}_barplot{suffix}.html"), selfcontained = TRUE)


