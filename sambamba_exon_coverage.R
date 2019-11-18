library(readr)
library(tidyverse)
library(plotly)
library(magrittr)
library(htmlwidgets)

# Purpose: This script....

# Usage:
# Rscript ./chanjo_exon_coverage.R --args "/home/graeme/Desktop/NGS300_coverage"

# Functions:

generate_coverage_plot <- function(df, panel) { 
  # Plot coverage data (A series of boxplots showing coverage for each region, ordered by median)
  p <- df %>%
    ggplot( aes(x=region, y=meanCoverage)
    ) +
    # geom_point(size = 0.25, shape = 1) +
    geom_boxplot(outlier.size = 0.5, aes(fill=gene)) +
    geom_jitter(color="grey", width = 0.01, size=1, alpha=0, shape = 1) +
    theme(
      #legend.position="none",
      plot.title = element_text(size=11),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6)
      #axis.text.x = NULL
    ) +
    ggtitle(paste("Run ", run_name,",  ", panel ,", Coverage over ", num_samples,"samples")) +
    xlab("Target Region") +
    ylab("Scaled average coverage")
  return(p)
}

# Script:

# Get directory location from commandline - directory should contain the Raw exon level coverage files
data_directory <- commandArgs(trailingOnly = TRUE)

# Import all files with the suffix "*.bed" from data directory
sambamba_files <- list.files(path = data_directory, pattern = "*.refined.sambamba_output.bed", full.names = TRUE)

# Import coverage data
tbl <- sapply(sambamba_files , read_tsv, col_types = "ciicicccinnc", simplify=FALSE) %>% 
  bind_rows(.id = "sample_id")

# Simplify & cleanup sample names
tbl$sample_id <- gsub(basename(tbl$sample_id), pattern = ".refined.sambamba_output.bed", replacement = "")
# Rename 2nd column to remove proceding '#'
colnames(tbl)[2] <- "chrom"

#  Convert dataframe to wide format:
tbl$sample_id <- factor(tbl$sample_id)
tbl <- tibble(sample_id = tbl$sample_id, 
                      region = paste(tbl$chrom,
                                     tbl$chromStart,
                                     tbl$chromEnd,
                                     tbl$F6,
                                     sep=";"),
                      meanCoverage = tbl$meanCoverage)

data_wide <- spread(tbl, sample_id, meanCoverage)

# Use scale to centre the data - scaling is done by dividing the columns of x by their root mean square.
data_wide[,-1] <- scale(data_wide[,-1], scale=TRUE, center = FALSE)

# Cast back to long format for visualisation
data_long <- gather(data_wide, sample, meanCoverage, 2:33)

# Visualise per exon metrics

num_samples <- length(unique(data_long$sample))
run_name <- strsplit(data_long$sample, "_")[[1]][1]

# Split column
data_long$region2 <- data_long$region
data_long <- separate(data = data_long,
                      col = region2,
                      into = c("chrom", "chromStart", "chromEnd", "gene", "transcript"),
                      sep = ";")
data_long$gene[is.na(data_long$transcript)] <- "SNP" # Any SNPs referenced by their RS accession number will not have a transcript
data_long$gene_transcript <- paste0(data_long$gene, ":", data_long$transcript)

# Reorder the factors in region by median (Ensures the boxplots are plotted in order form lowest to highest)
data_long$region <- fct_reorder(data_long$region, data_long$meanCoverage,.fun = median, .desc = FALSE)

# Create a color palette to highlight samples by gene/transcript
col=rainbow(length(levels(factor(data_long$gene))))[factor(data_long$gene)]

# Identify Pan number from sample name and add as additional column
data_long$pan_number <- stringr::str_split(string = data_long$sample, pattern = "_", simplify = TRUE)[,7]

for(panel in unique(data_long$pan_number)){
  df <- data_long[data_long$pan_number==panel,]
  
  # Generate static plot of data for each
  static_plot <- generate_coverage_plot(df, panel)
  
  # Add intactivity to plot:
  interactive_plot <- ggplotly(static_plot)
  
  # Save interactive plot as a single html file:
  filename <- paste0(run_name, "_", panel, "_coverage.html")
  saveWidget(ggplotly(interactive_plot), file = filename)
}

