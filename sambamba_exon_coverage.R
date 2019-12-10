library(tidyverse)
library(plotly)
library(htmlwidgets)

# Purpose: This script takes the RAW output from sambamba and produces summary tables and plots highlighting the uniformity of coverage

# Usage: Rscript chanjo_exon_coverage.R --args "/path_to_folder/exon_coverage"

# Functions:

generate_coverage_plot <- function(df, panel) { 
  # Remove rows with NAs caused by regions not included between panels  
  df <- df[complete.cases(df), ]
  # Reorder the factors in region by median (Ensures the boxplots are plotted in order form lowest to highest)
  df$region <- fct_reorder(df$region, df$scaled_meanCoverage,.fun = median, .desc = FALSE)
  # Create a color palette to highlight samples by gene/transcript
  col=rainbow(length(levels(factor(df$gene))))[factor(df$gene)]
  # Plot coverage data (A series of boxplots showing coverage for each region, ordered by median)
  p <- df %>%
    ggplot( aes(x=region, y=scaled_meanCoverage)
    ) +
    geom_boxplot(outlier.size = 0.5, aes(fill=gene)) +
    geom_jitter(color="grey", width = 0.01, size=1, alpha=0, shape = 1) +
    theme(
      plot.title = element_text(size=11),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6)
    ) +
    ggtitle(paste0("Run ", run_name,",  ", panel ," (", num_target_regions," target regions), Coverage over ", num_samples," samples")) +
    xlab("Target Region") +
    ylab("Scaled average coverage")
  return(p)
}

generate_simple_coverage_plot <- function(df, panel) { 
  # Remove rows with NAs caused by regions not included between panels  
  df <- df[complete.cases(df), ]
  # Group the tibble data structure by 'region'
  region_mean <- df %>% 
    group_by(region) %>% 
    summarise(gene = unique(gene),
              transcript = unique(transcript),
              genomicCoordinates = unique(genomicCoordinates),
              region_meanCoverage = mean(scaled_meanCoverage))
  # Order region factors by region_meanCoverage to produce plot in correct order
  region_mean$region <- as.factor(region_mean$region)
  region_mean$region <- fct_reorder(region_mean$region, region_mean$region_meanCoverage,.fun = median, .desc = FALSE)
  # Plot region 
  region_mean %>%
    ggplot(aes(x=region, y=region_meanCoverage)) +
    geom_point(col = "red", shape = 20, size = 0.1) + 
    theme(
      legend.position = "none",
      plot.title = element_text(size=11),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 1)
    ) +
    ggtitle(paste0("Run ", run_name,",  ", panel ," (", num_target_regions," target regions), Coverage over ", num_samples," samples")) +
    xlab("Target Region") +
    ylab("Scaled average coverage")
}

# Uses scale() function on data - scaling is done by dividing the columns for each sample by their root mean square.
# This allows easier comparison between samples.
scale_this <- function(x) as.vector(scale(x, scale=TRUE, center = FALSE))

# Main Script:

# Get directory location from commandline - directory should contain the Raw exon level coverage files
data_directory <- commandArgs(trailingOnly = TRUE)

# Get all files with the suffix "*.bed" from data directory
sambamba_files <- list.files(path = data_directory, pattern = "*.refined.sambamba_output.bed", full.names = TRUE)

# Import coverage data and add relevant sample ID to each imported row 
tbl <- sapply(sambamba_files , read_tsv, col_types = "ciicicccinnc", simplify=FALSE) %>% 
  bind_rows(.id = "sample_id")

# Simplify & cleanup sample names
tbl$sample_id <- gsub(basename(tbl$sample_id), pattern = ".refined.sambamba_output.bed", replacement = "")
# Rename 2nd column to remove proceding '#'
colnames(tbl)[2] <- "chrom"
# Replace F1:F6 labels with meaningful names
colnames(tbl)[5:9] <- c("genomicCoordinates", "score", "strand", "gene_transcript", "accessionNum")

# Add new column 'region' so that each target region is represented by unique ID
tbl$region <- paste(tbl$chrom,
                    tbl$chromStart,
                    tbl$chromEnd,
                    tbl$gene_transcript,
                    sep=";")

# Group the tibble data structure by 'region'
tbl <- tbl %>% 
  group_by(sample_id) %>%
  mutate(scaled_meanCoverage = scale_this(meanCoverage))

# Identify Run ID from sample name and add as additional column
tbl$run_name <- stringr::str_split(string = tbl$sample_id, pattern = "_", simplify = TRUE)[,1] 
# Extract gene and transcript names into separate columns:
tbl$gene <- stringr::str_split(string = tbl$gene_transcript, pattern = ";", simplify = TRUE)[,1]
tbl$transcript <- stringr::str_split(string = tbl$gene_transcript, pattern = ";", simplify = TRUE)[,2]
# Any SNPs referenced by their RS accession number will not have a transcript - label as 'dbSNP'
tbl$gene[tbl$transcript==""] <- "dbSNP" 
# Identify Pan number from sample name and add as additional column
tbl$pan_number <- stringr::str_split(string = tbl$sample_id, pattern = "_", simplify = TRUE)[,7]
# Produce separate output for each panel

# Extract meta data from sample name
for(run_name in unique(tbl$run_name)){
for(panel in unique(tbl$pan_number)){

  df <- tbl[tbl$pan_number==panel,]
  
  # Update number of samples to be plotted 
  num_samples <- length(unique(df$sample_id))
  # Update number of target regions for this panel
  num_target_regions <- length(unique(df$region))  
  # Generate static plot of data for each
  static_plot <- generate_coverage_plot(df, panel)
  
  # Add interactivity to plot:
  interactive_plot <- ggplotly(static_plot)
  
  # Create coverage plot of means for PDF
  
  simplified_plot <- generate_simple_coverage_plot(df, panel)
  
  # Save interactive plot as a single html file:
  filename <- paste0(run_name, "_", panel, "_coverage.html")
  saveWidget(ggplotly(interactive_plot), file = filename)
  # Save simplified plot to pdf:
  filename <- paste0(run_name, "_", panel, "_coverage.pdf")
  ggsave(filename = filename, 
         simplified_plot,
         device = "pdf",
         width = 297,
         height = 200,
         units = "mm")
  # Save table
  filename <- paste0(run_name, "_", panel, "_coverage.csv")
  summary_df <- df %>% 
    group_by(region) %>% 
  # Summarise data by region
      summarise(gene = unique(gene),
              run_name = unique(run_name),
              panel = unique(panel),
              transcript = unique(transcript),
              genomicCoordinates = unique(genomicCoordinates),
              accessionNum = unique(accessionNum),
              region_meanCoverage = mean(scaled_meanCoverage)) %>%
    arrange(region_meanCoverage)
  write_delim(summary_df, filename, delim = "\t")
}
}
