# Load necessary libraries
library(rvest)
library(dplyr)
library(stringr)
library(rentrez)
library(ggplot2)
library(ape)

# Function to extract genomic reference sequence ID and coordinates
get_reference_assembly_info <- function(protein_id) {
  url <- paste0("https://www.ncbi.nlm.nih.gov/gene/?term=", protein_id)
  webpage <- read_html(url)
  page_text <- webpage %>%
    html_text(trim = TRUE)
  lines <- unlist(strsplit(page_text, "\n"))
  range_index <- grep("Range", lines, ignore.case = TRUE)
  
  if (length(range_index) > 0) {
    nearby_lines <- lines[max(1, range_index - 5):min(length(lines), range_index + 5)]
    reference_id <- str_extract(nearby_lines, "NW_\\d+\\.\\d+")
    reference_id <- reference_id[!is.na(reference_id)]  # Remove NA values
    
    if (length(reference_id) > 0) {
      reference_id <- reference_id[1]  # Use the first match if multiple
    } else {
      reference_id <- NA
    }
    
    range_line <- str_trim(lines[range_index + 2])
    return(list(reference_id = reference_id, coordinates = range_line))
  } else {
    return(NULL)  # Return NULL if no data is found
  }
}

# Function to fetch the GenBank content based on reference ID and genomic range
fetch_genbank_content <- function(reference_id, start, end) {
  genbank_content <- entrez_fetch(db = "nuccore", id = reference_id, rettype = "gb", seq_start = start, seq_stop = end)
  return(genbank_content)
}

# Function to extract exon coordinates from GenBank content, including open-ended exons
extract_exon_coordinates <- function(genbank_content) {
  cds_line <- str_extract(genbank_content, "join\\([^\\)]+\\)")
  
  if (is.na(cds_line)) {
    return(NULL)  # Return NULL if no CDS coordinates are found
  }
  
  exon_pairs <- unlist(str_extract_all(cds_line, "\\>?\\d+\\.\\.\\>?\\d+"))
  
  exon_df <- do.call(rbind, lapply(exon_pairs, function(pair) {
    coords <- unlist(strsplit(pair, "\\.\\."))
    start <- as.numeric(gsub(">", "", coords[1]))
    end <- as.numeric(gsub(">", "", coords[2]))
    return(data.frame(start = start, end = end))
  }))
  
  return(exon_df)
}

# Function to map exons to a list of gene IDs
map_exons_to_genes <- function(gene_ids) {
  # Limit to the first 3 gene IDs
  #gene_ids <- head(gene_ids, 3)
  
  exon_data <- list()  # Initialize an empty list to store exon data for each gene
  
  for (protein_id in gene_ids) {
    message("Processing protein ID: ", protein_id)
    
    reference_info <- get_reference_assembly_info(protein_id)
    
    if (!is.null(reference_info) && !is.na(reference_info$reference_id)) {
      coords <- strsplit(reference_info$coordinates, "\\.\\.")[[1]]
      start <- as.numeric(coords[1])
      end <- as.numeric(coords[2])
      
      genbank_content <- fetch_genbank_content(reference_info$reference_id, start, end)
      exon_coordinates <- extract_exon_coordinates(genbank_content)
      
      exon_data[[protein_id]] <- exon_coordinates
      message("Exon data found for protein ID: ", protein_id)
    } else {
      exon_data[[protein_id]] <- NULL  # Leave it blank if no data is found
      message("No exon data found for protein ID: ", protein_id)
    }
  }
  
  return(exon_data)  # Return the list of exon data
}