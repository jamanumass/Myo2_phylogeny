# Load necessary libraries
if (!requireNamespace("seqinr", quietly = TRUE)) {
  install.packages("seqinr")
}

library(seqinr)

# Function to read FASTA alignment
read_fasta_alignment <- function(file) {
  sequences <- read.fasta(file, as.string = TRUE)
  return(sequences)
}

# Function to write FASTA alignment
write_fasta_alignment <- function(sequences, file) {
  write.fasta(sequences = sequences, names = names(sequences), file.out = file)
}

# Function to find the most similar synonymous codon
find_best_synonymous_codon <- function(guide_codon, candidate_codons) {
  best_codon <- candidate_codons[1]
  max_similarity <- -1
  
  for (candidate in candidate_codons) {
    similarity <- sum(unlist(strsplit(guide_codon, "")) == unlist(strsplit(candidate, "")))
    if (similarity > max_similarity) {
      max_similarity <- similarity
      best_codon <- candidate
    }
  }
  
  return(best_codon)
}

# Function to modify codons
modify_codons <- function(guide_seq, target_seq) {
  codons_guide <- substring(guide_seq, seq(1, nchar(guide_seq), 3), seq(3, nchar(guide_seq), 3))
  codons_target <- substring(target_seq, seq(1, nchar(target_seq), 3), seq(3, nchar(target_seq), 3))
  
  modified_codons <- codons_target
  for (i in seq_along(codons_target)) {
    # Get synonymous codons for the target codon
    synonymous_codons <- tryCatch(as.character(unlist(syncodons(codons_target[i]))), error = function(e) NULL)
    
    if (!is.null(synonymous_codons) && length(synonymous_codons) > 0) {
      # Find the best synonymous codon that matches the guide codon
      best_codon <- find_best_synonymous_codon(codons_guide[i], synonymous_codons)
      if (best_codon != codons_target[i]) {
        modified_codons[i] <- best_codon
      }
    }
  }
  
  modified_seq <- paste(modified_codons, collapse = "")
  return(modified_seq)
}

# Main function
main <- function(input_fasta, output_fasta) {
  sequences <- read_fasta_alignment(input_fasta)
  if (length(sequences) < 2) {
    stop("The input FASTA file must contain at least two sequences.")
  }
  
  guide_seq <- sequences[[1]]
  target_seq <- sequences[[2]]
  
  modified_target_seq <- modify_codons(guide_seq, target_seq)
  
  sequences[[2]] <- modified_target_seq
  
  write_fasta_alignment(sequences, output_fasta)
  cat("Modified alignment written to", output_fasta, "\n")
}

input_fasta <- "Translation alignment of NgA_optimised to Dd_native.fasta"
output_fasta <- "output1.fasta"

# Run the main function
main(input_fasta, output_fasta)
