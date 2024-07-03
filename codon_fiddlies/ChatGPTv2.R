# Load necessary libraries
if (!requireNamespace("seqinr", quietly = TRUE)) {
  install.packages("seqinr")
}

library(seqinr)

# Function to read FASTA alignment
read_fasta_alignment <- function(file) {
  sequences <- seqinr::read.fasta(file, as.string = TRUE, seqtype = "DNA")
  return(sequences)
}

# Function to write FASTA alignment
write_fasta_alignment <- function(sequences, file) {
  seqinr::write.fasta(sequences = sequences, names = names(sequences), file.out = file)
}

# Function to translate DNA sequence to protein
translate_sequence <- function(dna_seq) {
  protein_seq <- seqinr::translate(seqinr::s2c(dna_seq))
  return(protein_seq)
}

# Function to identify conserved peptides
identify_conserved_peptides <- function(protein1, protein2) {
  conserved <- protein1 == protein2
  return(conserved)
}

# Function to modify codons based on conserved peptides
modify_codons <- function(seq1, seq2, conserved_peptides) {
  codons1 <- substring(seq1, seq(1, nchar(seq1), 3), seq(3, nchar(seq1), 3))
  codons2 <- substring(seq2, seq(1, nchar(seq2), 3), seq(3, nchar(seq2), 3))
  
  modified_codons <- codons2
  for (i in seq_along(codons2)) {
    if (conserved_peptides[i]) {
      similar_codon <- codons1[i]
      modified_codons[i] <- similar_codon
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
  
  seq1 <- sequences[[1]]
  seq2 <- sequences[[2]]
  
  protein1 <- translate_sequence(seq1)
  protein2 <- translate_sequence(seq2)
  
  conserved_peptides <- identify_conserved_peptides(protein1, protein2)
  
  modified_seq2 <- modify_codons(seq1, seq2, conserved_peptides)
  
  sequences[[2]] <- modified_seq2
  
  write_fasta_alignment(sequences, output_fasta)
  cat("Modified alignment written to", output_fasta, "\n")
}

# Specify input and output file paths
input_fasta <- "Translation alignment of NG_MHCA_harmonized_gal.fasta"
output_fasta <- "output2.fasta"


# Run the main function
main(input_fasta, output_fasta)
