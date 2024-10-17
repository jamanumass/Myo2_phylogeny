#!/bin/bash

# Set the location of your protein databases
protein_database_directory="/work/pi_lfritzlaylin_umass_edu/users/jaman/sequence_databases/myo2_seqs_collect_v1/full_genomes_peptide_databases"

# File containing the list of gene IDs (one per line, make sure a carriage rturn after last ID, making blank line at end)
gene_ids_file="test.txt"  # Replace with your file

# Output file to store the extracted sequences
output_fasta_file="test_grep_extracted_sequences.fasta"

# Create or clear the output file
> "$output_fasta_file"

# Loop through each protein database file in the directory
for protein_database_file_path in "${protein_database_directory}"/*_protein.faa; do
    # Set the current protein genome being searched
    protein_database_file=$(basename "$protein_database_file_path")

    echo "Searching in database: $protein_database_file"

    # Loop through each gene ID in the file
    while IFS= read -r gene_id; do
        echo "Searching for gene ID: $gene_id"
        # Use awk to extract the sequence for each gene ID
        awk -v id="$gene_id" '
        BEGIN { found=0 }
        $0 ~ ">" && $0 ~ id { found=1; print; next }
        found && /^>/ { found=0 }
        found { print }
        ' "$protein_database_file_path" >> "$output_fasta_file"
    done < "$gene_ids_file"
done

echo "Extraction complete. Sequences saved to $output_fasta_file."

#quickly check if you got them all, these numbers should be equal
cat $gene_ids_file | wc -l
grep -c ">" "$output_fasta_file"

