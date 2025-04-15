#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16         # Reduced CPU cores per task
#SBATCH --mem=60G                   # Reduced memory
#SBATCH -t 01:00:00
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err

module load uri/main
module load HMMER/3.3.2-iimpi-2021b
module load MAFFT/7.505-GCC-11.3.0-with-extensions

# Source the variable file to load all variables into the environment
if [ -f "$shared_variables_file" ]; then
    source $shared_variables_file
    IFS=' ' read -r -a database_dirs <<< "$database_dirs"
else
    echo "Error: Variable file not found."
    exit 1
fi

echo "Using the following settings from the main script:"
echo "Input FASTA: $INPUT_FASTA_FOR_HMMSEARCH"
if [ "$use_hmm_score_threshold" = true ]; then
        echo "Using HMM score threshold: $hmm_score_threshold"
    elif [ "$use_hmm_score_threshold" = false ]; then
       echo "Number of hits from each genome to collect: $num_top_hmm_hits_per_genome"
    else
        echo "Error: Invalid value for use_hmm_score_threshold. Must be true or false."
        exit 1
    fi
 

# 2. Prepare input for HMM
#cp "$INPUT_FASTA_FOR_HMMSEARCH" "${run_dir}/${INPUT_FASTA_FOR_HMMSEARCH}" #this just to provide a copy for future
#INPUT_FASTA_FOR_HMMSEARCH="${run_dir}/${run_name}_input.fasta" # Use the copied version in the local results directory
#above is probably a mistake to add, probably delete later

# 3. Check if the input is aligned and align if necessary
    aligned_input_fasta="${run_dir}/${run_name}_input_aligned.fasta" # Define the aligned file this step will end with
if awk 'NR==2 {if (gsub(/-/, "&") >= 5) exit 1; else exit 0}' "$INPUT_FASTA_FOR_HMMSEARCH"; then
    echo "Not an alignment. Running MAFFT..."
    mafft --maxiterate 1000 --localpair --reorder --thread $((NUM_CORES - 1)) $INPUT_FASTA_FOR_HMMSEARCH > $aligned_input_fasta
else
    echo "Input is already an alignment. Proceeding"
    cat $INPUT_FASTA_FOR_HMMSEARCH > $aligned_input_fasta
fi

# 4. Build the HMM from the input
hmmbuild_file="${hmm_search_results_dir}/${run_name}.hmm"
hmmbuild "${hmmbuild_file}" "${aligned_input_fasta}"

#set up files for collected IDs and sequences
hmm_collected_seqs_file="${hmm_search_results_dir}/${run_name}_collected_seqs.fa"
hmm_collected_IDs_file="${hmm_search_results_dir}/${run_name}_collected_IDs.txt"
touch "${hmm_collected_seqs_file}"

# 5a. Define an array of protein database filenames to exclude
# Add the exact filenames you want to exclude from processing
excluded_databases=(
"EP00057_Strongylocentrotus_purpuratus.fasta"
"EP00058_Branchiostoma_floridae.fasta"
"EP00059_Ciona_intestinalis.fasta"
"EP00067_Danio_rerio.fasta"
"EP00074_Homo_sapiens.fasta"
"EP00076_Gallus_gallus.fasta"
"EP00081_Caenorhabditis_elegans.fasta"
"EP00090_Calanus_glacialis.fasta"
"EP00096_Caligus_rogercresseyi.fasta"
"EP00099_Drosophila_melanogaster.fasta"
"EP00103_Capitella_teleta.fasta"
"EP00110_Nematostella_vectensis.fasta"
"EP00114_Trichoplax_sp_H2.fasta"
"EP00115_Mnemiopsis_leidyi.fasta"
"EP00118_Oscarella_pearsei.fasta"
"EP00119_Amphimedon_queenslandica.fasta"
"EP00179_Madagascaria_erythrocladioides.fasta"
   
        # Add more filenames as needed
)

# Function to check if a file is in the excluded list
is_excluded() {
    local filename="$1"
    for excluded in "${excluded_databases[@]}"; do
        if [[ "$filename" == "$excluded" ]]; then
            return 0  # True: excluded
        fi
    done
    return 1  # False: not excluded
}

# 5b. Loop through each database file and perform the search
echo "Searching each genome database using HMM..."

# Loop through each database directory and search
for dir in "${database_dirs[@]}"; do
    for file_path in "${dir}/"*.fasta "${dir}/"*_protein.faa; do
        # Skip if no matching files
        if [ ! -f "$file_path" ]; then
            continue
        fi
        
        # Extract the protein database file name
        protein_database_file=$(basename "$file_path")
        
        # Check if the current database is in the exclusion list
        if is_excluded "$protein_database_file"; then
            echo "Skipping excluded database: $protein_database_file"
            continue
        fi
        
        # Handle both file types (.fasta and _protein.faa)
        if [[ "$protein_database_file" == *_protein.faa ]]; then
            protein_database_name="${protein_database_file%_protein.faa}"
        else
            protein_database_name="${protein_database_file%.fasta}"
        fi
    
        # Make a directory for the current genome's search results
        current_species_search_results_dir="${hmm_search_results_dir}/search_results_${protein_database_name}"
        mkdir -p "${current_species_search_results_dir}"
        echo "Searching in protein database: ${protein_database_file}"
    
        # Perform the HMM search
        hmmsearch_hits_file="${current_species_search_results_dir}/${run_name}_in_${protein_database_name}_hmmsearch_out.txt"
        hmmsearch --tblout "${hmmsearch_hits_file}" --cpu "${NUM_CORES}" --noali "${hmmbuild_file}" "${file_path}" > /dev/null # Suppress extra readouts to .out
    
        # Filter results based on either the HMM score threshold or total number of hits
        thresholded_IDs_file="${current_species_search_results_dir}/${protein_database_name}_thresholded_IDs.txt"
        if [ "$use_hmm_score_threshold" = true ]; then
            # Use hmm_score_threshold to filter hits
            grep -v '^#' "${hmmsearch_hits_file}" | awk -v threshold="${hmm_score_threshold}" '$6 >= threshold {print $1}' > "${thresholded_IDs_file}"
        elif [ "$use_hmm_score_threshold" = false ]; then
            # Use num_top_hmm_hits_per_genome to get top N hits
            grep -v '^#' "${hmmsearch_hits_file}" | head -n "${num_top_hmm_hits_per_genome}" | awk '{print $1}' > "${thresholded_IDs_file}"
        else
            echo "Error: Invalid value for use_hmm_score_threshold. Must be true or false."
            exit 1
        fi
    
        # Add the filtered hits gene IDs to the run's collected IDs file
        cat "${thresholded_IDs_file}" >> "${hmm_collected_IDs_file}"
        
        # Extract the fasta sequences of the filtered hits
        thresholded_seqs_file="${current_species_search_results_dir}/${protein_database_name}_thresholded_seqs.fasta"
        touch "${thresholded_seqs_file}"
        
        # Collect the sequences from the thresholded_IDs_file
        while IFS= read -r gene_id; do
            # For UniProt-style IDs (sp|ID|...), extract the second field
            if [[ "$gene_id" =~ ^sp\| ]]; then
                id=$(echo "$gene_id" | cut -d'|' -f2)
            # For NCBI-style IDs, assume the full gene_id (XP_... or NP_...) is used
            else
                id="$gene_id"
            fi
            
            awk -v id="$id" '
            BEGIN { found=0 }
            # Match either UniProt or NCBI format based on the extracted id
            $0 ~ ">" && $0 ~ id { found=1; print; next }
            found && /^>/ { found=0; exit }  # Stop processing after finding the sequence
            found { print }
            ' "${file_path}" >> "${thresholded_seqs_file}"
        done < "${thresholded_IDs_file}"
        
        # Append these sequences to the overall collected sequences file
        cat "${thresholded_seqs_file}" >> "${hmm_collected_seqs_file}"
        echo "Collected $(wc -l < "${thresholded_IDs_file}") sequences from ${protein_database_file}"
    
    done
done

echo "Finished collecting results from the HMM search. Gene IDs found are in ${hmm_collected_IDs_file} and their sequences are in ${hmm_collected_seqs_file}"


# Append HMMER-specific variables to the job variables file
echo "hmm_collected_seqs_file=${hmm_collected_seqs_file}" >> ${shared_variables_file}
echo "hmm_collected_IDs_file=${hmm_collected_IDs_file}" >> ${shared_variables_file}
echo "hmmbuild_file=${hmmbuild_file}" >> ${shared_variables_file}
