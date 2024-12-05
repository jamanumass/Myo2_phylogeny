#!/bin/bash
#SBATCH --partition=cpu        # Partition to use
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks=1                  # Number of tasks (processes)
#SBATCH --cpus-per-task=32           # Number of CPU cores per task
#SBATCH --mem=120G                   # Memory per node
SBATCH -t 04:00:00
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err

module load uri/main
module load HMMER/3.3.2-iimpi-2021b
module load mafft/7.481



# 1 Source the variable file to load all variables into the environment
if [ -f "$output_var_file" ]; then
    source $output_var_file
else
    echo "Error: Variable file not found."
    exit 1
fi

echo "Using the following settings from the main script:"
echo "Input FASTA: $INPUT_FASTA"
echo "Number of results to collect: $number_results_to_collect"




# 2 prepare input for HMM
cp "$INPUT_FASTA" "${run_dir}/${run_name}_input.fasta" # Copy the input FASTA to the run directory

# Update INPUT_FASTA to point to the copied version in the run directory
INPUT_FASTA="${run_dir}/${run_name}_input.fasta"

# Check if input is an alignment, and align if it is not
if awk 'NR==2 {if (gsub(/-/, "&") >= 5) exit 1; else exit 0}' "$INPUT_FASTA"; then
    echo "Not an alignment. Running MAFFT..."
    program_call="mafft --maxiterate 1000 --localpair --reorder --thread ${NUM_CORES} $INPUT_FASTA > ${run_dir}/${run_name}_input_aligned.fasta"
        echo "MAFFT to align input seqs:" >> ${timing_log}
        { time $program_call ; } 2>> ${timing_log}
    INPUT_FASTA="${run_dir}/${run_name}_input_aligned.fasta"  # Update INPUT_FASTA to point to the aligned file
else
    echo "File is already an alignment. No action taken."
fi

# 3 Construct the HMM file from the input
#set up files for collected IDs and sequences
collected_seqs_file="${hmm_search_results_dir}/${run_name}_collected_seqs.fa"
collected_IDs_file="${hmm_search_results_dir}/${run_name}_collected_IDs.txt"
touch "${collected_seqs_file}"
touch "${collected_IDs_file}"

#Create the HMM for the input
hmmbuild_file="${hmm_search_results_dir}/${run_name}.hmm"
hmmbuild "${hmmbuild_file}" "${INPUT_FASTA}"



# 4 Search the databases for best hits from the HMM file generated in previous step

#loop here through each protein database file and perform the search and sequence collection
echo "Find best ${number_results_to_collect} hits from each genome:" >> ${timing_log}
{ time for file_path in "${protein_database_directory}"/*_protein.faa; do #Set the current protein genome being searched
    protein_database_file=$(basename "$file_path")


    #Set database directories based on the current working protein genome
    protein_database_name="${protein_database_file%_protein.faa}"
    protein_database_file_path="${protein_database_directory}/${protein_database_file}"

    #Make a new sub-directory for the search results of the current working protein genome
    current_species_search_results_dir_name="search_results_${protein_database_name}"
    current_species_search_results_dir="${hmm_search_results_dir}/${current_species_search_results_dir_name}" 
    mkdir -p ${current_species_search_results_dir}
    
    #Dictate where the search result files will go and be named
    hmmsearch_hits_file="${current_species_search_results_dir}/${run_name}_in_${protein_database_name}_hmmsearch_out.txt"
    
    #Perform the HMM search
    hmmsearch --tblout "${hmmsearch_hits_file}" --cpu 32 --noali "${hmmbuild_file}" "${protein_database_file_path}"
    
    # Extract the top N hits from the output table and save to top hits file
    top_hits_file="${current_species_search_results_dir}/${protein_database_name}_top_hits_table.txt"
    #create the top hits file by extracting the info from the full HMMsearch hits file. Maintains same columns
    #In future look into HMM gathering thresholds
    grep -v '^#' ${hmmsearch_hits_file} | sort -k5,5g | head -n ${number_results_to_collect} > ${top_hits_file}
   
    # Add the top hits gene IDs and their scores to the run's collected IDs file
    while read -r line; do
        gene_id=$(echo $line | awk '{print $1}')
        score=$(echo $line | awk '{print $6}')
        printf "%-20s %s\n" "$gene_id" "$score"
        echo $gene_id >> ${collected_IDs_file}
    done < $top_hits_file
    
    # Add the top hits gene sequences to the current protein-database collected sequences file
	# Set up the protein database-specific hits fasta file
	current_species_top_hits_fasta_file="${current_species_search_results_dir}/${protein_database_name}_top_hits_sequences.fasta"
	touch $current_species_top_hits_fasta_file

    # Collect the fasta sequences for all top hits in the protein-database specific fasta file
    while IFS= read -r line; do
        gene_id=$(echo "$line" | awk '{print $1}')
        echo "Searching for gene ID: $gene_id"
        # Use awk to extract the sequence for each gene ID
        awk -v id="$gene_id" '
        BEGIN { found=0 }
        $0 ~ ">" && $0 ~ id { found=1; print; next }
        found && /^>/ { found=0 }
        found { print }
        ' $protein_database_file_path >> ${current_species_top_hits_fasta_file}
    done < ${top_hits_file}
    
    # Copy these to the run's collected hits fasta file
    cat ${current_species_top_hits_fasta_file} >> ${collected_seqs_file}
    
done; } 2>> ${timing_log} #entire for loop wrapped in timing log

# Append the following variables to the job_variables.txt file for downstream scripts to use

# HMMER specific files and directories
echo "collected_seqs_file=${collected_seqs_file}" >> ${output_var_file}
echo "collected_IDs_file=${collected_IDs_file}" >> ${output_var_file}
echo "hmmbuild_file=${hmmbuild_file}" >> ${output_var_file}


