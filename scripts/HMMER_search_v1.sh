#!/bin/bash
#SBATCH --partition=gpu-long         # Select partition
#SBATCH --gpus=1                     # Request 1 GPU
#SBATCH --cpus-per-gpu=16            # 16 CPU cores per GPU
#SBATCH --mem-per-gpu=20G            # 20 GB memory per GPU
#SBATCH -t 48:00:00
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err

module load uri/main
module load HMMER/3.3.2-iimpi-2021b
module load mafft/7.481


# Set the root directory for Myo searches
root_dir="/work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny"
cd $root_dir




# Check if INPUT_FASTA environment variable is set
if [ -z "$INPUT_FASTA" ]; then
  echo "Error: INPUT_FASTA environment variable is not set."
  exit 1
fi
  
run_name=$(basename $INPUT_FASTA | sed 's/_input\.fasta$//')

echo "Input fasta file is $INPUT_FASTA"
echo "Run name is $run_name"



#check if input is an alignment, and align if it is not
if awk 'NR==2 {if (gsub(/-/, "&") >= 5) exit 1; else exit 0}' "$INPUT_FASTA"; then
    echo "Not an alignment. Running MAFFT..."
    mafft --maxiterate 1000 --localpair --reorder "$INPUT_FASTA" > "${INPUT_FASTA}.aligned"
    cat "${INPUT_FASTA}.aligned" > "$INPUT_FASTA"
    rm "${INPUT_FASTA}.aligned"
else
    echo "File is already an alignment. No action taken."
fi


#Location of permanent protein databases
protein_database_directory="/work/pi_lfritzlaylin_umass_edu/users/jaman/sequence_databases/myo2_seqs_collect_v1/full_genomes_peptide_databases"

#define the run directory 
run_dir="/work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny/results_${run_name}"
mkdir ${run_dir}



#define the hmmsearch directory
hmm_search_dir="${run_dir}/hmm_search_results"
mkdir ${hmm_search_dir}


#set up files for collected IDs and sequences
collected_seqs_file="${hmm_search_dir}/${run_name}_collected_seqs.fa"
collected_IDs_file="${hmm_search_dir}/${run_name}_collected_IDs.txt"
touch "${collected_seqs_file}"
touch "${collected_IDs_file}"
number_results_to_collect=10# define the number of results from each protein database to collect

#Create the HMM for the input
hmmbuild_file="${hmm_search_dir}/${run_name}.hmm"
hmmbuild "${hmmbuild_file}" "${INPUT_FASTA}"

#loop here through each protein database file and perform the search and sequence collection
for file_path in "${protein_database_directory}"/*_protein.faa; do #Set the current protein genome being searched
    protein_database_file=$(basename "$file_path")


    #Set database directories based on the current working protein genome
    protein_database_name="${protein_database_file%_protein.faa}"
    protein_database_file_path="${protein_database_directory}/${protein_database_file}"

    #Make a new sub-directory for the search results of the current working protein genome
    search_results_dir_name="search_results_${protein_database_name}"
    search_results_dir="${hmm_search_dir}/${search_results_dir_name}" 
    mkdir ${search_results_dir}
    
    #Dictate where the search result files will go and be named
    hmmsearch_hits_file="${search_results_dir}/${run_name}_in_${protein_database_name}_hmmsearch_out.txt"
    
    #Perform the HMM search
    hmmsearch --tblout "${hmmsearch_hits_file}" "${hmmbuild_file}" "${protein_database_file_path}"
    
    # Extract the best 5 hits from the output table and save to best5 file
    #define file name
    best5_file="${search_results_dir}/${protein_database_name}_best5.txt"
    
    #create the best5 file
    grep -v '^#' ${hmmsearch_hits_file} | sort -k5,5g | head -n 5 > ${best5_file}
    
   
    # Output the list of gene IDs and their scores
    echo "Top 5 HMMsearch Hits:"
    echo "---------------------"
    printf "%-20s %s\n" "Gene Name" "Full Gene Score"
    echo "---------------------"
   
    while read -r line; do
        gene_id=$(echo $line | awk '{print $1}')
        score=$(echo $line | awk '{print $6}')
        printf "%-20s %s\n" "$gene_id" "$score"
        echo $gene_id >> ${collected_IDs_file}
    done < $best5_file
    
    #collect the sequence of the best5 hits
    best5_fasta_file="${search_results_dir}/${protein_database_name}_best5_sequences.fasta"
    
    # Read each line from the best5 file and extract the gene ID (first column)
    while IFS= read -r line; do
        gene_id=$(echo "$line" | awk '{print $1}')
        echo "Searching for gene ID: $gene_id"
        # Use awk to extract the sequence for each gene ID
        awk -v id="$gene_id" '
        BEGIN { found=0 }
        $0 ~ ">" && $0 ~ id { found=1; print; next }
        found && /^>/ { found=0 }
        found { print }
        ' $protein_database_file_path >> ${best5_fasta_file}
    done < ${best5_file}
    
    echo "Sequences have been extracted to ${best5_fasta_file}"
    cat ${best5_fasta_file} >> ${collected_seqs_file}
    
done




