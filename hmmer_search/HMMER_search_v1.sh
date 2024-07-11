#!/bin/bash
#SBATCH --partition=gpu,gpu-preempt,gpu-long
#SBATCH --constraint=vram8
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-gpu=8
#SBATCH --mem-per-gpu=20G
#SBATCH -t 8:00:00
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err

module load iq-tree/2.1.3-noMPI 
module load mafft/7.481
module load uri/main
module load HMMER/3.3.2-iimpi-2021b


#Location of permanent databases
protein_database_directory="/work/pi_lfritzlaylin_umass_edu/users/jaman/sequence_databases/myo2_seqs_collect_v1/full_genomes_peptide_databases"

#define the results directory
results_dir="/work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny/hmmer_search/results"

#define the input directory (should contain sequence alignments that end in "_query.fasta")
input_dir="/work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny/hmmer_search/inputs"

#Name run after the input and create directory
protein_input_file="${input_dir}/$(ls *_query.fasta)" #search if only using 1 input
protein_input_name=${protein_input_file%_query.fasta} #remove _query.fasta
run_name=$(basename "$protein_input_name")  #remove directory info
run_dir="${results_dir}/search_for_${run_name}"
mkdir ${run_dir}

#set up files for collected IDs and sequences
collected_seqs_file="${run_dir}/${run_name}_collected_seqs.fa"
collected_IDs_file="${run_dir}/${run_name}_collected_IDs.txt"
touch "${collected_seqs_file}"
touch "${collected_IDs_file}"

#loop here through each protein database file and perform the search and sequence collection
for file_path in "${protein_database_directory}"/*_protein.faa; do #Set the current protein genome being searched
    protein_database_file=$(basename "$file_path")


    #Set database directories based on the current working protein genome
    protein_database_name="${protein_database_file%_protein.faa}"
    protein_database_file_path="${protein_database_directory}/${protein_database_file}"
    #cat $protein_database_file_path | head -n 10 #check that the file is called properly
    
    #Make a new sub-direcotry for the search results of the current working protein genome
    search_results_dir_name="search_results_${protein_database_name}"
    search_results_dir="${run_dir}/${search_results_dir_name}" 
    mkdir ${search_results_dir}
    
    #Dictate where the search result files will go and be named
    hmmbuild_file="${search_results_dir}/${run_name}.hmm"
    hmmsearch_hits_file="${search_results_dir}/${run_name}_in_${protein_database_name}_hmmsearch_out.txt"
    
    #build HMM and perform the search
    hmmbuild "${hmmbuild_file}" "${protein_input_file}"
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




