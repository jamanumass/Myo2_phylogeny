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

#Name run after the input and create directory
protein_input_file="$(ls *_query.fasta)" # all inputs should end in "_query.fasta"
run_name=${protein_input_file%_query.fasta}
run_dir="search_for_${run_name}"
mkdir "${run_dir}"
collected_seqs_file="${run_dir}/${run_name}_collected_seqs.fa"
touch "${collected_seqs_file}"

#loop here
#Set the current protein genome being searched
protein_database_file="GCF_000004985.1_V1.0_protein.faa"  #replace this with for loop

#Set database directories based on the current working protein genome
protein_database_name="${protein_database_file%_protein.faa}"
protein_database_file_path="${protein_database_directory}/${protein_database_file}"
#cat $protein_database_file_path | head -n 10 #check that the file is called properly

#Make a new direcotry for the search results of the current working protein genome
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
    echo $gene_id >> all_hits_IDs.txt
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






