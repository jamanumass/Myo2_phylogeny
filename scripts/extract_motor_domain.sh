#!/bin/bash
#SBATCH --partition=gpu,gpu-preempt,gpu-long
#SBATCH --constraint=vram8
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-gpu=8
#SBATCH --mem-per-gpu=20G
#SBATCH -t 8:00:00
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err


#This script is designed to perform domain extraction and sequence alignment on a collection of sequences discovered using HMMER. 


#load required modules
module load iq-tree/2.1.3-noMPI 
module load mafft/7.481
module load uri/main
module load HMMER/3.3.2-iimpi-2021b


#working directory is always started at Myo2 project /work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny
root_dir="/work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny"
cd $root_dir

#Set the project directory
proj_dir="${root_dir}/hmmer_search"

#Set the input name, which will also be the name of the folder in results directory
input_dir="search_for_Amoebozoan_MyoII_motors"

#Find the collected seqs file
run_dir="${proj_dir}/results/${input_dir}"
collected_seqs_file=$(ls $run_dir/*_collected_seqs.fa)

#make new folder to contain the results from the domain extraction process
extraction_dir="${run_dir}/domain_extraction_result"
mkdir $extraction_dir

#set the HMM file that describes just the motor domain
domain_hmm_dir="/work/pi_lfritzlaylin_umass_edu/users/jaman/Pfam_HMMs"
domain_hmm_file_name="myosin_head_motor_domain.hmm"
domain_hmm_file="${domain_hmm_dir}/${domain_hmm_file_name}"



# Extract the collected seqs and HMM base names without the directories
collected_seqs_base=$(basename "$collected_seqs_file" .fa)
domain_hmm_base=$(basename "$domain_hmm_file" .hmm)

# Define the output file
extraction_output_file="$extraction_dir/${collected_seqs_base}_${domain_hmm_base}_extraction.fa"

# Run hmmsearch and extract the sequences
hmmsearch --domtblout "${extraction_output_file}.domtblout" "$domain_hmm_file" "$collected_seqs_file"

#Parse the output and extract the domains found within each of the sequences
#Reads specific fields (columns) from the domtblout file
#Removes duplicates.
#For each unique entry, it extracts the corresponding sequence segment from another file.
#The extracted sequences are saved to the specified output file.
awk '$1 !~ /^#/ { print $1, $20, $21 }' "${extraction_output_file}.domtblout" | sort | uniq | while read -r seqid start end; do
    awk -v seqid="$seqid" -v start="$start" -v end="$end" '
    BEGIN { found = 0; header = ""; seq = "" }
    /^>/ {
        if (found) { exit }
        if ($0 ~ ">" seqid) { found = 1; header = $0 }
    }
    found && !/^>/ {
        seq = seq $0
    }
    END {
        if (found) {
            print header
            print substr(seq, start, end - start + 1)
        }
    }
    ' "$collected_seqs_file"
done > "${extraction_output_file}.temp"

echo "Domain extraction complete. The sequences are saved in "${extraction_output_file}.temp"



####Output file can have duplicates. Remove. 
# Define input and output files
input_file="${extraction_output_file}.temp"
output_file="${extraction_output_file}"

# Create a temporary file to store unique sequences
temp_file=$(mktemp)

# Initialize an associative array to track seen sequences
declare -A seen

# Process the input file
{
  while read -r line; do
    # If the line starts with '>', it's a header line
    if [[ $line == '>'* ]]; then
      header="$line"
    else
      # It's a sequence line
      sequence="$line"
      # Check if the sequence has been seen before
      if [[ -z "${seen[$sequence]}" ]]; then
        # If not, mark it as seen and print the header and sequence to the temp file
        seen[$sequence]=1
        echo "$header" >> "$temp_file"
        echo "$sequence" >> "$temp_file"
      fi
    fi
  done
} < "$input_file"

# Move the temp file to the output file
mv "$temp_file" "${extraction_output_file}"

echo "Duplicates removed. Unique sequences saved to ${extraction_output_file}"



#align the extracted domains
extraction_output_file_name="${extraction_output_file%.fa}"
aligned_file="${extraction_output_file_name}_aligned.fa"
mafft ${extraction_output_file} > ${aligned_file}