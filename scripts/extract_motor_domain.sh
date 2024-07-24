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



# Set the root directory
root_dir="/work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny"
cd $root_dir

# Set the project directory
proj_dir="${root_dir}/hmmer_search"

# Set the input name, which will also be the name of the folder in the results directory
input_dir="search_for_Amoebozoan_MyoII_motors"

# Find the collected seqs file
run_dir="${proj_dir}/results/${input_dir}"
collected_seqs_file=$(ls $run_dir/*_collected_seqs.fa)

# Make a new folder to contain the results from the domain extraction process
extraction_dir="${run_dir}/domain_extraction_result"
mkdir -p $extraction_dir

# Set the HMM file that describes just the motor domain
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

# Define the simplified output file
simplified_output_file="${extraction_output_file}.simplified"

# Process the domain table to extract columns 1 (target name), 20 (env from), and 21 (env to)
awk '$1 !~ /^#/ { print $1, $20, $21 }' "${extraction_output_file}.domtblout" > "$simplified_output_file"

# Define the merged output file
merged_output_file="${simplified_output_file}.merged"

# Merge lines with the same target in column 1, taking the lower value of column 2 and the higher value of column 3
awk '
{
    if (!($1 in min_start)) {
        min_start[$1] = $2;
        max_end[$1] = $3;
    } else {
        if ($2 < min_start[$1]) min_start[$1] = $2;
        if ($3 > max_end[$1]) max_end[$1] = $3;
    }
}
END {
    for (target in min_start) {
        print target, min_start[target], max_end[target];
    }
}' "$simplified_output_file" > "$merged_output_file"

# Debug: Check the contents of the merged file
echo "Contents of $merged_output_file:"
head -n 20 "$merged_output_file"

# Define the output file with length
length_output_file="${merged_output_file}.length"

# Add the length as the 4th column
awk '{ print $1, $2, $3, $3 - $2 + 1 }' "$merged_output_file" > "$length_output_file"

# Debug: Check the contents of the length file
echo "Contents of $length_output_file:"
head -n 20 "$length_output_file"






## Calculate the median length from column 4
median_length=$(awk '{print $4}' "$length_output_file" | sort -n | awk '{
    count[NR] = $1;
}
END {
    if (NR % 2 == 1) {
        print count[(NR + 1) / 2];
    } else {
        print (count[(NR / 2)] + count[(NR / 2) + 1]) / 2;
    }
}')

# Calculate 80% of the median length
eighty_percent_median=$(echo "$median_length * 0.8" | bc)

# Calculate 120% of the median length
one_twenty_percent_median=$(echo "$median_length * 1.2" | bc)

# Define the filtered output file
filtered_output_file="${length_output_file}.filtered"

# Filter lines with length more than 80% of the median length and less than 120% of the median length
awk -v min_threshold="$eighty_percent_median" -v max_threshold="$one_twenty_percent_median" '$4 > min_threshold && $4 < max_threshold { print $0 }' "$length_output_file" > "$filtered_output_file"




# Extract and trim sequences based on the filtered file
while read -r seqid start end length; do
    # Extract the full sequence for the current seqid
    awk -v seqid="$seqid" -v start="$start" -v end="$end" '
    BEGIN { found = 0 }
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
    ' "$collected_seqs_file" >> "$final_fasta_output_file"
done < "$filtered_output_file"

# Cleanup intermediate files
rm "$simplified_output_file" "$merged_output_file" "$length_output_file" "$filtered_output_file"

#clean up name of output fasta file
mv "$final_fasta_output_file" "$(dirname "$final_fasta_output_file")/$(basename "${final_fasta_output_file%.fa.simplified.merged.length.filtered.fasta}.fasta")"





#align the extracted domains
extraction_output_file_name="${extraction_output_file%.fa}"
aligned_file="${extraction_output_file_name}_aligned.fa"
mafft ${extraction_output_file} > ${aligned_file}