#!/bin/bash
#SBATCH --partition=cpu        # Partition to use
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks=1                  # Number of tasks (processes)
#SBATCH --cpus-per-task=32           # Number of CPU cores per task
#SBATCH --mem=120G                   # Memory per node
#SBATCH -t 01:00:00
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err


#This script is designed to perform domain extraction and sequence alignment on a collection of sequences discovered using HMMER. 


#load required modules
module load uri/main
module load HMMER/3.3.2-iimpi-2021b


# Source the variables from the output_var_file
if [ -z "$output_var_file" ]; then
    echo "Error: output_var_file is not set."
    exit 1
fi
source $output_var_file





##this section to pull previously made HMM files to extract the portion of results fitting the HMM
# Set the HMM file that describes just the myosin motor domain
domain_hmm_dir="/work/pi_lfritzlaylin_umass_edu/users/jaman/Pfam_HMMs"
domain_hmm_file_name="myosin_head_motor_domain.hmm"
domain_hmm_file="${domain_hmm_dir}/${domain_hmm_file_name}"

# name of the HMM with the domain to be used
domain_hmm_base_name=$(basename "$domain_hmm_file" .hmm)

# Define the output file
domain_table_file="$domain_extraction_dir/${run_name}_${domain_hmm_base_name}.domtblout"

#1 - Generate a domain table with the target of each hit identified in its sequence
#Run hmmsearch and extract the sequences
hmmsearch --domtblout ${domain_table_file} --noali ${domain_hmm_file} ${collected_seqs_file}

#2 Generate a simplified output file of the domain coordinates
#Define the simplified output file
simplified_output_file="${domain_table_file}.simplified"

# Process the domain table to extract columns 1 (target name), 20 (env from), and 21 (env to) to make simplified table
awk '$1 !~ /^#/ { print $1, $20, $21 }' ${domain_table_file} > ${simplified_output_file}
#debug count the number of lines in the simplified table
echo "Number of lines in the simplified domain table:"
cat $simplified_output_file | wc -l

#3 Merge any domains that came from the same gene into one coordinate set
#Define the merged output file
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



#4 calculate and write the lengths of each sequence
# Define the output file with length
length_output_file="${merged_output_file}.length"

# Add the length as the 4th column
awk '{ print $1, $2, $3, $3 - $2 + 1 }' "$merged_output_file" > "$length_output_file"



#5 Remove sequences outside of length parameters to produce a filtered list of sequences
#Calculate the median length all sequences
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

# Calculate lower length threshold
lower_length_threshold=$(echo "$median_length * 0.6" | bc)

# Calculate upper length threshold
upper_length_threshold=$(echo "$median_length * 1.4" | bc)

# Define the filtered output file
filtered_output_file="${length_output_file}.filtered"

# remove lines with outside of the median length thresholds
awk -v min_threshold="$lower_length_threshold" -v max_threshold="$upper_length_threshold" '$4 > min_threshold && $4 < max_threshold { print $0 }' "$length_output_file" > "$filtered_output_file"
#debug: filtered output file should have fewer lines
echo "filtered output file with too long or too short seqs removed, lines"
cat $filtered_output_file | wc -l




#6 Extract the portion of sequences defined by the domains
# Define the final fasta output file
extracted_domains_fasta_file="${domain_extraction_dir}/${run_name}_${domain_hmm_base_name}_extraction.fasta"
touch $extracted_domains_fasta_file

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
    ' "$collected_seqs_file" >> "$extracted_domains_fasta_file"
done < "$filtered_output_file"

# Cleanup intermediate files
rm "$simplified_output_file" "$merged_output_file" "$length_output_file" "$filtered_output_file"

# Report 
echo "number of lines in the domain-extraction collected seqs file, should be 2X the number in the filtered output file:"
cat ${extracted_domains_fasta_file} | wc -l

# Append the following variables to the job_variables.txt file for downstream scripts to use
echo "extracted_domains_fasta_file=${extracted_domains_fasta_file}" >> ${output_var_file}
