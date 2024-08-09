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
module load uri/main
module load HMMER/3.3.2-iimpi-2021b

# Set the root directory for Myo searches
root_dir="/work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny"
cd $root_dir

#create new folder for the files from the run, which is based on the input file
# Set the input name, which will also be the name of the folder in the results directory
input_file="Amoebozoan_MyoII_motors_input.fasta"
input_path="/work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny/input_seqs"
input_name="${input_file%_input.fasta}"
run_name=$input_name

#Location of permanent protein databases
protein_database_directory="/work/pi_lfritzlaylin_umass_edu/users/jaman/sequence_databases/myo2_seqs_collect_v1/full_genomes_peptide_databases"

#define the run directory 
run_dir="/work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny/results_${input_name}"

#define the hmmsearch directory
hmm_search_dir="${run_dir}/hmm_search_results"







# Find the collected seqs file
collected_seqs_file=$(ls ${hmm_search_dir}/*_collected_seqs.fa)
    #debug: print lines in the collected seqs file
    echo "Lines in collected seqs file:"
    cat $collected_seqs_file | wc -l

# Make a new folder to contain the results from the domain extraction process
extraction_dir="${run_dir}/domain_extraction_result"
mkdir -p $extraction_dir

##this section to pull previously made HMM files to extract the portion of results fitting the HMM
# Set the HMM file that describes just the motor domain
domain_hmm_dir="/work/pi_lfritzlaylin_umass_edu/users/jaman/Pfam_HMMs"
domain_hmm_file_name="myosin_head_motor_domain.hmm"
domain_hmm_file="${domain_hmm_dir}/${domain_hmm_file_name}"

# Extract the collected seqs and HMM base names without the directories
collected_seqs_base=$(basename "$collected_seqs_file" .fa)
domain_hmm_base=$(basename "$domain_hmm_file" .hmm)

# Define the output file
extraction_output_file="$extraction_dir/${collected_seqs_base}_${domain_hmm_base}_extraction.fa"

#1 - Generate a domain table with the target of each hit identified in its sequence
#Run hmmsearch and extract the sequences
hmmsearch --domtblout "${extraction_output_file}.domtblout" "$domain_hmm_file" "$collected_seqs_file"
#debug: how many hits found. Should be about 5x the number of genomes
echo "number of lines in HMMsearch domain results"
cat ${extraction_output_file}.domtblout | wc -l


#2 Generate a simplified output file of the domain coordinates
#Define the simplified output file
simplified_output_file="${extraction_output_file}.simplified"

# Process the domain table to extract columns 1 (target name), 20 (env from), and 21 (env to) to make simplified table
awk '$1 !~ /^#/ { print $1, $20, $21 }' "${extraction_output_file}.domtblout" > "$simplified_output_file"
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

# Debug: Check the contents of the merged file
echo "number of lines in the gene-merged table"
cat "$merged_output_file" | wc -l


#4 calculate and write the lengths of each sequence
# Define the output file with length
length_output_file="${merged_output_file}.length"

# Add the length as the 4th column
awk '{ print $1, $2, $3, $3 - $2 + 1 }' "$merged_output_file" > "$length_output_file"

# Debug: Check the contents of the length file
echo "Contents of length_output_file: [geneID] [domain start] [domain end] [length]"
head -n 2 "$length_output_file"


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

# Calculate 80% of the median length
eighty_percent_median=$(echo "$median_length * 0.8" | bc)

# Calculate 120% of the median length
one_twenty_percent_median=$(echo "$median_length * 1.2" | bc)

# Define the filtered output file
filtered_output_file="${length_output_file}.filtered"

# Filter lines with length more than 80% of the median length and less than 120% of the median length
awk -v min_threshold="$eighty_percent_median" -v max_threshold="$one_twenty_percent_median" '$4 > min_threshold && $4 < max_threshold { print $0 }' "$length_output_file" > "$filtered_output_file"
#debug: filtered output file should have fewer lines
echo "filtered output file with too long or too short seqs removed, lines"
cat $filtered_output_file | wc -l





#6 Extract the portion of sequences defined by the domains
# Define the final fasta output file
hits_domain_extraction_file="${extraction_dir}/${input_name}_hits_domain_extraction.fa"
touch $hits_domain_extraction_file

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
    ' "$collected_seqs_file" >> "$hits_domain_extraction_file"
done < "$filtered_output_file"



# Cleanup intermediate files
rm "$simplified_output_file" "$merged_output_file" "$length_output_file" "$filtered_output_file"


