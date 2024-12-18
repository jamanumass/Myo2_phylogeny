#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

# Load required modules
module load uri/main
module load mpich/4.2.1 iq-tree/2.3.1
module load MAFFT/7.505-GCC-11.3.0-with-extensions
module load gcc/9.4.0

# Source the variable file to load all variables into the environment
if [ -f "$shared_variables_file" ]; then
    source "$shared_variables_file"
else
    echo "Error: Variable file not found."
    exit 1
fi

# ==============================
# Define File Variables
# ==============================
input1_seqs="${tree_inference_dir}/${run_name}_input1_seqs.fasta"
input2_cleaned_seqs="${tree_inference_dir}/${run_name}_input2_cleaned_seqs.fasta"
input3_deduplicated_seqs="${tree_inference_dir}/${run_name}_input3_deduplicated_seqs.fasta"
input4_aligned="${tree_inference_dir}/${run_name}_input4_aligned.fasta"
final_input_ready_for_iqtree=""

# ==============================
# Input File Selection
# ==============================
if [ -n "$extracted_domains_fasta_file" ] && [ -f "$extracted_domains_fasta_file" ]; then
    starting_data="$extracted_domains_fasta_file"
elif [ -n "$hmm_collected_seqs_file" ] && [ -f "$hmm_collected_seqs_file" ]; then
    starting_data="$hmm_collected_seqs_file"
elif [ -n "$INPUT_FASTA" ] && [ -f "$INPUT_FASTA" ]; then
    starting_data="$INPUT_FASTA"
else
    echo "Error: No valid input files found."
    exit 1
fi

# Copy the selected starting file to the tree inference directory
cp "$starting_data" "$input1_seqs"
echo "Starting from file: $starting_data"
echo "starting data file has $(grep -c ">" $starting_data) sequences"
echo "Copying to input1 file $input1_seqs"
echo "Genes in the input file: $(grep -c ">" $input1_seqs)"

# ==============================
# Cleanup Gene Names and Sequences
# ==============================
if [ ! -f "$input2_cleaned_seqs" ]; then
    echo "Cleaning gene names and sequences..."
    /home/jaman_umass_edu/jaman/general_tools/clean_fasta_names_and_seqs.sh "$input1_seqs" "$input2_cleaned_seqs"
    if [ ! -s "$input2_cleaned_seqs" ]; then
        echo "Error: Gene name and sequence cleanup failed."
        exit 1
    fi
    echo "Gene name and sequence cleanup complete: Genes in the cleaned file: $(grep -c ">" $input2_cleaned_seqs)"

else
    echo "Cleaned sequences file $input2_cleaned_seqs already exists. Skipping cleanup."
fi

# ==============================
# Add Outgroup Sequence
# ==============================
outgroup_name="Bnat_34821_MyoIV"
if ! grep -q ">${outgroup_name}" "$input2_cleaned_seqs"; then
    echo "Adding outgroup sequence..."
    echo ">${outgroup_name}" >> "$input2_cleaned_seqs"
    echo "MQRRFEQNKIYTNVGTILISVNPYQRLPLYTEQVLKKYTSRGLGVVDMPPHVFNIAHDAFYGVTSFSKGQSIIISGESGAGKTEATKQCLQYLAAIAGSTSDVEKKVLRANPILEAFGNAKTLRNDNSSRFGKYLEIYLDEKGRISSSATENYLLEKIRVVQPSLKERNFHIFYQLVKAASSKLRQKLKLKEDAGKYNYLKSCTDVPSIDDTRDYKEVIEAFRELGISDSEREQTFRICAAILHLGNCTFTDGKNHSGGCQVNEKAVLRDAAELLGVNGDKLLERLTTREIRVRGQSAAKAVMGAEEASDTRHALCKFVYGRMFDWIVARINKSMPGGGGGRSIGILDIFGFEIFEKNSFEQLCINFTNERLQQHFNRHTFKLEQNIYCSEGIDFDEIDYIDNQPMVDLITKKPHGVLPLLDEELRIPKGSDETFLAKLETKQCKNPVFKRQMKKRAHFAIKHYAGKVLYHCKGFLEKNRDTLTEDLVEILQTSNQPLLQELYPSDMQISSKQRKSSLATQFQQQLTRLMHSLNATQPHYIRCIKPNNDKAPMKFVAKNCHEQLLYSGVFEAVAIRKQGFPFRLSHEEFEKRYSICLGKTSALSNQSSVKGRCKLILNEMKCDPKNTRIGSSRVLYR" >> "$input2_cleaned_seqs"
    echo "Outgroup ${outgroup_name} and its sequence added."
else
    echo "Outgroup ${outgroup_name} is already present in the cleaned sequences."
fi

# ==============================
# Sequence Deduplication
# ==============================
if [ ! -f "$input3_deduplicated_seqs" ]; then
    similarity_threshold=0.95
    echo "Running duplicate sequence removal..."

    # Call the remove_duplicate_sequences_fasta.sh script to remove duplicate gene names
    /home/jaman_umass_edu/jaman/general_tools/remove_duplicate_sequences_fasta.sh "$input2_cleaned_seqs"

    # Check if the deduplication was successful by verifying that the file exists
    if [ ! -s "$input2_cleaned_seqs" ]; then
        echo "Error: Duplicate removal failed or no sequences left."
        exit 1
    fi

    echo "Running CD-HIT to remove any sequences with >${similarity_threshold} similarity..."

    export PATH=/home/jaman_umass_edu/bin:$PATH
    export LD_LIBRARY_PATH=/home/jaman_umass_edu/software/cdhit-master/lib:$LD_LIBRARY_PATH
    export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH

    /home/jaman_umass_edu/bin/cd-hit -i "$input2_cleaned_seqs" -o "$input3_deduplicated_seqs" -c "$similarity_threshold"

    # Check if CD-HIT deduplication was successful
    if [ ! -s "$input3_deduplicated_seqs" ]; then
        echo "Error: Deduplication failed."
        exit 1
    fi

    echo "Deduplication complete: Genes in the deduplicated file: $(grep -c ">" $input3_deduplicated_seqs)"
else
    echo "Deduplicated sequences file $input3_deduplicated_seqs already exists. Skipping deduplication."
fi

# ==============================
# Sequence Alignment
# ==============================
if [ ! -f "$input4_aligned" ]; then
    if awk 'NR==2 {if (gsub(/-/, "&") >= 1) exit 1; else exit 0}' "$input3_deduplicated_seqs"; then
        echo "Input data is not an alignment. Running MAFFT..."
        mafft --quiet --maxiterate 1000 --localpair --reorder --thread 16 "$input3_deduplicated_seqs" > "$input4_aligned"
        final_input_ready_for_iqtree="$input4_aligned"
        echo "Alignment complete: Genes in the aligned file: $(grep -c ">" $input4_aligned)"
    else
        echo "File is already aligned. No further alignment needed."
        final_input_ready_for_iqtree="$input3_deduplicated_seqs"
    fi
else
    echo "Aligned sequences file $input4_aligned already exists. Skipping alignment."
    final_input_ready_for_iqtree="$input4_aligned"
fi

# ==============================
# Tree Inference with IQ-TREE
# ==============================
echo "Starting IQ-TREE inference..."
echo "Genes in the tree input file: $(grep -c ">" $final_input_ready_for_iqtree)"
# Create a subdirectory for IQ-TREE analysis files
iqtree_dir="${tree_inference_dir}/iqtree_analysis"
mkdir -p "$iqtree_dir"
cd "$iqtree_dir" || exit #move to the directory or exit if this fails

# Run IQ-TREE on the final input file
iqtree2 -s "$final_input_ready_for_iqtree" -m Q.yeast+R7 -bb 1000 -bnni -nt ${NUM_CORES} -o "$outgroup_name" --keep-ident

echo "IQ-TREE inference complete. Results are in $iqtree_dir"
