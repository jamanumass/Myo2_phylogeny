#!/bin/bash
#SBATCH --partition=cpu        # Partition to use
#SBATCH --nodes=1              # Number of nodes
#SBATCH --ntasks=1             # Number of tasks (processes)
#SBATCH --cpus-per-task=32     # Number of CPU cores per task
#SBATCH --mem=120G             # Memory per node
#SBATCH -t 6:00:00
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err

# Load required modules
module load iq-tree/2.1.3-noMPI 
module load MAFFT/7.505-GCC-11.3.0-with-extensions

# Set up file paths
starting_data="${tree_inference_dir}/${run_name}_starting_data.fasta"
cleaned_data="${tree_inference_dir}/${run_name}_starting_data_cleaned.fasta"
aligned_data="${tree_inference_dir}/${run_name}_cleaned_aligned.fasta"
deduplicated_data="${tree_inference_dir}/${run_name}_deduplicated.fasta"
final_data="$aligned_data"

# Copy input data for processing
cp ${extracted_domains_fasta_file} $starting_data
echo "Copied input sequences from ${extracted_domains_fasta_file} to ${starting_data}"

# Deduplicate sequences using CD-HIT (remove sequences >95% similar)
echo "Running CD-HIT to remove sequences with >95% similarity..."
cd-hit -i ${starting_data} -o ${deduplicated_data} -c 0.95
echo "Deduplication complete. Output file is ${deduplicated_data}"

# Clean up sequence names for downstream compatibility
echo "Cleaning up gene names..."
awk '
    /^>/ {gsub(/ /, "_"); gsub(/\[|\]/, ""); gsub(/,/, ""); gsub(/putative_/, ""); gsub(/:/, "_"); print; next}
    {print}
' ${deduplicated_data} > ${cleaned_data}
echo "Name cleaning complete. Cleaned data saved to ${cleaned_data}"

# Add outgroup to the cleaned data
outgroup_name="Bnat_34821_MyoIV"
echo ">${outgroup_name}" >> $cleaned_data
echo "MQRRFEQNKIYTNVGTILISVNPYQRLPLYTEQVLKKYTSRGLGVVDMPPHVFNIAHDAFYGVTSFSKGQSIIISGESGAGKTEATKQCLQYLAAIAGSTSDVEKKVLRANPILEAFGNAKTLRNDNSSRFGKYLEIYLDEKGRISSSATENYLLEKIRVVQPSLKERNFHIFYQLVKAASSKLRQKLKLKEDAGKYNYLKSCTDVPSIDDTRDYKEVIEAFRELGISDSEREQTFRICAAILHLGNCTFTDGKNHSGGCQVNEKAVLRDAAELLGVNGDKLLERLTTREIRVRGQSAAKAVMGAEEASDTRHALCKFVYGRMFDWIVARINKSMPGGGGGRSIGILDIFGFEIFEKNSFEQLCINFTNERLQQHFNRHTFKLEQNIYCSEGIDFDEIDYIDNQPMVDLITKKPHGVLPLLDEELRIPKGSDETFLAKLETKQCKNPVFKRQMKKRAHFAIKHYAGKVLYHCKGFLEKNRDTLTEDLVEILQTSNQPLLQELYPSDMQISSKQRKSSLATQFQQQLTRLMHSLNATQPHYIRCIKPNNDKAPMKFVAKNCHEQLLYSGVFEAVAIRKQGFPFRLSHEEFEKRYSICLGKTSALSNQSSVKGRCKLILNEMKCDPKNTRIGSSRVLYR" >> $cleaned_data
echo "Outgroup ${outgroup_name} and its sequence added to ${cleaned_data}"

# Align sequences (MAFFT)
if awk 'NR==2 {if (gsub(/-/, "&") >= 1) exit 1; else exit 0}' "$cleaned_data"; then
    echo "Input data is not an alignment. Running MAFFT..."
    mafft --maxiterate 1000 --localpair --reorder --thread ${NUM_CORES} ${cleaned_data} > $aligned_data
    final_data=$aligned_data
    echo "MAFFT complete. Output file is $aligned_data"
else
    echo "File is already aligned. No further alignment needed."
    final_data=$cleaned_data
fi

run_noisy="true"
# Optional: Run Noisy to improve alignment
if [ "$run_noisy" = true ]; then
    echo "Running Noisy..."
    noisy "$final_data"
    final_data="${final_data%.fasta}_out.fas"
    echo "Noisy finished. Final aligned data: $final_data"
fi

# Tree inference using IQ-TREE
echo "Starting IQ-TREE inference using ${final_data}"
iqtree2 -s ${final_data} -m Q.yeast+R7 -bb 1000 -bnni -nt AUTO -o $outgroup_name
echo "IQ-TREE inference complete"