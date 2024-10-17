#!/bin/bash
#SBATCH --partition=cpu        # Partition to use
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks=1                  # Number of tasks (processes)
#SBATCH --cpus-per-task=32           # Number of CPU cores per task
#SBATCH --mem=120G                   # Memory per node
#SBATCH -t 6:00:00
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err

# Load required modules
module load iq-tree/2.1.3-noMPI 
module load mafft/7.481

run_noisy=true  # or run_noisy=1

# Source the variables from the output_var_file
if [ -z "$output_var_file" ]; then
    echo "Error: output_var_file is not set."
    exit 1
fi
source $output_var_file


# Select the output file from previous step to become the tree's input
	output_from_previous_step=${extracted_domains_fasta_file} #If coming from the domain extraction script
	# output_from_previous_step=${collected_seqs_file} # If coming directly from the HMMsearch result
	
# copy the data from previous step output to the tree directory and rename
starting_data_file="${tree_inference_dir}/${run_name}_starting_data.fasta"
cp ${output_from_previous_step} ${starting_data_file}
echo "Input sequences to be used are copied from ${output_from_previous_step}"
echo "New file created in this directory with starting data, named ${starting_data_file}"

# clean up names in the starting data
starting_data_cleaned_file="${tree_inference_dir}/${run_name}_starting_data_cleaned.fasta"	

#process the names in the input file to prevent downstream problems with the tree inference software
echo "Processing gene names in starting data to remove characters that programs don't like"
awk '
    # If the line starts with ">", it is a name line
    /^>/ {
        # Modify the name by replacing spaces with underscores, removing square brackets, commas, and "putative_"
        name = $0
        gsub(/ /, "_", name)                 # Replace spaces with underscores
        gsub(/\[|\]/, "", name)              # Remove square brackets
        gsub(/,/, "", name)                  # Remove commas
        gsub(/putative_/, "", name)          # Remove the word "putative_"
        gsub(/:/, "_", name)                 # Replace colons with underscores
        print name
        next
    }
    # If the line does not start with ">", it is a sequence line, modify nothing
    {
        print
    }
' ${starting_data_file} > ${starting_data_cleaned_file}
echo "finished processing gene names. New file with cleaned names is ${starting_data_cleaned_file}"

# Add outgroup to the starting_data_cleaned_file
outgroup_name="Bnat_34821_MyoIV"
echo ">${outgroup_name}" >> $starting_data_cleaned_file
echo "MQRRFEQNKIYTNVGTILISVNPYQRLPLYTEQVLKKYTSRGLGVVDMPPHVFNIAHDAFYGVTSFSKGQSIIISGESGAGKTEATKQCLQYLAAIAGSTSDVEKKVLRANPILEAFGNAKTLRNDNSSRFGKYLEIYLDEKGRISSSATENYLLEKIRVVQPSLKERNFHIFYQLVKAASSKLRQKLKLKEDAGKYNYLKSCTDVPSIDDTRDYKEVIEAFRELGISDSEREQTFRICAAILHLGNCTFTDGKNHSGGCQVNEKAVLRDAAELLGVNGDKLLERLTTREIRVRGQSAAKAVMGAEEASDTRHALCKFVYGRMFDWIVARINKSMPGGGGGRSIGILDIFGFEIFEKNSFEQLCINFTNERLQQHFNRHTFKLEQNIYCSEGIDFDEIDYIDNQPMVDLITKKPHGVLPLLDEELRIPKGSDETFLAKLETKQCKNPVFKRQMKKRAHFAIKHYAGKVLYHCKGFLEKNRDTLTEDLVEILQTSNQPLLQELYPSDMQISSKQRKSSLATQFQQQLTRLMHSLNATQPHYIRCIKPNNDKAPMKFVAKNCHEQLLYSGVFEAVAIRKQGFPFRLSHEEFEKRYSICLGKTSALSNQSSVKGRCKLILNEMKCDPKNTRIGSSRVLYR" >> $starting_data_cleaned_file
echo "Outgroup ${outgroup_name} and its sequence was added to input sequences"

### Process starting_data_cleaned_file to prepare for tree
echo "Data prep: align input sequences"
starting_data_cleaned_aligned_file="${tree_inference_dir}/${run_name}_cleaned_aligned.fasta"	
# Check if input is an alignment, and align if it is not
if awk 'NR==2 {if (gsub(/-/, "&") >= 1) exit 1; else exit 0}' "$starting_data_cleaned_file"; then
    echo "Input data is not an alignment. Running MAFFT..."
    mafft --maxiterate 1000 --localpair --reorder --thread ${NUM_CORES} ${starting_data_cleaned_file} > $starting_data_cleaned_aligned_file
    echo "MAFFT complete. Output file is $starting_data_cleaned_aligned_file"
else
    echo "File is already an alignment. No action taken."
    starting_data_cleaned_aligned_file=$starting_data_cleaned_file #update with original input which was an alignment
fi

# Run noisy to improve alignment and remove positions without phylogenetic utility
echo "Running Noisy"
if [ "$run_noisy" = true ]; then
    # Clean alignment with noisy
    cd "${tree_inference_dir}"
    ~/noisy/bin/noisy "${starting_data_cleaned_aligned_file}" #Output MSA from noisy is basename(input) + _out.fas  No way to alter this in the arguments.
    starting_data_cleaned_aligned_noisyed_file="${starting_data_cleaned_aligned_file%.fasta}_out.fas"  # Should be the output from noisy
    echo "Noisy finished. Result is ${starting_data_cleaned_aligned_noisyed_file}"
else
    echo "Noisy set to not run. Alignment file will be passed directly to tree inference"
    starting_data_cleaned_aligned_noisyed_file=$starting_data_cleaned_aligned_file
fi





# Run the tree inference program
echo "Starting IQtree inference using ${starting_data_cleaned_aligned_noisyed_file}"
#	iqtree2 -s ${tree_input_clean_file} -fast -m LG+F+R7 -B 1000 --wbtl -nstop 500 #quick version
#	iqtree2 -s ${tree_input_clean_file} -m MFP -bb 1000 -bnni -o $outgroup_name #full version, infer site rates
#	iqtree2 -s ${tree_input_clean_file} -m Q.yeast+R7 -bb 1000 -bnni -nt ${NUM_CORES} -o $outgroup_name  #full version, force Q.yeast+R7 model because Modelfinder keeps deciding on this one
	iqtree2 -s ${starting_data_cleaned_aligned_noisyed_file} -m Q.yeast+R7 -bb 1000 -bnni -nt AUTO -o $outgroup_name  #full version, force Q.yeast+R7 model because Modelfinder keeps deciding on this one
echo "Finished IQtree program"

#remove iqtree log files
cd ${tree_dir}
echo "cleaning up.."
	rm *splits.nex
	rm *ckp.gz
	rm *bionj
	rm *mldist
#	rm *.log
	rm *.iqtree
	rm *model.gz
	rm *uniqueseq.phy
