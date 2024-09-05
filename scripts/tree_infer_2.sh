#!/bin/bash
#SBATCH --partition=cpu        # Partition to use
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks=1                  # Number of tasks (processes)
#SBATCH --cpus-per-task=32           # Number of CPU cores per task
#SBATCH --mem=120G                   # Memory per node
SBATCH -t 4:00:00
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err

# Load required modules
module load iq-tree/2.1.3-noMPI 
module load mafft/7.481


# Source the variables from the output_var_file
if [ -z "$output_var_file" ]; then
    echo "Error: output_var_file is not set."
    exit 1
fi
source $output_var_file

# Select the input file for the tree inference and copy to the tree inference directory
	#If coming from the domain extraction script
	tree_input=${extracted_domains_fasta_file}
	# If coming directly from the HMMsearch result
	# tree_input=${collected_seqs_file}
	
tree_input_file="${tree_inference_dir}/${run_name}_tree_input_seqs.fasta"	

cp ${tree_input} ${tree_input_file}

# Make a clean name version of the tree input file 
tree_input_clean_file="${tree_inference_dir}/${run_name}_input_clean.fasta"	

#process the names in the input file to prevent downstream problems with the tree inference software
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
' ${tree_input_file} > ${tree_input_clean_file}

# Add outgroup to the clean tree input file
outgroup_name="Bnat_34821_MyoIV"
echo ">${outgroup_name}" >> $tree_input_clean_file
echo "MQRRFEQNKIYTNVGTILISVNPYQRLPLYTEQVLKKYTSRGLGVVDMPPHVFNIAHDAFYGVTSFSKGQSIIISGESGAGKTEATKQCLQYLAAIAGSTSDVEKKVLRANPILEAFGNAKTLRNDNSSRFGKYLEIYLDEKGRISSSATENYLLEKIRVVQPSLKERNFHIFYQLVKAASSKLRQKLKLKEDAGKYNYLKSCTDVPSIDDTRDYKEVIEAFRELGISDSEREQTFRICAAILHLGNCTFTDGKNHSGGCQVNEKAVLRDAAELLGVNGDKLLERLTTREIRVRGQSAAKAVMGAEEASDTRHALCKFVYGRMFDWIVARINKSMPGGGGGRSIGILDIFGFEIFEKNSFEQLCINFTNERLQQHFNRHTFKLEQNIYCSEGIDFDEIDYIDNQPMVDLITKKPHGVLPLLDEELRIPKGSDETFLAKLETKQCKNPVFKRQMKKRAHFAIKHYAGKVLYHCKGFLEKNRDTLTEDLVEILQTSNQPLLQELYPSDMQISSKQRKSSLATQFQQQLTRLMHSLNATQPHYIRCIKPNNDKAPMKFVAKNCHEQLLYSGVFEAVAIRKQGFPFRLSHEEFEKRYSICLGKTSALSNQSSVKGRCKLILNEMKCDPKNTRIGSSRVLYR" >> $tree_input_clean_file 


# Check if input is an alignment, and align if it is not
if awk 'NR==2 {if (gsub(/-/, "&") >= 5) exit 1; else exit 0}' "$tree_input_clean_file"; then
    echo "Not an alignment. Running MAFFT..."
    #make a temp file to store the alignment
    touch "${tree_input_clean_file}.temp"
    mafft --maxiterate 1000 --localpair --reorder --thread ${NUM_CORES} ${tree_input_clean_file} > "${tree_input_clean_file}.temp"
    cat "${tree_input_clean_file}.temp" > $tree_input_clean_file  # Overwrite input file with the alignment
    rm "${tree_input_clean_file}.temp"
else
    echo "File is already an alignment. No action taken."
fi

# Run the tree inference program
#	iqtree2 -s ${tree_input_clean_file} -fast -m LG+F+R7 -B 1000 --wbtl -nstop 500 #quick version
#	iqtree2 -s ${tree_input_clean_file} -m MFP -bb 1000 -bnni -o $outgroup_name #full version, infer site rates
	iqtree2 -s ${tree_input_clean_file} -m Q.yeast+R7 -bb 1000 -bnni -nt ${NUM_CORES} -o $outgroup_name  #full version, force Q.yeast+R7 model because Modelfinder keeps deciding on this one

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
