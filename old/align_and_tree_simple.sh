#!/bin/bash
#SBATCH --partition=cpu             # Partition to use
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks=1                  # Number of tasks (processes)
#SBATCH --cpus-per-task=17           # Number of CPU cores per task
#SBATCH --mem=120G                   # Memory per node
#SBATCH -t 10:00:00
#SBATCH -o /home/jaman_umass_edu/jaman/Myo2_phylogeny/other/deduplicate/tree_dedupe/slurm-%j.out
#SBATCH -e /home/jaman_umass_edu/jaman/Myo2_phylogeny/other/deduplicate/tree_dedupe/slurm-%j.err

#load required modules
#module load iq-tree/2.1.3-noMPI 
module load MAFFT/7.505-GCC-11.3.0-with-extensions

#define the run directory, will contain both the input and outputs
run_dir="/home/jaman_umass_edu/jaman/Myo2_phylogeny/other/deduplicate/tree_dedupe"


# Set the directory for alignment and tree results

#define the input 
input_seqs_file="${run_dir}/Myo2_Pedro_1_deduplicated.fa"

#Define clean_name input file
clean_input_name_file="${input_seqs_file}_cleanNames.fa" 
rm "${clean_input_name_file}"
touch "${clean_input_name_file}"

#add outgroup
#Alternate outgroup: use a Myo IV gene from Bigelowiella natans
outgroup_name="Bnat_34821_MyoIV"
echo ">${outgroup_name}" >> "${clean_input_name_file}"
echo "MQRRFEQNKIYTNVGTILISVNPYQRLPLYTEQVLKKYTSRGLGVVDMPPHVFNIAHDAFYGVTSFSKGQSIIISGESGAGKTEATKQCLQYLAAIAGSTSDVEKKVLRANPILEAFGNAKTLRNDNSSRFGKYLEIYLDEKGRISSSATENYLLEKIRVVQPSLKERNFHIFYQLVKAASSKLRQKLKLKEDAGKYNYLKSCTDVPSIDDTRDYKEVIEAFRELGISDSEREQTFRICAAILHLGNCTFTDGKNHSGGCQVNEKAVLRDAAELLGVNGDKLLERLTTREIRVRGQSAAKAVMGAEEASDTRHALCKFVYGRMFDWIVARINKSMPGGGGGRSIGILDIFGFEIFEKNSFEQLCINFTNERLQQHFNRHTFKLEQNIYCSEGIDFDEIDYIDNQPMVDLITKKPHGVLPLLDEELRIPKGSDETFLAKLETKQCKNPVFKRQMKKRAHFAIKHYAGKVLYHCKGFLEKNRDTLTEDLVEILQTSNQPLLQELYPSDMQISSKQRKSSLATQFQQQLTRLMHSLNATQPHYIRCIKPNNDKAPMKFVAKNCHEQLLYSGVFEAVAIRKQGFPFRLSHEEFEKRYSICLGKTSALSNQSSVKGRCKLILNEMKCDPKNTRIGSSRVLYR" >> "$clean_input_name_file" 


#process the names in the input file to prevent downstream problems.
awk '
    # If the line starts with ">", it is a name line
    /^>/ {
        # Modify the name by replacing spaces with underscores, removing square brackets, commas, and "putative_"
        name = $0
        gsub(/ /, "_", name)                 # Replace spaces with underscores
        gsub(/\[|\]/, "", name)              # Remove square brackets
        gsub(/,/, "", name)                  # Remove commas
        gsub(/putative_/, "", name)          # Remove the word "putative_"
        print name
        next
    }
    # If the line does not start with ">", it is a sequence line
    {
        print
    }
' $input_seqs_file >> ${clean_input_name_file}




#Define the alignment output file and perform alignment
aligned_file_name="${clean_input_name_file}_aligned.fasta"
# mafft ${clean_input_name_file} > ${aligned_file_name} #fast version
mafft ${clean_input_name_file} > ${aligned_file_name}


#infer tree
	iqtree2 -s ${aligned_file_name} -m Q.yeast+R7 -bb 1000 -bnni -nt AUTO  #full version, force Q.yeast+R7 model because Modelfinder keeps deciding on this one
	
#remove iqtree log files

echo "cleaning up.."
	rm *splits.nex
	rm *ckp.gz
	rm *bionj
	rm *mldist
#	rm *.log
	rm *.iqtree
	rm *model.gz
	rm *uniqueseq.phy
