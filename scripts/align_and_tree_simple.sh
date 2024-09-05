#!/bin/bash
#SBATCH --partition=cpu-long        # Partition to use
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks=1                  # Number of tasks (processes)
#SBATCH --cpus-per-task=128           # Number of CPU cores per task
#SBATCH --mem=120G                   # Memory per node
SBATCH -t 12:00:00
#SBATCH -o /home/jaman_umass_edu/jaman/Myo2_phylogeny/results_Amoebozoan_MyoII_motors/refine1/slurm-%j.out
#SBATCH -e /home/jaman_umass_edu/jaman/Myo2_phylogeny/results_Amoebozoan_MyoII_motors/refine1/slurm-%j.err

#load required modules
module load iq-tree/2.1.3-noMPI 
module load mafft/7.481

#define the run directory, will contain both the input and outputs
run_dir="/home/jaman_umass_edu/jaman/Myo2_phylogeny/results_Amoebozoan_MyoII_motors/refine1"


# Set the directory for alignment and tree results
tree_dir="${run_dir}/tree_results"
mkdir -p $tree_dir

#define the input 
input_seqs_file="refined1_all_seqs_input.fasta"

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
mafft --maxiterate 1000 --localpair --reorder ${clean_input_name_file} > ${aligned_file_name} #good version

#copy the alignment to the tree dir and run IQtree there.
cp ${aligned_file_name} ${tree_dir}/${aligned_file_name}
cd ${tree_dir}

#infer tree
#	iqtree2 -s ${tree_dir}/${aligned_file_name} -fast -m LG+F+R7 -B 1000 --wbtl -nstop 500 #quick version
#	iqtree2 -s ${tree_dir}/${aligned_file_name} -m MFP -bb 1000 -bnni -o $outgroup_name #full version, infer site rates
	iqtree2 -s ${tree_dir}/${aligned_file_name} -m Q.yeast+R7 -bb 1000 -bnni -nt AUTO  #full version, force Q.yeast+R7 model because Modelfinder keeps deciding on this one
	
#remove iqtree log files
cd ${tree_dir} 
echo "cleaning up.."
	rm *splits.nex
	rm *ckp.gz
	rm *bionj
	rm *mldist
#	rm *.log
#	rm *.iqtree
	rm *model.gz
	rm *uniqueseq.phy
