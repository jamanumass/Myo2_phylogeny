#!/bin/bash
#SBATCH --partition=cpu-long        # Partition to use
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks=1                  # Number of tasks (processes)
#SBATCH --cpus-per-task=32           # Number of CPU cores per task
#SBATCH --mem=120G                   # Memory per node
SBATCH -t 4:00:00
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err


#This script takes the extracted domains from the HMM search hits, and build a phylogenetic tree


#load required modules
module load iq-tree/2.1.3-noMPI 
module load mafft/7.481


# Set the root directory for Myo searches
root_dir="/work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny"
cd $root_dir

# define input and run names
input_file="MyoI_motors_input.fasta"
input_path="/work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny/input_seqs"
input_name="${input_file%_input.fasta}"
run_name=$input_name

#Location of permanent protein databases
protein_database_directory="/work/pi_lfritzlaylin_umass_edu/users/jaman/sequence_databases/myo2_seqs_collect_v1/full_genomes_peptide_databases"

#define the run directory 
run_dir="/work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny/results_${input_name}"

# Set the directory for alignment and tree results
tree_dir="${run_dir}/tree_results"
mkdir -p $tree_dir

#define the $hits_domain_extraction_file
extraction_dir="${run_dir}/domain_extraction_result"
hits_domain_extraction_file="${extraction_dir}/${input_name}_hits_domain_extraction.fa"

#Define clean_name output file
clean_name_hits="${tree_dir}/${input_name}_hits_domainExtraction_cleanNames.fa" #

#populate this with the list of extracted domains
cat $hits_domain_extraction_file > $clean_name_hits.temp #copies original to a temp file

#add outgroup


#Use Dicty Myo A1  XP_641363.1 (no relation to the -A1 option in the grep function). grep first 9 lines, roughly corresponds to the motor domain
#delete if already present, if not pull the sequence from database and add to the list
#if ! grep -q "XP_641363.1" "$clean_name_hits.temp"; then
#    grep -A9 "XP_641363.1" "/work/pi_lfritzlaylin_umass_edu/users/jaman/sequence_databases/myo2_seqs_collect_v1/full_genomes_peptide_databases/GCF_000004695.1_dicty_2.7_protein.faa" >> "$clean_name_hits.temp"
#fi
#Alternate outgroup: use a Myo IV gene from Bigelowiella natans
outgroup_name="Bnat_34821_MyoIV"
echo ">${outgroup_name}" >> "${clean_name_hits}.temp"
echo "MQRRFEQNKIYTNVGTILISVNPYQRLPLYTEQVLKKYTSRGLGVVDMPPHVFNIAHDAFYGVTSFSKGQSIIISGESGAGKTEATKQCLQYLAAIAGSTSDVEKKVLRANPILEAFGNAKTLRNDNSSRFGKYLEIYLDEKGRISSSATENYLLEKIRVVQPSLKERNFHIFYQLVKAASSKLRQKLKLKEDAGKYNYLKSCTDVPSIDDTRDYKEVIEAFRELGISDSEREQTFRICAAILHLGNCTFTDGKNHSGGCQVNEKAVLRDAAELLGVNGDKLLERLTTREIRVRGQSAAKAVMGAEEASDTRHALCKFVYGRMFDWIVARINKSMPGGGGGRSIGILDIFGFEIFEKNSFEQLCINFTNERLQQHFNRHTFKLEQNIYCSEGIDFDEIDYIDNQPMVDLITKKPHGVLPLLDEELRIPKGSDETFLAKLETKQCKNPVFKRQMKKRAHFAIKHYAGKVLYHCKGFLEKNRDTLTEDLVEILQTSNQPLLQELYPSDMQISSKQRKSSLATQFQQQLTRLMHSLNATQPHYIRCIKPNNDKAPMKFVAKNCHEQLLYSGVFEAVAIRKQGFPFRLSHEEFEKRYSICLGKTSALSNQSSVKGRCKLILNEMKCDPKNTRIGSSRVLYR" >> "${clean_name_hits}.temp" 


#process the names in the temp file to prevent downstream problems.
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
' "$clean_name_hits.temp" > "${clean_name_hits}"

#remove the temp file
rm $clean_name_hits.temp



#Define the alignment output file and perform alignment
aligned_file_name="${input_name}_hits_aligned.fasta"
# mafft "${clean_name_hits}" > "${tree_dir}/${aligned_file_name}" #fast version
mafft --maxiterate 1000 --localpair --reorder "${clean_name_hits}" > "${tree_dir}/${aligned_file_name}" #good version




#infer tree
#	iqtree2 -s ${tree_dir}/${aligned_file_name} -fast -m LG+F+R7 -B 1000 --wbtl -nstop 500 #quick version
#	iqtree2 -s ${tree_dir}/${aligned_file_name} -m MFP -bb 1000 -bnni -o $outgroup_name #full version, infer site rates
	iqtree2 -s ${tree_dir}/${aligned_file_name} -m Q.yeast+R7 -bb 1000 -bnni -nt AUTO -o $outgroup_name #full version, force Q.yeast+R7 model because Modelfinder keeps deciding on this one
	
#remove iqtree log files
cd ${tree_dir} when 
echo "cleaning up.."
	rm *splits.nex
	rm *ckp.gz
	rm *bionj
	rm *mldist
#	rm *.log
#	rm *.iqtree
	rm *model.gz
	rm *uniqueseq.phy
