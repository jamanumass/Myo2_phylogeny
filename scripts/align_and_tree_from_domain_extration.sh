#!/bin/bash
#SBATCH --partition=gpu,gpu-preempt,gpu-long
#SBATCH --constraint=vram8
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-gpu=8
#SBATCH --mem-per-gpu=20G
#SBATCH -t 8:00:00
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
input_file="Amoebozoan_MyoII_motors_input.fasta"
input_path="/work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny/input_seqs"
input_name="${input_file%_input.fasta}"
run_name=$input_name

#Location of permanent protein databases
protein_database_directory="/work/pi_lfritzlaylin_umass_edu/users/jaman/sequence_databases/myo2_seqs_collect_v1/full_genomes_peptide_databases"

#define the run directory 
run_dir="/work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny/results_${input_name}"

# Set the directory for alignment and tree results
tree_dir="${run_dir}/tree_results"
mkdir $tree_dir

#define the $hits_domain_extraction_file
extraction_dir="${run_dir}/domain_extraction_result"
hits_domain_extraction_file="${extraction_dir}/${input_name}_hits_domain_extraction.fa"


#Define the alignment output file and perform alignment
aligned_file_name="${input_name}_hits_aligned.fasta"
mafft ${hits_domain_extraction_file} > ${tree_dir}/${aligned_file_name}

#infer tree
echo "Running IQtree in fast mode with auto find model of evo and inference of site rate evolution.."
	iqtree2 -s ${tree_dir}/${aligned_file_name} -fast -m LG+F+R7 #quick version
#	iqtree -s ${tree_dir}/${aligned_file_name} -fast -nt 16 -m LG+F+R7 #old quick version
	#iqtree -s *out.fas -bb 1000 -wsr -nt 16 -m LG+F+R7 #full version, infer site rates
	
#remove iqtree log files
cd ${tree_dir}
echo "cleaning up.."
	rm *splits.nex
	rm *ckp.gz
	rm *bionj
	rm *mldist
	rm *.log
	rm *.iqtree
	rm *model.gz
	rm *uniqueseq.phy
