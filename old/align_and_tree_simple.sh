#!/bin/bash
#SBATCH --partition=cpu             # Partition to use
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks=1                  # Number of tasks (processes)
#SBATCH --cpus-per-task=17           # Number of CPU cores per task
#SBATCH --mem=120G                   # Memory per node
#SBATCH -t 2:00:00
#SBATCH --output=slurm-%j.out    # Output log file (%j will be replaced by the job ID)
#SBATCH --error=slurm-%j.err     # Error log file

# Load required modules for MPI, IQ-TREE, and MAFFT
module load uri/main
module load mpich/4.2.1 iq-tree/2.3.1
module load MAFFT/7.505-GCC-11.3.0-with-extensions
module load gcc/9.4.0


##define the run directory, will contain both the input and outputs
#run_dir="/home/jaman_umass_edu/jaman/Myo2_phylogeny/results/sup_figure_run"
#
#
## Set the directory for alignment and tree results
#
##define the input 
#input_seqs_file="${run_dir}/sup_fig_run_combined.fasta"
#clean_input_name_file="${input_seqs_file}_cleaned.fasta"
#
##add outgroup
##Alternate outgroup: use a Myo IV gene from Bigelowiella natans
outgroup_name="Bnat_34821_MyoIV"

#clean the names in the input file to prevent downstream problems.
#/home/jaman_umass_edu/jaman/general_tools/clean_gene_names.sh $input_seqs_file >> ${clean_input_name_file}


##Define the alignment output file and perform alignment
#aligned_file_name="${clean_input_name_file}_aligned.fasta"
## mafft ${clean_input_name_file} > ${aligned_file_name} #fast version
#mafft ${clean_input_name_file} > ${aligned_file_name}
#
## ==============================
## Tree Inference with IQ-TREE
## ==============================
#


## Update the final_data variable to point to the file in the new directory

final_data=sup_run_input.fasta

# Run IQ-TREE with best model search, additional support tests, and outgroup assignment
iqtree2 -s "${final_data}" -m Q.yeast+R7 -bb 1000 -nt 16 -o "$outgroup_name" --keep-ident

echo "Finished IQ-TREE inference"

