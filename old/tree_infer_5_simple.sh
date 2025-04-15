#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=60G
#SBATCH --time=05:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

# Load required modules
module load uri/main
module load mpich/4.2.1 iq-tree/2.3.1
module load MAFFT/7.505-GCC-11.3.0-with-extensions
module load gcc/9.4.0




input_seqs="/home/jaman_umass_edu/jaman/Myo2_phylogeny/results/sup_figure_run_17/sup_run_17_input.fasta"
input_aligned="${input_seqs%.fasta}_aligned.fasta"

touch $input_aligned

mafft --quiet --localpair --thread -1 "$input_seqs" > "$input_aligned"

# Run IQ-TREE on the final input file
iqtree2 -s "$input_aligned" -m Q.yeast+R7 -bb 5000 -nt 47


