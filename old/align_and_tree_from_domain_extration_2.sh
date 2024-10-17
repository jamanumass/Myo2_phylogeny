#!/bin/bash
#SBATCH --partition=gpu-long         # Select partition
#SBATCH --gpus=1                     # Request 1 GPU
#SBATCH --cpus-per-gpu=32            # 32 CPU cores per GPU
#SBATCH --mem-per-gpu=20G            # 20 GB memory per GPU
#SBATCH -t 48:00:00
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err

# Source the variables from the output_var_file
if [ -z "$output_var_file" ]; then
    echo "Error: output_var_file is not set."
    exit 1
fi

source $output_var_file

echo "Using the following settings from the main script:"
echo "Input FASTA: $INPUT_FASTA"
echo "Number of results to collect: $number_results_to_collect"
echo "Run directory: $run_dir"

# Change to the run directory
cd $run_dir || exit

# Define the alignment output file and perform alignment using 32 CPUs for MAFFT
aligned_file_name="${run_name}_hits_aligned.fasta"
clean_name_hits="${run_name}_clean_hits.fasta"

echo "Running MAFFT alignment..."
mafft --maxiterate 1000 --localpair --reorder --thread 32 "${clean_name_hits}" > "${aligned_file_name}"

# Infer the phylogenetic tree using IQ-TREE with 32 threads
tree_dir="${run_dir}/tree_results"
mkdir -p ${tree_dir}

echo "Running IQ-TREE for phylogenetic inference..."
iqtree2 -s "${tree_dir}/${aligned_file_name}" -m Q.yeast+R7 -bb 1000 -bnni -nt 32 -o $outgroup_name

# Cleanup unnecessary IQ-TREE log files
cd ${tree_dir}
echo "Cleaning up unnecessary IQ-TREE files..."
rm *splits.nex
rm *ckp.gz
rm *bionj
rm *mldist
rm *model.gz
rm *uniqueseq.phy