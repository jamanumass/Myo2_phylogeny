#!/bin/bash
#SBATCH --partition=gpu,gpu-preempt,gpu-long
#SBATCH --constraint=vram8
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-gpu=8
#SBATCH --mem-per-gpu=20G
#SBATCH -t 1:00:00
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err

module load iq-tree/2.1.3-noMPI 
module load mafft/7.481
module load uri/main
module load  HMMER/3.3.2-gompi-2022a

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <protein_input_file> <protein_database_file>"
    exit 1
fi

# Assign input arguments to variables
protein_input_file=$1
protein_database_file=$2

# Output file for best hits
output_file="hmmer_best_hits.txt"

# Run HMMER to find the best hits
hmmscan --tblout $output_file $protein_database_file $protein_input_file

# Check if hmmscan command was successful
if [ $? -ne 0 ]; then
    echo "HMMER search failed"
    exit 1
fi

# Parse the output to get the best hits
grep -v "^#" $output_file | sort -k5,5g | awk '!seen[$1]++' > best_hits.txt

# Display the best hits
cat best_hits.txt

echo "Best hits have been saved to best_hits.txt"