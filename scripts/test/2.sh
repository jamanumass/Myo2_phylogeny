#!/bin/bash
#SBATCH --partition=cpu-long        # Partition to use
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks=1                  # Number of tasks (processes)
#SBATCH --cpus-per-task=1           # Number of CPU cores per task
#SBATCH --mem=2G                    # Memory per node
#SBATCH -t 00:10:00                 # Time limit (hh:mm:ss)
#SBATCH -o /work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny/scripts/test/slurm-2-%j.out
#SBATCH -e /work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny/scripts/test/slurm-2-%j.err

# Load required modules (if any)
# module load iq-tree/2.1.3-noMPI 
# module load mafft/7.481

# Check if the 1_output.fasta file exists
if [ ! -f "1_output.fasta" ]; then
  echo "Error: 1_output.fasta does not exist. Run 1.sh first."
  exit 1
fi

# Process the file: remove the first line
tail -n +2 1_output.fasta > 2_output.fasta

echo "2.sh: Processed 1_output.fasta and created 2_output.fasta"