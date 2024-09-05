#!/bin/bash
#SBATCH --partition=cpu-long        # Partition to use
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks=1                  # Number of tasks (processes)
#SBATCH --cpus-per-task=1           # Number of CPU cores per task
#SBATCH --mem=2G                    # Memory per node
#SBATCH -t 00:10:00                 # Time limit (hh:mm:ss)
#SBATCH -o /work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny/scripts/test/slurm-1-%j.out
#SBATCH -e /work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny/scripts/test/slurm-1-%j.err

# Load required modules (if any)
# module load iq-tree/2.1.3-noMPI 
# module load mafft/7.481

# Check if INPUT_FASTA environment variable is set
if [ -z "$INPUT_FASTA" ]; then
  echo "Error: INPUT_FASTA environment variable is not set."
  exit 1
fi

# Process the file: keep only the first two lines
head -n 2 "$INPUT_FASTA" > 1_output.fasta

echo "1.sh: Processed $INPUT_FASTA and created 1_output.fasta"