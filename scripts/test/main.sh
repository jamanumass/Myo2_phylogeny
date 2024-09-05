#!/bin/bash

# Check if an argument was passed
if [ -z "$1" ]; then
  echo "Usage: $0 <input_fasta_file>"
  exit 1
fi

# Set the environment variable
export INPUT_FASTA=$1

# Run the first script
job1=$(sbatch --export=ALL /work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny/scripts/test/1.sh)
job1_id=$(echo $job1 | awk '{print $4}')

# Run the second script after the first one completes
sbatch --dependency=afterok:$job1_id --export=ALL /work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny/scripts/test/2.sh