#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=60G
#SBATCH -t 01:00:00
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err

# Load required modules
module load uri/main
module load HMMER/3.3.2-iimpi-2021b
module load MAFFT/7.505-GCC-11.3.0-with-extensions

# ========= User Configurable Parameters =========
INPUT_FASTA="/home/jaman_umass_edu/jaman/Myo2_phylogeny/input_seqs/Myo2_Pedro_input.fasta"           # Fasta file with known gene group sequences
HMM_THRESHOLD=1000                            # HMM gathering threshold for filtering hits

# Define database directories
DB_DIRS=(
    "/home/jaman_umass_edu/jaman/sequence_databases/proteins/Hets_expanded/transdecoder_output/cleaned_names"
#    "/home/jaman_umass_edu/jaman/sequence_databases/proteins/Eukprot"
#    "/home/jaman_umass_edu/jaman/sequence_databases/proteins/Uniprot"
)
OUTPUT_FASTA="$(basename "${INPUT_FASTA%.fasta}_high_confidence_genes.fasta")"  # Final output: high-confidence genes (fasta) in working directory
NUM_CORES=4                                                      # Number of CPU cores for MAFFT/hmmsearch

# ========= Step 1: Align Input Sequences =========
# Skipping alignment step as input is already aligned
ALIGNED_FASTA="${INPUT_FASTA}"

# ========= Step 2: Build HMM Profile =========
HMM_PROFILE="${INPUT_FASTA%.fasta}.hmm" # Define the HMM profile output file name based on the input fasta file name
echo "Building HMM profile..."
hmmbuild "${HMM_PROFILE}" "${ALIGNED_FASTA}"

# ========= Step 3: Initialize Output File =========
: > "${OUTPUT_FASTA}"   # Truncate or create the output file

# ========= Step 4: Scan Databases and Collect Hits =========
mkdir -p scan_results  # Create the scan_results directory if it doesn't exist

for db_dir in "${DB_DIRS[@]}"; do
    find "${db_dir}" -name "*.fa" -type f | while read db_file; do
        echo "Scanning database: ${db_file}"
        
        TBL_FILE="scan_results/$(basename "${db_file%.fa}")_$(basename "${HMM_PROFILE}" .hmm)_tbl.txt"
        hmmsearch --tblout "${TBL_FILE}" --cpu "${NUM_CORES}" "${HMM_PROFILE}" "${db_file}" > /dev/null 2>&1
        
        # Filter hits above the threshold and extract their IDs
        grep -v '^#' "${TBL_FILE}" | awk -v thresh="${HMM_THRESHOLD}" '$6 >= thresh {print $1}' | while read hit_id; do
            awk -v id="$hit_id" '
            BEGIN {found=0}
            /^>/ {
                header = $0
                sub(/^>/, "", header)
                split(header, tokens, " ")
                if (tokens[1] == id) { found=1; print; next }
            }
            found && /^>/ { found=0; exit }
            found { print }
            ' "${db_file}" >> "${OUTPUT_FASTA}"
        done
        done
    done
done
