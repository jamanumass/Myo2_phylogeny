#!/bin/bash

# Set the parameters
hmm_score_threshold=825  # Set the HMM score threshold 
number_results_to_collect=10  # Alternatively, define a set number of results from each protein database to collect.
use_hmm_score_threshold="false"


INPUT_FASTA="/work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny/input_seqs/pedros_1F_input.fasta"  # Set your input FASTA file here
#NUM_CORES=20  # Set this variable as needed to dynamically control the number of cores requested
NUM_CORES=16  # Set this variable as needed to dynamically control the number of cores requested

# Location of permanent protein databases
protein_database_directory="/work/pi_lfritzlaylin_umass_edu/users/jaman/sequence_databases/myo2_seqs_collect_v1/full_genomes_peptide_databases"

# Location of project directories
root_dir="/work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny"
results_dir="${root_dir}/results"

# Check if the input FASTA file exists
if [ ! -f "$INPUT_FASTA" ]; then
    echo "Error: Input FASTA file '$INPUT_FASTA' not found."
    exit 1
fi

# Set up the base run name and directory
run_name=$(basename $INPUT_FASTA | sed 's/_input\.fasta$//')
run_dir="${results_dir}/results_${run_name}"

# Check if run_dir exists and adjust the run_name if needed
while [ -d "$run_dir" ]; do
    # Check if run_name ends with _[0-9]
    if [[ "$run_name" =~ _[0-9]+$ ]]; then
        # Extract the numeric suffix, increment it, and replace in run_name
        suffix=$(echo "$run_name" | sed 's/.*_\([0-9]\+\)$/\1/')  # Extract the number
        suffix=$((suffix + 1))  # Increment the suffix
        run_name=$(echo "$run_name" | sed 's/_[0-9]\+$//')  # Remove the old suffix
    else
        # Append _1 if no numeric suffix exists
        suffix=1
    fi
    run_name="${run_name}_${suffix}"  # Update run_name with new suffix
    run_dir="${results_dir}/results_${run_name}"  # Update run_dir with new run_name
done

# Create the unique run directory
mkdir -p "$run_dir"

# File to store all the variables co-used by different sbatch jobs
output_var_file="${run_dir}/job_variables.txt"

# Store the variables in the output variable file to be accessed by any sub-scripts
echo "INPUT_FASTA=${INPUT_FASTA}" >> ${output_var_file}
echo "number_results_to_collect=${number_results_to_collect}" >> ${output_var_file}
echo "run_name=${run_name}" >> ${output_var_file}
echo "run_dir=${run_dir}" >> ${output_var_file}
echo "results_dir=${results_dir}" >> ${output_var_file}
echo "protein_database_directory=${protein_database_directory}" >> ${output_var_file}
echo "NUM_CORES=${NUM_CORES}" >> ${output_var_file}
echo "hmm_score_threshold=${hmm_score_threshold}" >> ${output_var_file}
echo "use_hmm_score_threshold=${use_hmm_score_threshold}" >> "${output_var_file}"


# Explicitly export the output_var_file variable for the downstream scripts
export output_var_file

#### DOWNSTREAM SCRIPT 1 - Search ####

# Submit the HMMER_search_v2.sh script
HMMER_script="/work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny/scripts/HMMER_search_v3.sh"

# Define the hmmsearch directory
hmm_search_results_dir="${run_dir}/hmm_search_results"
mkdir -p ${hmm_search_results_dir}
echo "hmm_search_results_dir=${hmm_search_results_dir}" >> ${output_var_file} # Store this directory in the common variable file

# Set the log files to the results directory for this script
job1_out="${hmm_search_results_dir}/slurm-${run_name}_%j.out"
job1_err="${hmm_search_results_dir}/slurm-${run_name}_%j.err"

# Run the script, and store the job details
job1=$(sbatch --cpus-per-task=$((NUM_CORES + 1)) --export=output_var_file -o $job1_out -e $job1_err $HMMER_script)
job1_id=$(echo $job1 | awk '{print $4}')

# Check if the job submission succeeded
if [ -z "$job1_id" ]; then
    echo "Error: HMMER search job submission failed."
    exit 1
fi

# Report
echo "Search job submitted with Job ID: $job1_id"




#### DOWNSTREAM SCRIPT 2 - Domain Extraction ####

# Submit the extract_motor_domain_v2.sh script only after the first job (job1) completes
extract_single_domain_script="/work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny/scripts/extract_motor_domain_v2.sh"  # Path to extract_motor_domain_v2.sh

# Define the domain extraction directory
domain_extraction_dir="${run_dir}/domain_extraction_result"
mkdir -p $domain_extraction_dir
echo "domain_extraction_dir=${domain_extraction_dir}" >> ${output_var_file} # Store this directory in the common variable file

# Set the log files to the results directory for this script
job2_out="${domain_extraction_dir}/slurm-${run_name}_extract_motor_domain_%j.out"
job2_err="${domain_extraction_dir}/slurm-${run_name}_extract_motor_domain_%j.err"

# Run the script, and store the job details
job2=$(sbatch --cpus-per-task=$((NUM_CORES + 1)) --dependency=afterok:$job1_id --export=output_var_file -o $job2_out -e $job2_err ${extract_single_domain_script})
job2_id=$(echo $job2 | awk '{print $4}')

# Check if the job submission succeeded
if [ -z "$job2_id" ]; then
    echo "Error: Extract motor domain job submission failed."
    exit 1
fi

# Report
echo "Extract motor domain job submitted with Job ID: $job2_id"




#### DOWNSTREAM SCRIPT 3 - Tree Inference ####

# Submit the tree inference script only after job2 completes
tree_inference_script="/work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny/scripts/tree_infer_2.sh"  # Path to tree inference script

# Define the directory for tree inference results
tree_inference_dir="${run_dir}/tree_inference_result"
mkdir -p $tree_inference_dir
echo "tree_inference_dir=${tree_inference_dir}" >> ${output_var_file} # Store this directory in the common variable file

# Set the log files to the results directory for this script
job3_out="${tree_inference_dir}/slurm-${run_name}_infer_tree_%j.out"
job3_err="${tree_inference_dir}/slurm-${run_name}_infer_tree_%j.err"

# Run the script, and store the job details
job3=$(sbatch --cpus-per-task=$((NUM_CORES + 1)) --dependency=afterok:$job2_id --export=output_var_file -o $job3_out -e $job3_err ${tree_inference_script})
#job3=$(sbatch --partition=gpu --gpus=1 --cpus-per-task=$((NUM_CORES + 1)) --dependency=afterok:$job2_id --export=output_var_file -o $job3_out -e $job3_err ${tree_inference_script})
job3_id=$(echo $job3 | awk '{print $4}')

# Check if the job submission succeeded
if [ -z "$job3_id" ]; then
    echo "Error: Tree inference job submission failed."
    exit 1
fi

# Report
echo "Tree inference job submitted with Job ID: $job3_id"
