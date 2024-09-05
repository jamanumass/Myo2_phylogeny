#!/bin/bash

# Set the number of results to collect and the input FASTA file
number_results_to_collect=10  # Define the number of results from each protein database to collect
INPUT_FASTA="/work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny/input_seqs/myo1s_2_repeat_input.fasta"  # Set your input FASTA file here
NUM_CORES=20  # Set this variable as needed to dynamically control the number of cores


# Check if the input FASTA file exists
if [ ! -f "$INPUT_FASTA" ]; then
    echo "Error: Input FASTA file '$INPUT_FASTA' not found."
    exit 1
fi

# Extract run name from the input FASTA file
run_name=$(basename $INPUT_FASTA | sed 's/_input\.fasta$//')

# Set up directories
root_dir="/work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny"
results_dir="${root_dir}/results"
run_dir="${results_dir}/results_${run_name}"
mkdir -p ${run_dir}

# Location of permanent protein databases
protein_database_directory="/work/pi_lfritzlaylin_umass_edu/users/jaman/sequence_databases/myo2_seqs_collect_v1/full_genomes_peptide_databases"

# Output file for timing info
timing_log="${run_dir}/timing.log" # Output file for timing info
    echo "Run times for ${run_name} components" >> ${timing_log}
output_var_file="${run_dir}/job_variables.txt" # File to store all the variables co-used by different sbatch jobs

# Store the variables in a the output variable file to be accessed by any sub-scripts
echo "timing_log=${timing_log}" > ${output_var_file}
echo "INPUT_FASTA=${INPUT_FASTA}" > ${output_var_file}
echo "number_results_to_collect=${number_results_to_collect}" >> ${output_var_file}
echo "run_name=${run_name}" >> ${output_var_file}
echo "run_dir=${run_dir}" >> ${output_var_file}
echo "results_dir=${results_dir}" >> ${output_var_file}
echo "protein_database_directory=${protein_database_directory}" >> ${output_var_file}
echo "NUM_CORES=${NUM_CORES}" >> ${output_var_file}

# Explicitly export the output_var_file variable for the downstream scripts
export output_var_file



####DOWNSTREAM SCRIPTS####

# Submit the HMMER_search_v2.sh script 
HMMER_script="/work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny/scripts/HMMER_search_v2.sh"

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

# Report
echo "HMMER job submitted with Job ID: $job1_id"
#echo "All files will be in it's results directory ${hmm_search_results_dir}"



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

# Report
echo "Extract motor domain job submitted with Job ID: $job2_id"
#echo "All files will be in it's results directory ${domain_extraction_dir}"





# Submit the tree inference script only after job2 completes
tree_inference_script="/work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny/scripts/tree_infer_2.sh"  # Path to tree inference script

# Define the directory
tree_inference_dir="${run_dir}/tree_inference_result"
mkdir -p $tree_inference_dir
echo "tree_inference_dir=${tree_inference_dir}" >> ${output_var_file} # Store this directory in the common variable file

# Set the log files to the results directory for this script
job3_out="${tree_inference_dir}/slurm-${run_name}_extract_motor_domain_%j.out"
job3_err="${tree_inference_dir}/slurm-${run_name}_extract_motor_domain_%j.err"

# Run the script, and store the job details
job3=$(sbatch --cpus-per-task=$((NUM_CORES + 1)) --dependency=afterok:$job2_id --export=output_var_file -o $job3_out -e $job3_err ${tree_inference_script})
#job3=$(sbatch --dependency=afterok:$job2_id --export=output_var_file -o $job3_out -e $job3_err ${tree_inference_script})
job3_id=$(echo $job3 | awk '{print $4}')

# Report
echo "Tree inference job submitted with Job ID: $job3_id"
#echo "All files will be in it's results directory ${tree_inference_dir}"


