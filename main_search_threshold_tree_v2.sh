#!/bin/bash

# =======================================================
# CONFIGURATION: Step Flags to Control Workflow
# Set to "true" or "false" to enable or disable each step
# =======================================================
run_step_1=false    # Step 1: HMMER Search
run_step_2=false   # Step 2: Domain Extraction
run_step_3=true    # Step 3: Tree Inference

# =======================================================
# GENERAL PARAMETERS AND PATHS
# =======================================================
NUM_CORES=48 #Set to the number of CPUs to request. Downstream programs are set to use this number minus 1 threads. Note that systematic testing resulted in optimisation at 48 cores on Unity cluster
hmm_score_threshold=1400
num_top_hmm_hits_per_genome=3
use_hmm_score_threshold="true"
INPUT_FASTA_FOR_HMMSEARCH="/work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny/results/sup_figure_run_8/sup_run_8_input.fasta"

# Specify additional FASTA files as a space-separated string
#ADDITIONAL_FASTAS_TO_RUN_IN_TREE="/home/jaman_umass_edu/jaman/Myo2_phylogeny/results/results_Myo2_Pedro_9/myo2hits_vahl.fasta"

# Specify the previous run directory if reusing data from a previous run
# Leave as an empty string "" to create a new run
PREVIOUS_RUN_DIR="/work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny/results/sup_figure_run_8"
#PREVIOUS_RUN_DIR=""
# =======================================================
# DATABASE SELECTION: Define Directories for Protein Databases
# =======================================================
PROTEIN_DB_OPTION="all"
db1="/home/jaman_umass_edu/jaman/sequence_databases/proteins/Uniprot"
db2="/home/jaman_umass_edu/jaman/sequence_databases/proteins/Eukprot"
database_dirs=()
[ "$PROTEIN_DB_OPTION" = "1" ] && database_dirs=("$db1")
[ "$PROTEIN_DB_OPTION" = "2" ] && database_dirs=("$db2")
[ "$PROTEIN_DB_OPTION" = "all" ] && database_dirs=("$db1" "$db2")

# =======================================================
# SETUP: Define Run Directory
# =======================================================
root_dir="/work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny"
results_dir="${root_dir}/results"


if [ -n "$PREVIOUS_RUN_DIR" ]; then
    if [ ! -d "$PREVIOUS_RUN_DIR" ]; then
        echo "Error: Specified previous run directory '${PREVIOUS_RUN_DIR}' does not exist."
        exit 1
    fi
    echo "Restarting using data from previous run in directory: ${PREVIOUS_RUN_DIR}"
    run_dir="$PREVIOUS_RUN_DIR"
    run_name=$(basename "$run_dir")
else
    # Generate a unique run name and directory
    run_name=$(basename $INPUT_FASTA_FOR_HMMSEARCH | sed 's/_input\.fasta$//')
    run_dir="${results_dir}/results_${run_name}"
    while [ -d "$run_dir" ]; do
        if [[ "$run_name" =~ _[0-9]+$ ]]; then
            suffix=$(echo "$run_name" | sed 's/.*_\([0-9]\+\)$/\1/')
            suffix=$((suffix + 1))
            run_name=$(echo "$run_name" | sed 's/_[0-9]\+$//')
        else
            suffix=1
        fi
        run_name="${run_name}_${suffix}"
        run_dir="${results_dir}/results_${run_name}"
    done
    echo "Initiating new run directory named ${run_name} in the results folder."
    mkdir -p "$run_dir"
fi

# =======================================================
# SHARED VARIABLES FILE
# =======================================================
shared_variables_file="${run_dir}/shared_variables.txt"
if [ ! -f "${shared_variables_file}" ]; then
    {
        echo "INPUT_FASTA_FOR_HMMSEARCH=${INPUT_FASTA_FOR_HMMSEARCH}"
        echo "ADDITIONAL_FASTAS_TO_RUN_IN_TREE=\"${ADDITIONAL_FASTAS_TO_RUN_IN_TREE}\""
        echo "num_top_hmm_hits_per_genome=${num_top_hmm_hits_per_genome}"
        echo "run_name=${run_name}"
        echo "run_dir=${run_dir}"
        echo "results_dir=${results_dir}"
        echo "database_dirs=\"${database_dirs[@]}\""
        echo "NUM_CORES=${NUM_CORES}"
        echo "hmm_score_threshold=${hmm_score_threshold}"
        echo "use_hmm_score_threshold=${use_hmm_score_threshold}"
    } > "${shared_variables_file}"
    echo "Created new shared variables file: $(basename $shared_variables_file) in $run_name directory"
else
    echo "Using existing shared variables file: ${shared_variables_file}"
fi

# =======================================================
# DEFINE OUTPUT DIRECTORIES FOR EACH STEP (Conditional)
# =======================================================
if [ "$run_step_1" = true ]; then
    hmm_search_results_dir="${run_dir}/hmm_search_results"
    mkdir -p "$hmm_search_results_dir"
    echo "hmm_search_results_dir=${hmm_search_results_dir}" >> "${shared_variables_file}"
else
    echo "HMMER Search is set to be skipped. Directory not created."
fi

if [ "$run_step_2" = true ]; then
    domain_extraction_dir="${run_dir}/domain_extraction_result"
    mkdir -p "$domain_extraction_dir"
    echo "domain_extraction_dir=${domain_extraction_dir}" >> "${shared_variables_file}"
else
    echo "Domain Extraction is set to be skipped. Directory not created."
fi

if [ "$run_step_3" = true ]; then
    tree_inference_dir="${run_dir}/tree_inference_result"
    mkdir -p "$tree_inference_dir"
    echo "tree_inference_dir=${tree_inference_dir}" >> "${shared_variables_file}"
else
    echo "Tree Inference is set to be skipped. Directory not created."
fi

# Define paths to scripts
HMMER_script="/work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny/scripts/HMMER_search_v3.sh"
extract_single_domain_script="/work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny/scripts/extract_motor_domain_v2.sh"
tree_inference_script="/work/pi_lfritzlaylin_umass_edu/users/jaman/Myo2_phylogeny/scripts/tree_infer_4.sh"

# =======================================================
# STEP 1: HMMER Search
# =======================================================
if [ "$run_step_1" = true ]; then
    job1_out="${hmm_search_results_dir}/slurm-${run_name}_hmm_search_%j.out"
    job1_err="${hmm_search_results_dir}/slurm-${run_name}_hmm_search_%j.err"
    job1=$(sbatch --cpus-per-task=${NUM_CORES} --export=ALL,shared_variables_file="${shared_variables_file}" -o "$job1_out" -e "$job1_err" "$HMMER_script")    job1_id=$(echo "$job1" | awk '{print $4}')
    echo "Step 1: HMMER Search job submitted with Job ID: $job1_id"
else
    echo "Step 1: HMMER Search is skipped."
fi

# =======================================================
# STEP 2: Domain Extraction
# =======================================================
if [ "$run_step_2" = true ]; then
    job2_out="${domain_extraction_dir}/slurm-${run_name}_extract_domain_%j.out"
    job2_err="${domain_extraction_dir}/slurm-${run_name}_extract_domain_%j.err"
    job2=$(sbatch --dependency=afterok:$job1_id --cpus-per-task=${NUM_CORES} --export=ALL,shared_variables_file="${shared_variables_file}" -o "$job2_out" -e "$job2_err" "$extract_single_domain_script")    job2_id=$(echo "$job2" | awk '{print $4}')
    echo "Step 2: Domain Extraction job submitted with Job ID: $job2_id"
else
    echo "Step 2: Domain Extraction is skipped."
fi

# =======================================================
# STEP 3: Tree Inference
# =======================================================
if [ "$run_step_3" = true ]; then
    job3_out="${tree_inference_dir}/slurm-${run_name}_infer_tree_%j.out"
    job3_err="${tree_inference_dir}/slurm-${run_name}_infer_tree_%j.err"
    
    # Set dependency for Tree Inference based on whether Step 2 or Step 1 ran
    if [ "$run_step_2" = true ]; then
        dependency="--dependency=afterok:$job2_id"
    elif [ "$run_step_1" = true ]; then
        dependency="--dependency=afterok:$job1_id"
    else
        dependency=""
    fi

# Submit the job with or without dependencies
if [ -n "$dependency" ]; then
    job3=$(sbatch $dependency --cpus-per-task=${NUM_CORES} --export=ALL,shared_variables_file="${shared_variables_file}" -o "$job3_out" -e "$job3_err" "$tree_inference_script")
else
    job3=$(sbatch --cpus-per-task=${NUM_CORES} --export=ALL,shared_variables_file="${shared_variables_file}" -o "$job3_out" -e "$job3_err" "$tree_inference_script")
fi

    job3_id=$(echo "$job3" | awk '{print $4}')
    echo "Step 3: Tree Inference job submitted with Job ID: $job3_id"
else
    echo "Step 3: Tree Inference is skipped."
fi