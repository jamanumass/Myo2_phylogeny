#!/bin/bash

# Script to filter sequences from a FASTA file based on sequence similarity
# either within species or globally.
#
# Usage: ./script.sh -i input_fasta -t threshold -o output_fasta -m mode [-s species_length]
#
# Arguments:
# -i input_fasta: The input FASTA file containing sequences
# -t threshold: The percent identity threshold for removing sequences (e.g., 90 for 90% identity)
# -o output_fasta: The output FASTA file where filtered sequences will be saved
# -m mode: 'species' for species-specific comparisons, 'global' for global comparisons
# -s species_length: (Optional) Number of characters from the beginning of each sequence header to use as the species identifier (only required in 'species' mode)
# -h: Show help and usage
#
# Description:
# In 'species' mode, the script compares sequences within the same species based
# on the species_length. If two sequences from the same species have a percent identity
# greater than or equal to the threshold, the shorter sequence is removed.
# In 'global' mode, all sequences are compared regardless of species, and shorter sequences
# with high similarity are removed. The remaining sequences are saved to the output file.
#
# Example:
# ./script.sh -i input.fasta -t 90 -o output.fasta -m species -s 2
# In this example, the script will use the first 2 characters of each sequence header to
# identify the species and remove shorter sequences within the same species that have 90% or
# higher identity. Use 'global' mode for cross-species comparison without requiring the species length.
#

# Function to display usage
show_help() {
    echo "Usage: $0 -i input_fasta -t threshold -o output_fasta -m mode [-s species_length]"
    echo ""
    echo "Arguments:"
    echo "  -i input_fasta: The input FASTA file containing sequences"
    echo "  -t threshold: The percent identity threshold for removing sequences (e.g., 90 for 90% identity)"
    echo "  -o output_fasta: The output FASTA file where filtered sequences will be saved"
    echo "  -m mode: 'species' for species-specific comparisons, 'global' for global comparisons"
    echo "  -s species_length: (Optional) Number of characters from the beginning of each sequence header to use as the species identifier (only required in 'species' mode)"
    echo "  -h: Show help"
    exit 1
}

# Parse command line arguments
while getopts ":i:t:o:m:s:h" opt; do
    case ${opt} in
        i ) input_fasta=$OPTARG ;;
        t ) threshold=$OPTARG ;;
        o ) output_fasta=$OPTARG ;;
        m ) mode=$OPTARG ;;
        s ) species_length=$OPTARG ;;
        h ) show_help ;;
        \? ) echo "Invalid option: -$OPTARG" 1>&2; show_help ;;
        : ) echo "Invalid option: -$OPTARG requires an argument" 1>&2; show_help ;;
    esac
done

# Ensure required arguments are provided
if [ -z "$input_fasta" ] || [ -z "$threshold" ] || [ -z "$output_fasta" ] || [ -z "$mode" ]; then
    echo "Error: Missing required arguments." 1>&2
    show_help
fi

# Ensure the mode is valid
if [[ "$mode" != "species" && "$mode" != "global" ]]; then
    echo "Error: Mode must be either 'species' or 'global'." 1>&2
    show_help
fi

# Ensure species length is provided in species mode
if [[ "$mode" == "species" && -z "$species_length" ]]; then
    echo "Error: Species length is required in 'species' mode." 1>&2
    show_help
fi

# Arrays to store headers, sequences, lengths, species identifiers, and removal flags
declare -a headers
declare -a sequences
declare -a lengths
declare -a species_ids
declare -a removed

# Variables for current header and sequence
current_header=""
current_sequence=""

# Read the input FASTA file
while read -r line; do
    if [[ $line == \>* ]]; then
        # Store the previous sequence if it exists
        if [ -n "$current_header" ]; then
            headers+=("$current_header")
            sequences+=("$current_sequence")
            seq_length=$(echo "$current_sequence" | tr -d '-' | wc -c)
            lengths+=("$seq_length")
            # Extract species identifier (first n characters after '>')
            if [[ "$mode" == "species" ]]; then
                species_id=$(echo "$current_header" | cut -c2-$((species_length+1)))
            else
                species_id="global"  # In global mode, treat all as the same species
            fi
            species_ids+=("$species_id")
            removed+=(0)
            current_sequence=""
        fi
        current_header="$line"
    else
        current_sequence+="$line"
    fi
done < "$input_fasta"

# Add the last sequence
if [ -n "$current_header" ]; then
    headers+=("$current_header")
    sequences+=("$current_sequence")
    seq_length=$(echo "$current_sequence" | tr -d '-' | wc -c)
    lengths+=("$seq_length")
    if [[ "$mode" == "species" ]]; then
        species_id=$(echo "$current_header" | cut -c2-$((species_length+1)))
    else
        species_id="global"
    fi
    species_ids+=("$species_id")
    removed+=(0)
fi

num_sequences=${#headers[@]}

# Pairwise sequence comparison
for (( i=0; i<$num_sequences; i++ )); do
    if [ ${removed[$i]} -eq 0 ]; then
        seq_i="${sequences[$i]}"
        len_i="${lengths[$i]}"
        species_i="${species_ids[$i]}"
        for (( j=i+1; j<$num_sequences; j++ )); do
            if [ ${removed[$j]} -eq 0 ]; then
                species_j="${species_ids[$j]}"
                # Only compare sequences within species in 'species' mode, compare all in 'global' mode
                if [[ "$mode" == "global" || "$species_i" == "$species_j" ]]; then
                    seq_j="${sequences[$j]}"
                    len_j="${lengths[$j]}"
                    # Calculate percent identity using awk
                    percent_identity=$(awk -v seq1="$seq_i" -v seq2="$seq_j" '
                        BEGIN {
                            len = length(seq1);
                            identical = 0;
                            total = 0;
                            for (k = 1; k <= len; k++) {
                                char1 = substr(seq1, k, 1);
                                char2 = substr(seq2, k, 1);
                                if (char1 != "-" && char2 != "-") {
                                    total++;
                                    if (char1 == char2) {
                                        identical++;
                                    }
                                }
                            }
                            if (total > 0) {
                                pid = (identical / total) * 100;
                                printf "%.2f", pid;
                            } else {
                                printf "0.00";
                            }
                        }
                    ')
                    # Remove the shorter sequence if threshold is met
                    if (( $(echo "$percent_identity >= $threshold" | bc -l) )); then
                        if [ "$len_i" -ge "$len_j" ]; then
                            removed[$j]=1
                        else
                            removed[$i]=1
                            break
                        fi
                    fi
                fi
            fi
        done
    fi
done

# Output the filtered sequences
> "$output_fasta"
for (( i=0; i<$num_sequences; i++ )); do
    if [ ${removed[$i]} -eq 0 ]; then
        echo "${headers[$i]}" >> "$output_fasta"
        echo "${sequences[$i]}" >> "$output_fasta"
    fi
done