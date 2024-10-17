#!/bin/bash

# Usage check
if [ $# -ne 3 ]; then
    echo "Usage: $0 input_fasta threshold output_fasta"
    exit 1
fi

input_fasta="$1"
threshold="$2"
output_fasta="$3"

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
            # Extract species identifier (first 7 characters after '>')
            species_id=$(echo "$current_header" | cut -c2-8)
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
    species_id=$(echo "$current_header" | cut -c2-8)
    species_ids+=("$species_id")
    removed+=(0)
fi

num_sequences=${#headers[@]}

# Pairwise sequence comparison within species
for (( i=0; i<$num_sequences; i++ )); do
    if [ ${removed[$i]} -eq 0 ]; then
        seq_i="${sequences[$i]}"
        len_i="${lengths[$i]}"
        species_i="${species_ids[$i]}"
        for (( j=i+1; j<$num_sequences; j++ )); do
            if [ ${removed[$j]} -eq 0 ]; then
                species_j="${species_ids[$j]}"
                # Only compare sequences from the same species
                if [ "$species_i" == "$species_j" ]; then
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