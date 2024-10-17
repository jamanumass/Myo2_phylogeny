#!/bin/bash

# Step 1: Download and Extract the Pfam-A HMM Profiles
wget -O Pfam-A.hmm.gz "http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
gunzip Pfam-A.hmm.gz

# Step 2: Extract the Specific HMM Profile for Myosin head, motor domain (PF00063)
awk '/^ACC   PF00063/,/^\/\//' Pfam-A.hmm > myosin_head_motor_domain.hmm