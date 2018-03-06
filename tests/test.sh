#!/usr/bin/env bash
set -ev

# Test script functionality
Rscript seqConsolidate.R tests/testSeq-1.R2.filt.fastq.gz \
    -o tests/testSeq-1.R2.consol.fasta -k tests/testSeq-1.R2.key.csv \
    -l testSeq1. --compress

# Check output for correct processing, 50 keys, 25 unique sequences
test_key_len=$(cat tests/testSeq-1.R2.key.csv | sed '/readNames/d' | wc -l)
test_seq_len=$(zcat tests/testSeq-1.R2.consol.fasta.gz | sed -n '2~2p' | wc -l)

if [ ! $test_key_len = 50 ] | [ ! $test_seq_len = 25 ]; then
    exit 1
fi

# Beginning sequences from test output
zcat tests/testSeq-1.R2.consol.fasta.gz | sed -n '2~2p' | head -n 5

# Beginning of key file for sequences
cat tests/testSeq-1.R2.key.csv | head -n 6

# Cleanup 
rm tests/*.consol.fasta* tests/*.key.csv

echo "Passed all tests."
exit
