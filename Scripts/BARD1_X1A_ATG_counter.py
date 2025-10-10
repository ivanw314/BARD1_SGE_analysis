#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  2 14:29:59 2025

@author: ivan
"""

#!/usr/bin/env python3

import argparse
import re

def read_sequences_file(sequences_file):
    """Read the file containing sequences to search for."""
    with open(sequences_file, 'r') as f:
        return [line.strip() for line in f if line.strip()]

def count_matching_reads(fastq_file, sequences):
    """Count reads in FASTQ file that contain any of the sequences."""
    # Compile regular expressions for each sequence
    patterns = [re.compile(re.escape(seq)) for seq in sequences]
    
    # Track counts for each sequence
    counts = {seq: 0 for seq in sequences}
    total_reads = 0
    matching_reads = 0
    
    # For non-gzipped files
    with open(fastq_file, 'r') as f:
        line_num = 0
        current_read = ""
        
        for line in f:
            line_num += 1
            # In FASTQ format, sequence lines are line 2, 6, 10, etc. (1-based)
            if line_num % 4 == 2:
                total_reads += 1
                current_read = line.strip()
                
                # Check if the read contains any of our sequences
                read_matched = False
                for i, pattern in enumerate(patterns):
                    if pattern.search(current_read):
                        counts[sequences[i]] += 1
                        read_matched = True
                
                if read_matched:
                    matching_reads += 1
    
    return total_reads, matching_reads, counts

def main():
    parser = argparse.ArgumentParser(description='Count reads containing specific sequences')
    parser.add_argument('fastq_file', help='Input FASTQ file')
    parser.add_argument('sequences_file', help='File containing sequences to search for')
    parser.add_argument('--output', '-o', help='Output TSV file for sequence counts', default='sequence_counts.tsv')
    
    args = parser.parse_args()
    
    # Read sequences to search for
    sequences = read_sequences_file(args.sequences_file)
    print(f"Loaded {len(sequences)} sequences to search for")
    
    # Count matching reads
    total_reads, matching_reads, counts = count_matching_reads(args.fastq_file, sequences)
    
    # Print summary results to console
    print(f"Total reads: {total_reads}")
    print(f"Reads matching any sequence: {matching_reads} ({matching_reads/total_reads*100:.2f}%)")
    
    # Write sequence counts to TSV
    with open(args.output, 'w') as f:
        # Write header
        f.write("Sequence\tCount\n")
        
        # Write each sequence and its count
        for seq, count in counts.items():
            f.write(f"{seq}\t{count}\n")
    
    print(f"Sequence counts written to {args.output}")
    
if __name__ == "__main__":
    main()