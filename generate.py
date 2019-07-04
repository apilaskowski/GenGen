#!/usr/bin/env python3

from random import choice
from random import random
from random import randint
import numpy as np
import argparse
import json
from pprint import pprint
import sys

nucleotides = ['A', 'C', 'T', 'G']
nucleotides_swaps = {
    'A': ['C', 'T', 'G'], 
    'C': ['A', 'T', 'G'],
    'T': ['A', 'C', 'G'],
    'G': ['A', 'C', 'T']
}

complementary = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}

line_length = 70

def generate_nucleotide():
    return choice(nucleotides)

def generate_genome(length):
    genome = ''
    for _ in range(length):
        genome += generate_nucleotide()
    return genome

def generate_repetitive_genome(genome, repetition_length, repetition_number):
    repetition_beg = randint(0, len(genome) - repetition_length)
    print("Repetition number:", repetition_number, ", of length:", repetition_length)
    print("Repetition pos:", repetition_beg)
    repetition = genome[repetition_beg:repetition_beg+repetition_length]
    for _ in range(repetition_number):
        random_pos = randint(0, len(genome) - repetition_length)
        print(_, "th repetition: ", random_pos) 
        genome = genome[:random_pos] + repetition + genome[random_pos+repetition_length:]
    return genome

def calc_read_number(genome_length, read_length, coverage):
    max_coverage = ((genome_length-read_length)*read_length)/genome_length
    print("Max coverage: " + str(max_coverage))
    print("Your coverage: " + str(coverage))
    if coverage > max_coverage:
        coverage = max_coverage
    return int(coverage * genome_length / read_length)

def choice_reads_without_replacement(genome, read_number, read_length):
    reads = []
    positions = np.random.choice(len(genome) - read_length, read_number, replace=False)
    for pos in positions:
        reads += [genome[pos:pos+read_length]]
    return reads, positions

def insert_errors(reads, error_rate):
    reads_with_errors = []
    for read in reads:
        read_with_errors = ''
        for nucleotide in read:
            rullete = random()
            if rullete < error_rate:
                index = randint(0, 2)
                read_with_errors += nucleotides_swaps[nucleotide][index]
            else:
                read_with_errors += nucleotide
        reads_with_errors.append(read_with_errors)
    return reads_with_errors

def generate_pairs(genome, positions, read_length, insert):
    paired_reads = []
    for pos in positions:
        sign = np.random.choice([-1, 1])
#        print(pos, "\t", sign, file=sys.stderr, end='\t')
        if pos + sign * insert + read_length > len(genome) or pos + sign * insert < 0:
            sign *= -1
#        print(sign, file=sys.stderr, end='\t')
        pair_pos = pos + sign * insert
#        print(pair_pos, file=sys.stderr, end='\t')
#        print(genome[pair_pos:pair_pos+read_length], file=sys.stderr)
        paired_reads.append(genome[pair_pos:pair_pos+read_length])
    return paired_reads

def save_element_to_file(file, id, element):
    file.write('>' + str(id) + '\n')
    for i, nucleotide in enumerate(element):
        file.write(nucleotide)
        if (i + 1) % line_length == 0:
            file.write('\n')
    file.write('\n')

def reverse_complementary_reads(reads):
	result = []
	for read in reads:
		new_read = ""
		for let in read:
			new_read += complementary[let]
		result.append(new_read[::-1])
	return result
			

def save_genome(genome_file_path, genome, suffix='.fasta'):
    genome_file = open(genome_file_path + suffix, 'w')
    save_element_to_file(genome_file, 0, genome)
    genome_file.close()

def save_reads(reads_file_path, reads, suffix='.fasta'):
    reads_file = open(reads_file_path + suffix, 'w')
    for i, r in enumerate(reads):
        save_element_to_file(reads_file, i, r)
    reads_file.close()

def main():
    parser = argparse.ArgumentParser(description='Generate synthetic genome with associated reads.')
    parser.add_argument('--configuration-file', help='Configuration file path', required=True)
    parser.add_argument('--genome-file', help='Where to write genome to', required=True)
    parser.add_argument('--reads-file', help='Where to write reads to', required=True)
    args = parser.parse_args()

    config = json.load(open(args.configuration_file))

    global line_length
    line_length = config["LINE_LENGTH"]

    genome = generate_genome(config["GENOME_LENGTH"])
    is_repetitive = bool(config["IS_REPETITIVE"])
    if is_repetitive:
        repetition_length = config["REPETITION_LENGTH"]
        repetition_number = config["REPETITION_NUMBER"]
        for idx in range(len(repetition_length)):
            genome = generate_repetitive_genome(genome, int(repetition_length[idx]), int(repetition_number[idx]))

    read_number = calc_read_number(config["GENOME_LENGTH"], config["READ_LENGTH"], float(config["COVERAGE"]))
    reads, positions = choice_reads_without_replacement(genome, read_number, config["READ_LENGTH"])
    reads_with_errors = insert_errors(reads, float(config["ERROR_RATE"]))

    if bool(config["ORDER_READS"]):
        reads = [read for _, read in sorted(zip(positions, reads))]
        reads_with_errors = [read for _, read in sorted(zip(positions, reads_with_errors))]
        positions = sorted(positions)

    save_genome(args.genome_file, genome)
    save_reads(args.reads_file, reads, '_clean_1' + config["SUFFIX"])
    save_reads(args.reads_file, reads_with_errors, '_1' + config["SUFFIX"])

    if bool(config["PAIRED"]):
        paired_reads = generate_pairs(genome, positions, config["READ_LENGTH"], int(config["INSERT_SIZE"]))
        paired_reads_with_errors = insert_errors(paired_reads, float(config["ERROR_RATE"]))
        if bool(config["REVCOMP"]):
            paired_reads_revcomp = reverse_complementary_reads(paired_reads)
            paired_reads_with_errors_revcomp = reverse_complementary_reads(paired_reads_with_errors)
			
            save_reads(args.reads_file, paired_reads_revcomp, '_clean_2' + config["SUFFIX"])
            save_reads(args.reads_file, paired_reads_with_errors_revcomp, '_2' + config["SUFFIX"])
        else:
            save_reads(args.reads_file, paired_reads, '_clean_2' + config["SUFFIX"])
            save_reads(args.reads_file, paired_reads_with_errors, '_2' + config["SUFFIX"])

if "__main__":
    main()
