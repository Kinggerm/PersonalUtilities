#!/usr/bin/env python
import os
import time
import sys
from optparse import OptionParser


def require_options():
    usage = "python this_script.py --fa1 fasta1 --fa2 fasta2"
    parser = OptionParser(usage=usage)
    parser.add_option('--fa1', dest='fasta1',
                      help='')
    parser.add_option('--fa2', dest='fasta2',
                      help='')
    parser.add_option('--phy1', dest='phylip1',
                      help='')
    parser.add_option('--phy2', dest='phylip2',
                      help='')
    options, args = parser.parse_args()
    if not ((options.fasta1 and options.fasta2) or (options.phylip1 and options.phylip2)) or \
            (options.fasta1 and options.fasta2) and (options.phylip1 and options.phylip2):
        parser.print_help()
        sys.stdout.write('\n######################################\nERROR: Invalid arguments!\n\n')
        exit()
    return options, args


def read_fasta(fasta_file):
    fasta_handler = open(fasta_file)
    names = []
    seqs = []
    this_line = fasta_handler.readline()
    interleaved = 0
    while this_line:
        if this_line.startswith('>'):
            names.append(this_line[1:].strip())
            this_seq = ''
            this_line = fasta_handler.readline()
            seq_line_count = 0
            while this_line and not this_line.startswith('>'):
                if seq_line_count == 1:
                    interleaved = len(this_seq)
                this_seq += this_line.strip()
                this_line = fasta_handler.readline()
                seq_line_count += 1
            seqs.append(this_seq)
        else:
            this_line = fasta_handler.readline()
    fasta_handler.close()
    return [names, seqs, interleaved]


def read_phylip(phylip_file):
    phylip_handler = open(phylip_file)
    names = []
    seqs = []
    head = phylip_handler.readline()
    for line in phylip_handler:
        this_name, this_seq = line.split()[:2]
        names.append(this_name)
        seqs.append(this_seq.strip())
    return [names, seqs, 0]


def transpose(matrix):
    return [[matrix[j][i] for j in range(len(matrix))] for i in range(len(matrix[0]))]


def sort_column_mode(matrix, built=lambda x: x):
    new_matrix = transpose(matrix)
    new_matrix.sort(key=built)
    return transpose(new_matrix)


def compare_site_patterns(matrix_1, matrix_2):
    if not matrix_1[0]:
        print("no species")
        sys.exit()
    if set(matrix_1[0]) != set(matrix_2[0]):
        print("taxa not match")
        sys.exit()
    if not (len(set([len(seq) for seq in matrix_1[1]])) == len(set([len(seq) for seq in matrix_2[1]])) == 1):
        print("not aligned")
        sys.exit()
    new_matrix_1 = sort_column_mode(matrix_1[:2])
    new_matrix_2 = sort_column_mode(matrix_2[:2])
    patterns_1 = {"id": [], "pattern": {}}
    patterns_2 = {"id": [], "pattern": {}}
    for go_to_site in range(len(new_matrix_1[1][0])):
        pattern_1 = []
        pattern_2 = []
        for go_to_sp in range(len(new_matrix_1[0])):
            pattern_1.append(new_matrix_1[1][go_to_sp][go_to_site])
            pattern_2.append(new_matrix_2[1][go_to_sp][go_to_site])
        pattern_1 = tuple(pattern_1)
        pattern_2 = tuple(pattern_2)
        patterns_1["id"].append(pattern_1)
        patterns_2["id"].append(pattern_2)
        if pattern_1 in patterns_1["pattern"]:
            patterns_1["pattern"][pattern_1].append(go_to_site)
        else:
            patterns_1["pattern"][pattern_1] = [go_to_site]
        if pattern_2 in patterns_2["pattern"]:
            patterns_2["pattern"][pattern_2].append(go_to_site)
        else:
            patterns_2["pattern"][pattern_2] = [go_to_site]
    all_patterns = {}
    for this_pattern in set(patterns_1["pattern"]) | set(patterns_2["pattern"]):
        all_patterns[this_pattern] = {1: len(patterns_1["pattern"].get(this_pattern, [])),
                                      2: len(patterns_2["pattern"].get(this_pattern, []))}
    return all_patterns, new_matrix_1[0]


def color_it(print_it, determine):
    if determine:
        return '\033[93m' + print_it + '\033[0m'
    else:
        return '\033[91m' + print_it + '\033[0m'


def main():
    # time0 = time.time()
    options, args = require_options()
    if options.fasta1:
        matrix1 = read_fasta(options.fasta1)
        matrix2 = read_fasta(options.fasta2)
    else:
        matrix1 = read_phylip(options.phylip1)
        matrix2 = read_phylip(options.phylip2)
    statistics, name_list = compare_site_patterns(matrix1, matrix2)
    difference = 0
    patterns = sorted(list(statistics.keys()))
    total = len(patterns)
    result = [["#pattern", "count in data1", "count in data2", "pattern"]]
    for i in range(total):
        count_1 = statistics[patterns[i]][1]
        count_2 = statistics[patterns[i]][2]
        if count_1 != count_2:
            difference += 1
        result.append([str(i + 1), str(count_1), str(count_2), "".join(patterns[i])])
    sys.stdout.write("\n")
    print_width = [max([len(row[column]) for row in result]) + 1 for column in range(len(result[0]))]
    for col in range(len(result[0])):
        sys.stdout.write(result[0][col].rjust(print_width[col]))
    sys.stdout.write("\n")
    for row in result[1:]:
        match = row[1] == row[2]
        for col in range(len(result[0])):
            sys.stdout.write(color_it(row[col].rjust(print_width[col]), match))
        sys.stdout.write("\n")
    sys.stdout.write("\nOrder for pattern: "+",".join(name_list)+"\n")
    sys.stdout.write("Matched/Unmatched: "+str(total - difference)+"/"+str(difference)+"\n\n")
    # print("Cost: "+str(time.time()-time0))


if __name__ == '__main__':
    main()
