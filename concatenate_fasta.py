#!/usr/bin/env python
import os
import time
import sys
import re
from optparse import OptionParser


def require_options(print_title):
    usage = "python this_script.py a_list_of_fasta_files -o output.fasta"
    parser = OptionParser(usage=usage)
    parser.add_option('-o', dest='output',
                      help='output file fasta file')
    parser.add_option('--separate', dest='aligned', default=True, action='store_false',
                      help='By default the input fasta is treated as alignment. Choose to treat as separate sequences.')
    parser.add_option('--sort', dest='sort_file_name', default=False, action='store_true',
                      help='By default this script would concatenate the file by the original argument order, '
                           'which may sometimes cause disorder when using wildcard (*). Choose to reorder the '
                           'concatenate order by file names.')
    parser.add_option('--config', dest='configuration_file',
                      help='If chosen and the input are aligned, output a configuration file recording the locations.')
    options, args = parser.parse_args()
    if not (options.output and args):
        parser.print_help()
        sys.stdout.write('\n######################################\nERROR: Insufficient REQUIRED arguments!\n\n')
        exit()
    print(print_title)
    if options.sort_file_name:
        min_len = min([len(file_name) for file_name in args])
        go_to = 0
        while go_to < min_len:
            this_character = args[0][go_to]
            is_equal = True
            for go_file in range(1, len(args)):
                if args[go_file][go_to] != this_character:
                    is_equal = False
            if not is_equal:
                break
            go_to += 1
        go_back = 1
        while go_back <= min_len:
            this_character = args[0][-go_back]
            is_equal = True
            for go_file in range(1, len(args)):
                if args[go_file][-go_back] != this_character:
                    is_equal = False
            if not is_equal:
                break
            go_back += 1
        try:
            args.sort(key=lambda x: int(x[go_to: len(x)-go_back+1]))
        except ValueError:
            args.sort()
        sys_args = sys.argv
        count_argv = 0
        sorted_argvs = set(args)
        while count_argv < len(sys_args):
            if sys_args[count_argv] in sorted_argvs:
                del sys_args[count_argv]
            else:
                count_argv += 1
        print(' '.join(sys_args)+' '+' '.join(args) + '\n')
    else:
        print(' '.join(sys.argv) + '\n')
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


def write_fasta(out_file, matrix, overwrite):
    if not overwrite:
        while os.path.exists(out_file):
            out_file = '.'.join(out_file.split('.')[:-1]) + '_.' + out_file.split('.')[-1]
    out_handler = open(out_file, 'w')
    if matrix[2]:
        for i in range(len(matrix[0])):
            out_handler.write('>'+matrix[0][i]+'\n')
            j = matrix[2]
            while j < len(matrix[1][i]):
                out_handler.write(matrix[1][i][(j-matrix[2]):j]+'\n')
                j += matrix[2]
            out_handler.write(matrix[1][i][(j-matrix[2]):j]+'\n')
    else:
        for i in range(len(matrix[0])):
            out_handler.write('>'+matrix[0][i] + '\n')
            out_handler.write(matrix[1][i] + '\n')
    out_handler.close()


def main():
    time0 = time.time()
    print_title = ""
    options, args = require_options(print_title)
    fasta_matrices = []
    seq_names_list = []
    seq_names_set = set()
    if options.aligned:
        lengths = []
    for fasta_file in args:
        this_matrix = read_fasta(fasta_file)
        fasta_matrices.append(this_matrix)
        for seq_name in this_matrix[0]:
            if seq_name not in seq_names_set:
                seq_names_list.append(seq_name)
                seq_names_set.add(seq_name)
        if options.aligned and this_matrix[1]:
            lengths.append(len(this_matrix[1][0]))
            for i in range(1, len(this_matrix[1])):
                if len(this_matrix[1][i]) != lengths[-1]:
                    print("Error: Unequal length between "+this_matrix[0][0]+" and "+this_matrix[0][i] +
                          " in "+fasta_file+"!")
                    exit()
    out_dict = {in_seq_name: '' for in_seq_name in seq_names_list}
    if options.aligned:
        if options.configuration_file:
            config_file = open(options.configuration_file, 'w')
            go_base = 0
        for i in range(len(fasta_matrices)):
            for j in range(len(fasta_matrices[i][0])):
                this_seq_name = fasta_matrices[i][0][j]
                out_dict[this_seq_name] += fasta_matrices[i][1][j]
                seq_names_set.remove(this_seq_name)
            for add_seq_name in seq_names_set:
                out_dict[add_seq_name] += "-"*lengths[i]
            if options.configuration_file:
                config_file.write(re.sub(".fasta$", '', args[i]) +
                                  '\t' + str(go_base + 1) + '-' + str(go_base + lengths[i]) + '\n')
                go_base += lengths[i]
            seq_names_set = set(seq_names_list)
    else:
        for this_matrix in fasta_matrices:
            for i in range(len(this_matrix[0])):
                out_dict[this_matrix[0][i]] += this_matrix[1][i]
    out_matrix = [seq_names_list, [out_dict[this_name] for this_name in seq_names_list], fasta_matrices[0][2]]
    write_fasta(options.output, out_matrix, False)
    print("Cost: "+str(time.time()-time0))


if __name__ == '__main__':
    main()
