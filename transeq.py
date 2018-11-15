#!/usr/bin/env python
import sys
import os
import time
from Bio import SeqIO
from Bio import Alphabet
time0 = time.time()
postfix_dict = {'phylip': 'phy', 'phylip-relaxed': 'phy', 'nexus': 'nex', 'fastq': 'fq', 'genbank': 'gb'}
if len(sys.argv) < 4: # or sys.argv[1] in {'-h', '--help'}:
    print("Usage:   transeq.py from_format to_format input_file(s)/input_folder(s)\
          \nFormat:  phylip, phylip-relaxed, fasta, nexus, gb, fastq ...\n")
else:
    # import glob
    # for f in glob.glob(sys.argv[3]):
    count_pass = 0
    count_fall = 0
    input_format, out_format = sys.argv[1], sys.argv[2]
    in_postfix = postfix_dict.get(input_format, input_format)
    out_postfix = postfix_dict.get(out_format, out_format)
    for f in sys.argv[3:]:
        if os.path.isfile(f):
            new_f = os.path.join(os.path.split(f)[0], '.'.join(os.path.split(f)[1].split('.')[:-1]) + '.' + out_postfix)
            try:
                seqs = SeqIO.parse(open(f, 'rU'), input_format, alphabet=Alphabet.generic_dna)
                SeqIO.write(seqs, open(new_f, 'w'), out_format)
                count_pass += 1
            except ValueError as e:
                count_fall += 1
                if 'Unknown format' in str(e):
                    print("\nError:", str(e))
                    os.remove(new_f)
                    break
                else:
                    os.remove(new_f)
                    raise e
        elif os.path.isdir(f):
            for sub_f in [os.path.join(f, x) for x in os.listdir(f) if x.endswith(in_postfix)]:
                new_f = os.path.join(os.path.split(sub_f)[0],
                                     '.'.join(os.path.split(sub_f)[1].split('.')[:-1]) + '.' + out_postfix)
                try:
                    seqs = SeqIO.parse(open(sub_f, 'rU'), input_format, alphabet=Alphabet.generic_dna)
                    SeqIO.write(seqs, open(new_f, 'w'), out_format)
                    count_pass += 1
                except ValueError as e:
                    count_fall += 1
                    if 'Unknown format' in str(e):
                        print("\nError:", str(e))
                        os.remove(new_f)
                        break
                    else:
                        os.remove(new_f)
                        raise e
        else:
            count_fall += 1
    print(str(count_pass) + "/" + str(count_pass + count_fall), "file(s) cost", time.time()-time0)
