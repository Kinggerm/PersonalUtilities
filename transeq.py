#!/usr/bin/env python
import sys
import os
import time
from Bio import SeqIO
from Bio import Alphabet
time0 = time.time()
postfix_dict = {'phylip':'phy', 'phylip-relaxed':'phy', 'nexus':'nex', 'fastq':'fq'}
if len(sys.argv) < 4: # or sys.argv[1] in {'-h', '--help'}:
    print("Usage:   transeq.py from_format to_format file\
          \nFormat:  phylip, phylip-relaxed, fasta, nexus, gb, fastq ...\n")
else:
    # import glob
    # for f in glob.glob(sys.argv[3]):
    for f in sys.argv[3:]:
        try:
            input_format, out_format = sys.argv[1], sys.argv[2]
            postfix = postfix_dict.get(out_format, out_format)
            seqs = SeqIO.parse(open(f, 'rU'), input_format, alphabet=Alphabet.generic_dna)
            new_f = os.path.join(os.path.split(f)[0],'.'.join(os.path.split(f)[1].split('.')[:-1])+'.'+postfix)
            SeqIO.write(seqs, open(new_f, 'w'), out_format)
        except ValueError as e:
            if 'Unknown format' in str(e):
                print("\nError:", str(e))
                os.remove(new_f)
                break
            else:
                os.remove(new_f)
                raise e
print(len(sys.argv[3:]), "file(s) cost", time.time()-time0)
