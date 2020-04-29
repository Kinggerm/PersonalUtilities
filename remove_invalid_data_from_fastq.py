#!/usr/bin/env python
import sys
import os
import time
time0 = time.time()
if len(sys.argv) != 3:
    print("Usage: python remove_invalid_data_from_fastq.py raw.fastq new.fastq")
assert os.path.isfile(sys.argv[1]), sys.argv[1] + " not found!"
# assert not os.path.isfile(sys.argv[2]), sys.argv[2] + " existed!"

valid_bases = set("ATGCRMYKHBDVN")
# assume Illumina 1.8+ Phred+33 format
# valid_quals = set("""!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ""")
# assume Illumina Sanger Phred+33 format
# valid_quals = set("""!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ""")
valid_quals = set("""!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghi""")

with open(sys.argv[2] + ".dumped.txt", "w") as error_f:
    with open(sys.argv[2], "w") as output_f:
        with open(sys.argv[1]) as input_f:
            go = 1
            print_s = "Checked lines: " + str(go)
            line = input_f.readline()
            lines = []
            while line:
                if line.startswith("@"):  # valid head
                    lines.append(line)
                    line = input_f.readline()
                    go += 1
                    if set(line.strip()).issubset(valid_bases):  # valid sequence
                        lines.append(line)
                        line = input_f.readline()
                        go += 1
                        if line.strip() == "+":
                            lines.append(line)
                            line = input_f.readline()
                            go += 1
                            if set(line.strip("\n").strip("\r")).issubset(valid_quals) and len(line.strip("\n").strip("\r")) == len(lines[1].strip("\n").strip("\r")):
                                lines.append(line)
                                output_f.writelines(lines)
                                line = input_f.readline()
                                go += 1
                            else:
                                error_f.writelines(["Line " + str(go + back_s - len(lines)) + "\t" + lines[back_s] for back_s in range(len(lines))])
                        else:
                            error_f.writelines(["Line " + str(go + back_s - len(lines)) + "\t" + lines[back_s] for back_s in range(len(lines))])
                    else:
                        error_f.write("Line " + str(go - 1) + "\t" + lines[-1])
                else:
                    error_f.write("Line " + str(go) + "\t" + line)
                    line = input_f.readline()
                    go += 1
                lines = []
                if go % 10000:
                    print_s = "Checked lines: " + str(go)
                    sys.stdout.write(print_s + "\b" * len(print_s))
                    sys.stdout.flush()
            print_s = "Checked lines: " + str(go)
            sys.stdout.write(print_s + "\n")
print("finished in " + str(round(time.time() - time0, 2)) + "s")
    