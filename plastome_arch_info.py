#!/usr/bin/env python
import sys
import os
from optparse import OptionParser
import time


def get_options():
    parser = OptionParser(usage="plastome_arch_info.py fasta_format_sequence_file(s)")
    parser.add_option("-o", dest="output",
                      help="output file.")
    parser.add_option("-r", dest="min_ir_length", default=5000, type=int,
                      help="The minimum repeat length treated as the IR region of plastome. Default: [%default]")
    parser.add_option("-v", dest="valid_bases", default="ATGCRMYKHBDVatgcrmykhbdv",
                      help="Valid bases. Default: ATGCRMYKHBDVatgcrmykhbdv")
    options, argv = parser.parse_args()
    if not len(argv):
        parser.print_help()
        sys.exit()
    else:
        for f in argv:
            if not os.path.isfile(f):
                raise FileNotFoundError(f + " not found/valid!")
        options.valid_bases = set(list(options.valid_bases))
        return options, argv


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
            seqs.append(this_seq.replace(" ", ""))
        else:
            this_line = fasta_handler.readline()
    fasta_handler.close()
    return [names, seqs, interleaved]

try:
    # python2
    import string
    translator = string.maketrans("ATGCRMYKHBDVatgcrmykhbdv", "TACGYKRMDVHBtacgykrmdvhb")

    def complementary_seq(input_seq):
        return string.translate(input_seq, translator)[::-1]
except AttributeError:
    # python3
    translator = str.maketrans("ATGCRMYKHBDVatgcrmykhbdv", "TACGYKRMDVHBtacgykrmdvhb")

    def complementary_seq(input_seq):
        return str.translate(input_seq, translator)[::-1]


# Hashing methods.
# I naively wrote this function by myself and latter found the name and description of this algorithm
# in Kurtz et al. 2001. REPuter: the manifold applications of repeat analysis on a genomic scale.
def find_exact_repeats(sequence_string, min_repeat_length, circular,
                       accepted_char=set(list("ATGCRMYKHBDVatgcrmykhbdv"))):
    word_size = min(13, min_repeat_length)
    if len(sequence_string) < min_repeat_length:
        return []
    if circular:
        long_sequence = sequence_string + sequence_string[:word_size - 1]
        here_seq = long_sequence
    else:
        here_seq = sequence_string
    here_seq_r = complementary_seq(here_seq)
    here_seq_length = len(here_seq)
    original_seq_length = len(sequence_string)
    """create accepted id set"""
    if accepted_char:
        accepted_id = set()
        count_cal = 0
        for base_count in range(here_seq_length):
            if here_seq[base_count] not in accepted_char:
                count_cal = 0
            else:
                count_cal += 1
                if count_cal >= word_size:
                    accepted_id.add(base_count - word_size + 1)
    else:
        accepted_id = set(list(range(here_seq_length - word_size + 1)))
    #
    """initialization"""
    words_to_index = {}
    index_to_words = {}

    def add_to_words(add_index, this_forward, this_reverse):
        if this_forward in words_to_index:
            words_to_index[this_forward].add((add_index, 1))
            words_to_index[this_reverse].add((add_index, -1))
        else:
            words_to_index[this_forward] = {(add_index, 1)}
            if this_reverse == this_forward:
                words_to_index[this_reverse].add((add_index, -1))
            else:
                words_to_index[this_reverse] = {(add_index, -1)}

    for i in range(0, here_seq_length):
        if i in accepted_id:
            forward_s = here_seq[i:i + word_size]
            reverse_s = here_seq_r[here_seq_length - i - word_size: here_seq_length - i]
            add_to_words(i, forward_s, reverse_s)
            index_to_words[i] = forward_s

    """find repeats"""
    repeat_indices = set()
    for repeat_tuples in words_to_index.values():
        if len(repeat_tuples) >= 2:
            for repeat_index in repeat_tuples:
                repeat_indices.add(repeat_index[0])
    repeat_indices = sorted(list(repeat_indices))
    repeats = []
    points_to_repeats = {}
    last_connection = set()
    len_indices = len(repeat_indices)
    for i in range(len_indices):
        this_index = repeat_indices[i]
        this_word = index_to_words[this_index]
        this_connection = words_to_index[this_word]

        """part 1: dealing with old connection"""
        # Loop 1: find repeats_to_stop
        # Loop 2: delete the pointers pointing to the stopped repeats
        # Loop 3: update the pointers
        # Loop 4: add new pointers if shorter repeats should be continued with less alias
        # Loop 5: update the repeats according to pointers
        repeats_to_stop = {}
        kinds_del_from_points = set()
        for one_connection in last_connection:
            here_id, direction_trans = one_connection
            candidate_new_connect = (here_id + direction_trans, direction_trans)
            if candidate_new_connect not in this_connection:
                # if one_connection in points_to_repeats:
                for repeat_kind, repeat_num in points_to_repeats[one_connection]:
                    if repeat_kind in repeats_to_stop:
                        repeats_to_stop[repeat_kind][repeat_num] = one_connection
                    else:
                        repeats_to_stop[repeat_kind] = {repeat_num: one_connection}
                    kinds_del_from_points.add(repeat_kind)
        for repeat_kind in kinds_del_from_points:
            for now_start, now_go_to, now_direction in repeats[repeat_kind]:
                connection_del_from_points = (now_go_to - (word_size - 1) * (now_direction == 1),
                                              now_direction)
                if connection_del_from_points in points_to_repeats:
                    count_this_group = 0
                    while count_this_group < len(points_to_repeats[connection_del_from_points]):
                        if points_to_repeats[connection_del_from_points][count_this_group][0] == repeat_kind:
                            del points_to_repeats[connection_del_from_points][count_this_group]
                        else:
                            count_this_group += 1
                    if not len(points_to_repeats[connection_del_from_points]):
                        del points_to_repeats[connection_del_from_points]
        for one_connection in last_connection:
            here_id, direction_trans = one_connection
            candidate_new_connect = (here_id + direction_trans, direction_trans)
            if candidate_new_connect in this_connection:
                if one_connection in points_to_repeats:
                    points_to_repeats[candidate_new_connect] = []
                    for one_repeat_id in points_to_repeats[one_connection]:
                        points_to_repeats[candidate_new_connect].append(one_repeat_id)
                    del points_to_repeats[one_connection]
        for repeat_kind in repeats_to_stop:
            if len(repeats[repeat_kind]) - len(repeats_to_stop[repeat_kind]) >= 2:
                repeat_to_be_continued = False
                for repeat_num in range(len(repeats[repeat_kind])):
                    if repeat_num not in repeats_to_stop[repeat_kind]:
                        start_id, go_to_id, go_to_direction = repeats[repeat_kind][repeat_num]
                        new_connect = (go_to_id - (word_size - 1) * (go_to_direction == 1) + go_to_direction,
                                       go_to_direction)
                        if new_connect in this_connection and new_connect not in points_to_repeats:
                            repeat_to_be_continued = True
                            break
                if repeat_to_be_continued:
                    repeats.append([])
                    for inside_repeat_num in range(len(repeats[repeat_kind])):
                        if inside_repeat_num not in repeats_to_stop[repeat_kind]:
                            start_id, go_to_id, go_to_direction = repeats[repeat_kind][inside_repeat_num]
                            repeats[-1].append([start_id, go_to_id, go_to_direction])
                            new_connect = (go_to_id - (word_size - 1) * (go_to_direction == 1) + go_to_direction,
                                           go_to_direction)
                            if new_connect in points_to_repeats:
                                points_to_repeats[new_connect].append((len(repeats) - 1, len(repeats[-1]) - 1))
                            else:
                                points_to_repeats[new_connect] = [(len(repeats) - 1, len(repeats[-1]) - 1)]
        for one_connection in points_to_repeats:
            for repeat_kind, repeat_num in points_to_repeats[one_connection]:
                start_id, previous_id, this_direction = repeats[repeat_kind][repeat_num]
                repeats[repeat_kind][repeat_num][1] += this_direction
        """part 2: dealing with new connection"""
        for one_connection in this_connection:
            here_id, direction_trans = one_connection
            candidate_last_connect = (here_id - direction_trans, direction_trans)
            if candidate_last_connect not in last_connection:
                    repeats.append([])
                    for inside_connection in this_connection:
                        inside_id, inside_direction = inside_connection
                        repeats[-1].append([inside_id + (word_size - 1) * (inside_direction == -1),
                                            inside_id + (word_size - 1) * (inside_direction == 1),
                                            inside_direction])
                        if (inside_id, inside_direction) in points_to_repeats:
                            points_to_repeats[inside_connection].append((len(repeats) - 1, len(repeats[-1]) - 1))
                        else:
                            points_to_repeats[inside_connection] = [(len(repeats) - 1, len(repeats[-1]) - 1)]
                    break

        if i + 1 < len_indices:
            next_index = repeat_indices[i + 1]
        else:
            next_index = None
        if next_index != this_index + 1:
            points_to_repeats = {}
            last_connection = set()
        else:
            last_connection = this_connection

    """aftertreatment"""
    # delete repeated repeats
    final_repeat = []
    repeat_dicts = set()
    for repeat_group in repeats:
        if tuple(repeat_group[0]) not in repeat_dicts:
            for single_repeat in repeat_group:
                here_start, here_end, here_direction = single_repeat
                repeat_dicts.add((here_start, here_end, here_direction))
                repeat_dicts.add((here_end, here_start, -here_direction))
            final_repeat.append(repeat_group)
        else:
            continue
    if circular:
        # merge repeats of two ends

        def from_start_or_not(repeat_group_inside):
            for c_from, c_to, c_direction in repeat_group_inside:
                if c_from == 0:
                    return True
            return False

        def to_end_or_not(repeat_group_inside):
            for c_from, c_to, c_direction in repeat_group_inside:
                if (c_from, c_to)[c_direction == 1] == here_seq_length - 1:
                    return True
            return False

        def be_merging_or_not(candidate_repeat_1, candidate_repeat_2):
            inside_1_from, inside_1_to, inside_1_direction = candidate_repeat_1
            inside_2_from, inside_2_to, inside_2_direction = candidate_repeat_2
            if inside_1_direction != inside_2_direction:
                return False
            else:
                if inside_1_to == here_seq_length - 1 and (inside_2_from, inside_2_to)[inside_2_direction == -1] == 0:
                    return True
                else:
                    # print((inside_1_to - inside_2_from+inside_1_direction)*inside_1_direction)
                    if (inside_1_to - inside_2_from+inside_1_direction)*inside_1_direction == word_size - 1:
                        return True
                    else:
                        return False

        count_from_start = 0
        while count_from_start < len(final_repeat) and from_start_or_not(final_repeat[count_from_start]):
            count_from_start += 1
        possible_ends = set()
        for candidate_group_id in range(len(final_repeat)):
            if to_end_or_not(final_repeat[candidate_group_id]):
                possible_ends.add(candidate_group_id)
        # transform the direction of end to positive
        for i in possible_ends:
            for e_from, e_to, e_direction in final_repeat[i]:
                if (e_from, e_to)[e_direction == 1] == here_seq_length - 1:
                    break
            if e_direction == -1:
                for j in range(len(final_repeat[i])):
                    e_from, e_to, e_direction = final_repeat[i][j]
                    final_repeat[i][j] = [e_to, e_from, -e_direction]

        repeat_group_to_del = set()
        for i in range(count_from_start):
            for j in possible_ends:
                matches = {}
                for candidate_id_1 in range(len(final_repeat[j])):
                    for candidate_id_2 in range(len(final_repeat[i])):
                        if be_merging_or_not(final_repeat[j][candidate_id_1], final_repeat[i][candidate_id_2]):
                            matches[candidate_id_1] = candidate_id_2
                            break
                # print('matches', matches)
                if len(matches) < 2:
                    # no effective matches
                    continue
                elif len(matches) == len(final_repeat[j]) == len(final_repeat[i]):
                    # all matched, merge
                    repeat_group_to_del.add(i)
                    for merge_1_id, merge_2_id in matches.items():
                        final_repeat[j][merge_1_id][1] = final_repeat[i][merge_2_id][1]
                elif len(matches) == len(final_repeat[j]):
                    # all of one side matches, extend this side
                    for merge_1_id, merge_2_id in matches.items():
                        final_repeat[j][merge_1_id][1] = final_repeat[i][merge_2_id][1]
                elif len(matches) == len(final_repeat[i]):
                    # all of next side matches, extend next side
                    for merge_1_id, merge_2_id in matches.items():
                        final_repeat[i][merge_2_id][0] = final_repeat[j][merge_1_id][0]
                else:
                    # some of the repeat matches, create a new repeat group
                    final_repeat.append([])
                    for merge_1_id, merge_2_id in matches.items():
                        final_repeat[-1].append([final_repeat[j][merge_1_id][0],
                                                 final_repeat[i][merge_2_id][1],
                                                 final_repeat[i][merge_2_id][2]])
        # print('repeat_group_to_del', repeat_group_to_del)
        repeat_group_to_del = sorted(list(repeat_group_to_del), reverse=True)
        for id_to_del in repeat_group_to_del:
            del final_repeat[id_to_del]

        # circularize number
        for count_group_ in range(len(final_repeat)):
            for j in range(len(final_repeat[count_group_])):
                here_start, here_end, here_direction = final_repeat[count_group_][j]
                final_repeat[count_group_][j] = [here_start % original_seq_length,
                                                 here_end % original_seq_length,
                                                 here_direction]

    # delete small repeats
    count_group__ = 0
    while count_group__ < len(final_repeat):
        here_start, here_end, here_direction = final_repeat[count_group__][0]
        if ((here_end - here_start + here_direction)*here_direction) % original_seq_length < min_repeat_length:
            del final_repeat[count_group__]
        else:
            count_group__ += 1
    # sort according to occurrence
    for group_to_sort in range(len(final_repeat)):
        start, end, direction = final_repeat[group_to_sort][0]
        if start == end:
            this_len = 1
        else:
            if (end - start) * direction > 0:
                this_len = (end - start) * direction + 1
            else:
                this_len = (start, end)[direction == 1] + original_seq_length - (start, end)[direction != 1] + 1
        # transform into dict
        final_repeat[group_to_sort] = [{"start": start, "end": end, "direction": direction, "length": this_len}
                                       for start, end, direction in final_repeat[group_to_sort]]
    # sort according to length: from longest to shortest
    final_repeat.sort(key=lambda x: -x[0]["length"])
    return final_repeat


def reverse_repeats_info(repeats):
    new_repeats = []
    for rep in repeats:
        new_repeats.append({"start": rep["end"], "end": rep["start"],
                            "direction": -rep["direction"], "length": rep["length"]})
    return new_repeats


def detect_architecture(sequence, min_repeat_length, accepted_char):
    # assume the longest is ir
    all_repeats = find_exact_repeats(sequence, min_repeat_length, True, accepted_char)
    if all_repeats:
        # Sorting makes:
        # direct1==1 and direct2==-1
        # start1 be the smallest forward start
        ir_locations_1 = sorted(all_repeats[0], key=lambda x: -x["direction"])
        ir_locations_2 = sorted(reverse_repeats_info(ir_locations_1), key=lambda x: -x["direction"])
        ir_locations = sorted([ir_locations_1, ir_locations_2], key=lambda x: x[0]["start"])[0]
        if len(ir_locations) != 2 or ir_locations[0]["direction"] == ir_locations[1]["direction"]:
            return "-", "-", "-", "no IR found"
        else:
            start1, end1, direct1 = ir_locations[0]["start"], ir_locations[0]["end"], ir_locations[0]["direction"]
            start2, end2, direct2 = ir_locations[1]["start"], ir_locations[1]["end"], ir_locations[1]["direction"]
            # cross the end, meaning site:seq_len in (IR1)
            if end1 < start1:
                # seq_len >= start2 >= start1
                if start2 >= start1:
                    return 0, 0, ir_locations[0]["length"], "IR overlaps"
                else:
                    return start1 - start2 - 1, 0, ir_locations[0]["length"], "IR overlaps"
            elif start2 < end2:
                if start2 >= start1:
                    if end1 >= end2:
                        return 0, 0, ir_locations[0]["length"], "IR overlaps"
                    else:
                        return end2 - end1 - 1, 0, ir_locations[0]["length"], "IR overlaps"
                else:
                    ssc, lsc = sorted([end2 - end1 - 1, start1 - start2 - 1])
                    return lsc, ssc, ir_locations[0]["length"], "-"
            else:
                ssc, lsc = sorted([end2 - end1 - 1, len(sequence) + start1 - start2 - 1])
                return lsc, ssc, ir_locations[0]["length"], "-"
    else:
        return "-", "-", "-", "no IR found"


def main():
    time0 = time.time()
    sys.stdout.write("\n"
                     "## This script helps you count the LSC/SSC/IR lengths from a batch of plastome sequences.\n"
                     "## by jinjianjun@mail.kib.ac.cn\n\n")
    options, argv = get_options()
    sys.stdout.write("file\tsequence_name\ttotal\tLSC\tSSC\tIR\tNotes\n")
    if options.output:
        out_handler = open(options.output, "w")
        out_handler.close()
    for this_f in argv:
        this_matrix = read_fasta(this_f)
        for i in range(len(this_matrix[0])):
            arch = detect_architecture(this_matrix[1][i], options.min_ir_length, options.valid_bases)
            out_line = "\t".join([this_f, this_matrix[0][i], str(len(this_matrix[1][i]))] +
                                 [str(x) for x in arch]) + "\n"
            sys.stdout.write(out_line)
            if options.output:
                open(options.output, "a").write(out_line)
    sys.stdout.write("\n## Cost: " + str(round(time.time() - time0, 2)) + "s\n\n")

if __name__ == '__main__':
    main()
