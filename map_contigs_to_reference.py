#!/usr/bin/env python

from optparse import OptionParser, OptionGroup
import xml.etree.cElementTree as ET
import subprocess
try:
    import commands
except ImportError:
    pass
import logging
import sys
import os
import time
import re
import math

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
keep_n_len = 0


def simple_log(log, output_base):
    log_simple = log
    for handler in list(log_simple.handlers):
        log_simple.removeHandler(handler)
    log_simple.setLevel(logging.NOTSET)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.NOTSET)
    logfile = logging.FileHandler(os.path.join(output_base, 'Map.contigs.log'), mode='a')
    logfile.setFormatter(logging.Formatter('%(message)s'))
    logfile.setLevel(logging.NOTSET)
    log_simple.addHandler(console)
    log_simple.addHandler(logfile)
    return log_simple


def timed_log(log, output_base):
    log_timed = log
    for handler in list(log_timed.handlers):
        log_timed.removeHandler(handler)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s: %(message)s'))
    console.setLevel(logging.NOTSET)
    logfile = logging.FileHandler(os.path.join(output_base, 'Map.contigs.log'), mode='a')
    logfile.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s: %(message)s'))
    logfile.setLevel(logging.NOTSET)
    log_timed.addHandler(console)
    log_timed.addHandler(logfile)
    return log_timed


def require_options(print_title):
    try:
        # python3
        blast_in_path = subprocess.getstatusoutput('blastn')
    except AttributeError:
        # python2
        blast_in_path = commands.getstatusoutput('blastn')
    if blast_in_path[0] == 32512:
        sys.stdout.write('\nError: blastn not in the path!')
        exit()
    try:
        # python3
        makeblastdb_in_path = subprocess.getstatusoutput('makeblastdb')
    except AttributeError:
        # python2
        makeblastdb_in_path = commands.getstatusoutput('makeblastdb')
    if makeblastdb_in_path[0] == 32512:
        sys.stdout.write('\nError: makeblastdb not in the path!')
        exit()
    usage = "python this_script.py [*.fasta| *.fastg] -r reference.fasta -o output"
    parser = OptionParser(usage=usage)
    group_basic = OptionGroup(parser, "BASIC OPTIONS", "")
    group_basic.add_option('-r', dest='reference_fasta',
                           help='Input reference sequence or matrix in fasta format.'
                                'Only the first sequence would be taken as the reference.')
    group_basic.add_option('-o', dest='output_directory',
                           help='Output directory.')

    group_extension = OptionGroup(parser, "EXTENSION OPTIONS", "")
    group_extension.add_option('--merge-max-gap', dest='max_gap', default=2000, type=int,
                               help='max gap for merging. Default: 2000.')
    group_extension.add_option('--merge-max-dif', dest='max_dif', default=0.5, type=float,
                               help='max difference in length for merging. Default: 0.5.')
    group_extension.add_option('--add-gap-disconnect', dest='add_gap_disconnect', default=10, type=int,
                               help='By default, if separate contigs have successive hits in reference, '
                                    'add extra 10-N there to imply.')
    group_extension.add_option('--min-overlap', dest='min_overlap', default=20, type=int,
                               help='Minimum overlap served as criteria of merging unknown nodes. Default: 20.')
    group_extension.add_option('--expand-percent', dest='expand_percent', default=1.0, type=float,
                               help='Expanding unmerged contigs with certain ratio. Default: 1.00')
    group_extension.add_option('--expand-limit', dest='expand_limit', default=1000000, type=int,
                               help='Expanding unmerged contigs with no more than certain length. Default: 1000000.')

    group_pretreatment = OptionGroup(parser, "PRETREATMENT OPTIONS", "")
    group_pretreatment.add_option('--rm-low', dest='rm_low_coverage', default=False, action='store_true',
                                  help='Call rm_low_coverage_duplicated_contigs.py to pre-treat input fastg file. '
                                       'Default: False.')
    group_pretreatment.add_option('--word-size', dest='word_size', default=50, type=int,
                                  help='Word size used to detect repeats. Default: 50.')
    group_pretreatment.add_option('--min-repeat', dest='min_repeat', default=1000, type=int,
                                  help="Minimum repeat to detect and delete from query and reference fasta file. "
                                       "Typically used to remove IR. Choose 0 to disable, choose >= 50 to enable. "
                                       "Default: 1000.")
    group_pretreatment.add_option('--add-gap-repeat', dest='add_gap_repeat', default=0, type=int,
                                  help='Gap to add to mark the removed repeats. Default: 0.')
    group_pretreatment.add_option('--keep-repeat-ends', dest='keep_repeat_ends', default=25, type=int,
                                  help='Keep ends of certain length when removing repeats areas. '
                                       'A value smaller than half minimum repeat is required. '
                                       'Also a value larger than min_overlap is suggested. Default: 25.')
    group_pretreatment.add_option('--linear-refer', dest='circular_refer', default=True, action='store_false',
                                  help='By default, this script assumes the reference sequence is circular. '
                                       'Choose to make it linear.')
    group_pretreatment.add_option('--linear-query', dest='circular_query', default=True, action='store_false',
                                  help='By default, this script assumes all the fasta query sequence is circular. '
                                       'Choose to make it linear.')

    group_blast = OptionGroup(parser, "BLAST OPTIONS", "")
    group_blast.add_option('--blast-w-1', dest='blast_word_size', default=10, type=int,
                           help='Blastn word size (-word_size) for fasta input. Default: 10.')
    group_blast.add_option('--blast-e-1', dest='blast_evalue', default=1E-20, type=float,
                           help="Blastn evalue (-evalue) for fasta input. Default: 1E-20.")
    group_blast.add_option('--blast-w-2', dest='blast_word_size_fg', default=9, type=int,
                           help='Blastn word size (-word_size) for fastg input. Default: 9.')
    group_blast.add_option('--blast-e-2', dest='blast_evalue_fg', default=1E-15, type=float,
                           help="Blastn evalue (-evalue) for fastg input. Default: 1E-15.")
    # group_blast.add_option('--blastn', dest='other_blastn_argv', default='',
    #                        help='Other blastn options. Use double quotation marks to include all the arguments '
    #                             'and parameters, such as "-evalue 1E-30 -word_size 11"')

    group_aftertreatment = OptionGroup(parser, "AFTERTREATMENT OPTIONS", "")
    group_aftertreatment.add_option('--raw', dest='mapping_like', default=False, action='store_true',
                                    help='Produce mapping like alignment for raw blast result.')
    group_aftertreatment.add_option('--aligned', dest='no_final_aligned', default=True, action='store_false',
                                    help='By default, this script produce aligned seqs, which would cost more '
                                         'memory and CPU. Choose to disable aligning for large-scale data. '
                                         'Anyway, you need to use mafft/muscle to get a more accurate alignment.')
    group_aftertreatment.add_option('--concatenate', dest='separate_result', default=True, action='store_false',
                                    help='By default, this script produce clusters of sequences separated by '
                                         'conservative sites. One file for each cluster. Choose to produce '
                                         'only one file of alignment/sequences.')
    group_aftertreatment.add_option('--gap-to-keep', dest='gap_to_keep', default=20,
                                    help='By default, this script would replace long gap with 20-N gap at last.')
    group_aftertreatment.add_option('--conserve-length', dest='conserve_len', default=40, type=int,
                                    help='Twice the minimum conservative site length for each cluster. Default: 40. ')

    group_filter = OptionGroup(parser, "FILTER OPTIONS", "")
    group_filter.add_option('--match-score', dest='match_score', default=1.0, type=float,
                            help='match score for evaluating best aligned block. Default: 1.0')
    group_filter.add_option('--mismatch-score', dest='mismatch_score', default=0.0, type=float,
                            help='mismatch penalty for evaluating best aligned block. Default: 0.0')
    group_filter.add_option('--gap-score', dest='gap_score', default=1.0, type=float,
                            help='gap open/extension penalty for evaluating best aligned block. Default: 1.0')
    group_filter.add_option('--h-cutoff', dest='hit_cut_off', default=0.75, type=float,
                            help='Discontinuous hit with overlap percent above which would be discarded.')
    group_filter.add_option('--q-cutoff', dest='query_cut_off', default=0.75, type=float,
                            help='Discontinuous hit with overlap percent above which would be discarded. '
                                 '(Beta, not applicable when multiple queries overlap)')

    group_others = OptionGroup(parser, "OTHER OPTIONS", "")
    group_others.add_option('--continue', dest='resume', default=False, action='store_true',
                            help="If blast result exist, skipped making blast.")
    group_others.add_option('--verbose', dest='verbose', default=False, action='store_true',
                            help='Verbose logging.')
    group_others.add_option('--keep-temp', dest='keep_temp', default=False, action='store_true',
                            help='Choose to keep all temp files.')
    # group_others.add_option('--keep-ref', dest='keep_ref', default=False, action='store_true',
    #                         help='Choose to keep modified reference.')
    parser.add_option_group(group_basic)
    parser.add_option_group(group_pretreatment)
    parser.add_option_group(group_blast)
    parser.add_option_group(group_filter)
    parser.add_option_group(group_extension)
    parser.add_option_group(group_aftertreatment)
    parser.add_option_group(group_others)
    options, args = parser.parse_args()

    if not (options.reference_fasta and len(args) and options.output_directory):
        parser.print_help()
        sys.stdout.write('\n######################################\nERROR: Insufficient REQUIRED arguments!\n\n')
        exit()
    if 0 < options.min_repeat < 50 or options.min_repeat < 0:
        parser.print_help()
        sys.stdout.write("\nIllegal minimum repeat length input!")
        exit()
    if options.keep_repeat_ends > int(options.min_repeat/2):
        sys.stdout.write("The value of keep_repeat_ends is required to be smaller than half minimum repeat.")
        exit()
    if not os.path.exists(options.output_directory):
        os.mkdir(options.output_directory)
    if not os.path.exists(os.path.join(options.output_directory, 'Blast')):
        os.mkdir(os.path.join(options.output_directory, 'Blast'))
    if options.mapping_like and not os.path.exists(os.path.join(options.output_directory, 'Raw')):
        os.mkdir(os.path.join(options.output_directory, 'Raw'))
    global keep_n_len
    keep_n_len = options.gap_to_keep
    log = simple_log(logging.getLogger(), options.output_directory)
    log.info(print_title)
    log.info(' '.join(sys.argv) + '\n')
    log = timed_log(log, options.output_directory)
    return options, args, log


def log_repeats(repeats_to_report, min_repeat_length, total_len, log):
    repeat_lengths = []
    if repeats_to_report:
        repeat_report = "Repeats (>=" + str(min_repeat_length) + "): "
        for repeat_item in repeats_to_report:
            x, y, z = repeat_item[0]
            if x == y:
                repeat_lengths.append(str(1))
            else:
                if (y - x)*z > 0:
                    repeat_lengths.append(str((y - x) * z + 1))
                else:
                    repeat_lengths.append(str((x, y)[z == 1] + total_len - (x, y)[z != 1] + 1))
        p_len = max([len(k) for k in repeat_lengths]) + 1
        count = 0
        for repeat_item in repeats_to_report:
            string_repeat_item = [(x + 1, y + 1, 'forward' if z == 1 else 'reverse') for x, y, z in sorted(repeat_item)]
            repeat_report += '\n' + 32 * ' ' + repeat_lengths[count].ljust(p_len) + ":" + str(tuple(string_repeat_item))
            count += 1
    else:
        repeat_report = "Repeats (>=" + str(min_repeat_length) + "): None"
    log.info(repeat_report)
    return repeat_lengths


def detect_repeats(sequence_string, min_repeat_length, circular, log,
                   seq_out=False, word_size=50, accepted_char=set(list("ATGCRMYKHBDVatgcrmykhbdv"))):
    log.info("Detecting repeats ...")
    word_size = min(word_size, min_repeat_length)
    if len(sequence_string) < min_repeat_length:
        if seq_out:
            return [None, None, None]
        else:
            return [None, None]
    if circular:
        long_sequence = sequence_string + sequence_string[:word_size - 1]
        here_seq = long_sequence
    else:
        here_seq = sequence_string
    here_seq_r = complementary_seq(here_seq)
    sequence_string_r = here_seq_r[word_size-1:]
    here_seq_length = len(here_seq)
    raw_seq_len = len(sequence_string)
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
    if circular:
        if len(repeat_indices) != raw_seq_len and len(repeat_indices):
            while (repeat_indices[0] - repeat_indices[-1]) % raw_seq_len == 1:
                repeat_indices.insert(0, repeat_indices.pop(-1))
        for i in range(len_indices):
            this_index = repeat_indices[i]
            this_word = index_to_words[this_index]
            this_connection = words_to_index[this_word]

            # test_id = 27753
            # if this_index == test_id:
            #     print(test_id)
            #     print('last connection', last_connection)
            #     print('this connection', this_connection)
            #     print('points to repeats', points_to_repeats)
            #     print(repeats)
            #     print('\n')

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
                candidate_new_connect = ((here_id + direction_trans) % raw_seq_len, direction_trans)
                if candidate_new_connect not in this_connection:
                    # if one_connection in points_to_repeats:
                    for repeat_kind, repeat_num in points_to_repeats[one_connection]:
                        if repeat_kind in repeats_to_stop:
                            repeats_to_stop[repeat_kind][repeat_num] = one_connection
                        else:
                            repeats_to_stop[repeat_kind] = {repeat_num: one_connection}
                        kinds_del_from_points.add(repeat_kind)
            # if this_index == test_id:
            #     print(test_id)
            #     print('repeats to stop', repeats_to_stop)
            #     print('points to repeats', points_to_repeats)
            #     print('kinds del from points', kinds_del_from_points)
            #     print()
            for repeat_kind in kinds_del_from_points:
                for now_start, now_go_to, now_direction in repeats[repeat_kind]:
                    connection_del_from_points = ((now_go_to - (word_size - 1) * (now_direction == 1)) % raw_seq_len,
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

            # if this_index == test_id:
            #     print('points to repeats2', points_to_repeats)
            #     print()
            for one_connection in last_connection:
                here_id, direction_trans = one_connection
                candidate_new_connect = ((here_id + direction_trans) % raw_seq_len, direction_trans)
                if candidate_new_connect in this_connection:
                    if one_connection in points_to_repeats:
                        points_to_repeats[candidate_new_connect] = []
                        for one_repeat_id in points_to_repeats[one_connection]:
                            points_to_repeats[candidate_new_connect].append(one_repeat_id)
                        del points_to_repeats[one_connection]
            # if this_index == test_id:
            #     print('points to repeats3', points_to_repeats)
            #     print()
            for repeat_kind in repeats_to_stop:
                if len(repeats[repeat_kind]) - len(repeats_to_stop[repeat_kind]) >= 2:
                    repeat_to_be_continued = False
                    for repeat_num in range(len(repeats[repeat_kind])):
                        if repeat_num not in repeats_to_stop[repeat_kind]:
                            start_id, go_to_id, gt_direction = repeats[repeat_kind][repeat_num]
                            new_connect = (
                                (go_to_id - (word_size - 1) * (gt_direction == 1) + gt_direction) % raw_seq_len,
                                gt_direction)
                            if new_connect in this_connection and new_connect not in points_to_repeats:
                                repeat_to_be_continued = True
                                break
                    if repeat_to_be_continued:
                        repeats.append([])
                        for inside_repeat_num in range(len(repeats[repeat_kind])):
                            if inside_repeat_num not in repeats_to_stop[repeat_kind]:
                                start_id, go_to_id, gt_direction = repeats[repeat_kind][inside_repeat_num]
                                repeats[-1].append([start_id, go_to_id, gt_direction])
                                new_connect = (
                                    (go_to_id - (word_size - 1) * (gt_direction == 1) + gt_direction) % raw_seq_len,
                                    gt_direction)
                                if new_connect in points_to_repeats:
                                    points_to_repeats[new_connect].append((len(repeats) - 1, len(repeats[-1]) - 1))
                                else:
                                    points_to_repeats[new_connect] = [(len(repeats) - 1, len(repeats[-1]) - 1)]
            for one_connection in points_to_repeats:
                for repeat_kind, repeat_num in points_to_repeats[one_connection]:
                    start_id, previous_id, this_direction = repeats[repeat_kind][repeat_num]
                    repeats[repeat_kind][repeat_num][1] += this_direction
                    repeats[repeat_kind][repeat_num][1] %= raw_seq_len
            # if this_index == test_id:
            #     print('points to repeats4', points_to_repeats)
            #     print('repeats', repeats)
            #     print()

            """part 2: dealing with new connection"""
            for one_connection in this_connection:
                here_id, direction_trans = one_connection
                candidate_last_connect = ((here_id - direction_trans) % raw_seq_len, direction_trans)
                if candidate_last_connect not in last_connection:
                    # if this_index == 27791:
                    repeats.append([])
                    # started_words.add(this_word)
                    for inside_connection in this_connection:
                        inside_id, inside_direction = inside_connection
                        repeats[-1].append([(inside_id + (word_size - 1) * (inside_direction == -1)) % raw_seq_len,
                                            (inside_id + (word_size - 1) * (inside_direction == 1)) % raw_seq_len,
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
            if next_index != (this_index + 1) % raw_seq_len:
                points_to_repeats = {}
                last_connection = set()
            else:
                last_connection = this_connection
    else:
        for i in range(len_indices):
            this_index = repeat_indices[i]
            this_word = index_to_words[this_index]
            this_connection = words_to_index[this_word]

            # test_id = 27753
            # if this_index == test_id:
            #     print(test_id)
            #     print('last connection', last_connection)
            #     print('this connection', this_connection)
            #     print('points to repeats', points_to_repeats)
            #     print(repeats)
            #     print('\n')

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
            # if this_index == test_id:
            #     print(test_id)
            #     print('repeats to stop', repeats_to_stop)
            #     print('points to repeats', points_to_repeats)
            #     print('kinds del from points', kinds_del_from_points)
            #     print()
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

            # if this_index == test_id:
            #     print('points to repeats2', points_to_repeats)
            #     print()
            for one_connection in last_connection:
                here_id, direction_trans = one_connection
                candidate_new_connect = (here_id + direction_trans, direction_trans)
                if candidate_new_connect in this_connection:
                    if one_connection in points_to_repeats:
                        points_to_repeats[candidate_new_connect] = []
                        for one_repeat_id in points_to_repeats[one_connection]:
                            points_to_repeats[candidate_new_connect].append(one_repeat_id)
                        del points_to_repeats[one_connection]
            # if this_index == test_id:
            #     print('points to repeats3', points_to_repeats)
            #     print()
            for repeat_kind in repeats_to_stop:
                if len(repeats[repeat_kind]) - len(repeats_to_stop[repeat_kind]) >= 2:
                    repeat_to_be_continued = False
                    for repeat_num in range(len(repeats[repeat_kind])):
                        if repeat_num not in repeats_to_stop[repeat_kind]:
                            start_id, go_to_id, gt_direction = repeats[repeat_kind][repeat_num]
                            new_connect = (go_to_id - (word_size - 1) * (gt_direction == 1) + gt_direction,
                                           gt_direction)
                            if new_connect in this_connection and new_connect not in points_to_repeats:
                                repeat_to_be_continued = True
                                break
                    if repeat_to_be_continued:
                        repeats.append([])
                        for inside_repeat_num in range(len(repeats[repeat_kind])):
                            if inside_repeat_num not in repeats_to_stop[repeat_kind]:
                                start_id, go_to_id, gt_direction = repeats[repeat_kind][inside_repeat_num]
                                repeats[-1].append([start_id, go_to_id, gt_direction])
                                new_connect = (go_to_id - (word_size - 1) * (gt_direction == 1) + gt_direction,
                                               gt_direction)
                                if new_connect in points_to_repeats:
                                    points_to_repeats[new_connect].append((len(repeats) - 1, len(repeats[-1]) - 1))
                                else:
                                    points_to_repeats[new_connect] = [(len(repeats) - 1, len(repeats[-1]) - 1)]
            for one_connection in points_to_repeats:
                for repeat_kind, repeat_num in points_to_repeats[one_connection]:
                    start_id, previous_id, this_direction = repeats[repeat_kind][repeat_num]
                    repeats[repeat_kind][repeat_num][1] += this_direction
            # if this_index == test_id:
            #     print('points to repeats4', points_to_repeats)
            #     print('repeats', repeats)
            #     print()

            """part 2: dealing with new connection"""
            for one_connection in this_connection:
                here_id, direction_trans = one_connection
                candidate_last_connect = (here_id - direction_trans, direction_trans)
                if candidate_last_connect not in last_connection:
                        # if this_index == 27791:
                        repeats.append([])
                        # started_words.add(this_word)
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

    log_repeats(final_repeat, word_size, raw_seq_len, log)
    # delete small repeats
    count_group__ = 0
    while count_group__ < len(final_repeat):
        here_start, here_end, here_direction = final_repeat[count_group__][0]
        if ((here_end - here_start + here_direction)*here_direction) % raw_seq_len < min_repeat_length:
            del final_repeat[count_group__]
        else:
            count_group__ += 1
    # out seq
    for group_to_sort in range(len(final_repeat)):
        final_repeat[group_to_sort].sort(key=lambda x:min(x[:2]))
    repeat_lengths = log_repeats(final_repeat, min_repeat_length, raw_seq_len, log)
    if seq_out:
        final_seqs = []
        for count_group___ in range(len(final_repeat)):
            final_seqs.append([])
            for count_num_ in range(len(final_repeat[count_group___])):
                here_start, here_end, here_direction = final_repeat[count_group___][count_num_]
                if here_direction == 1:
                    if here_end + 1 > here_start:
                        final_seqs[-1].append(sequence_string[here_start: here_end + 1])
                    else:
                        final_seqs[-1].append(sequence_string[here_start:] + sequence_string[:here_end + 1])
                else:
                    if here_end < here_start + 1:
                        final_seqs[-1].append(
                            sequence_string_r[raw_seq_len - here_start - 1: raw_seq_len - here_end])
                    else:
                        final_seqs[-1].append(sequence_string_r[raw_seq_len - here_start - 1:] +
                                              sequence_string_r[:raw_seq_len - here_end])
        log.info("Detecting repeats finished.")
        return final_repeat, repeat_lengths, final_seqs
    else:
        # print(final_repeat)
        log.info("Detecting repeats finished.")
        return final_repeat, repeat_lengths


def remove_repeats(seq_to_modify, repeats_info, keep_repeat_ends, del_randomly, circular, add_gap, log):
    # to improve:
    # 1. do not consider those without any overlap with other repeats
    # 2. randomly delete
    log.info("Removing repeats ...")
    repeats, repeat_lengths = repeats_info
    seq_len = len(seq_to_modify)
    candidates = set(range(seq_len))

    def generate_site_from_repeats(this_repeat):
        here_start, here_end, here_direction = this_repeat
        if here_direction == 1:
            if here_end + 1 > here_start:
                for in_site in range(here_start, here_end + 1):
                    yield in_site
            else:
                for in_site in range(here_start, seq_len):
                    yield in_site
                for in_site in range(0, here_end + 1):
                    yield in_site
        else:
            if here_end < here_start + 1:
                for in_site in range(here_end, here_start + 1):
                    yield in_site
            else:
                for in_site in range(here_end, seq_len):
                    yield in_site
                for in_site in range(0, here_start + 1):
                    yield in_site

    if del_randomly:
        # I'm tired to write
        pass
    else:
        # # choose the global minimum length to keep
        # find repeat groups that overlap
        taken_sites = {}
        overlapping_repeats = []
        recorded_overlap = set()
        for count_g in range(len(repeats)):
            for candidate_single_repeat in repeats[count_g]:
                find_overlap = False
                for coming_site in generate_site_from_repeats(candidate_single_repeat):
                    if coming_site in taken_sites:
                        find_overlap = True
                        taking_group_num = taken_sites[coming_site]
                        if taking_group_num not in recorded_overlap:
                            overlapping_repeats.append(sorted(repeats[taking_group_num], key=lambda x: min(x[:2])))
                            recorded_overlap.add(taking_group_num)
                        if count_g not in recorded_overlap:
                            overlapping_repeats.append(sorted(repeats[count_g], key=lambda x: min(x[:2])))
                            recorded_overlap.add(count_g)
                        break
                    else:
                        taken_sites[coming_site] = count_g
                if find_overlap:
                    break
        # choose the first repeat from those groups that do not overlap
        for count_g in range(len(repeats)):
            if count_g not in recorded_overlap:
                for repeat_to_del in repeats[count_g][1:]:
                    for site_to_del in generate_site_from_repeats(repeat_to_del):
                        candidates.discard(site_to_del)

        def create_recombination_list(previous_list, item_length_list):
            if item_length_list:
                if previous_list:
                    new_list = []
                    for j_ in range(item_length_list.pop(0)):
                        for k_ in range(len(previous_list)):
                            new_list.append(previous_list[k_]+[j_])
                    return create_recombination_list(new_list, item_length_list)
                else:
                    return create_recombination_list([[j_] for j_ in range(item_length_list.pop(0))], item_length_list)
            else:
                return previous_list
        combinations = create_recombination_list([], [len(repeat_group) for repeat_group in overlapping_repeats])
        best_choice = [seq_len, set(), ()]
        count_choice = 1
        total_choice = len(combinations)
        for choice in combinations:
            this_waive = set()
            for count_g_ in range(len(overlapping_repeats)):
                for this_site in generate_site_from_repeats(overlapping_repeats[count_g_][choice[count_g_]]):
                    this_waive.add(this_site)
            len_this_waive = len(this_waive)
            if len_this_waive < best_choice[0]:
                best_choice = [len_this_waive, this_waive, choice]
            to_print = str(count_choice)+'/'+str(total_choice)+':'+str(choice)
            sys.stdout.write(to_print+'\b'*len(to_print))
            sys.stdout.flush()
            count_choice += 1
        if total_choice:
            log.info("Best removing choice: "+str(best_choice[2])+" of "+str(total_choice)+" choices")
        for count_g_2 in range(len(overlapping_repeats)):
            for count_n_ in range(len(overlapping_repeats[count_g_2])):
                for this_site in generate_site_from_repeats(overlapping_repeats[count_g_2][count_n_]):
                    if this_site not in best_choice[1]:
                        candidates.discard(this_site)
    resumed = True
    removing_ranges = []
    for here_site in range(seq_len):
        if here_site in candidates:
            if not resumed:
                removing_ranges[-1].append(here_site)
                resumed = True
        elif resumed:
            removing_ranges.append([here_site+1])
            resumed = False
    if not resumed:
        removing_ranges[-1].append(here_site+1)

    if keep_repeat_ends:
        removing_report = 'Real Repeat Range:'
        for removing_item in removing_ranges:
            removing_report += '\n' + 32 * ' ' + str(tuple(removing_item))
        log.info(removing_report)
        if circular:
            if len(removing_ranges) == 1:
                removing_ranges[0][0] += keep_repeat_ends
                removing_ranges[0][1] -= keep_repeat_ends
            elif int(removing_ranges[0][0] == 1) + int(removing_ranges[-1][-1] == seq_len) <= 1:
                for k in range(len(removing_ranges)):
                    removing_ranges[k][0] += keep_repeat_ends
                    removing_ranges[k][1] -= keep_repeat_ends
            else:
                for k in range(1, len(removing_ranges)-1):
                    removing_ranges[k][0] += keep_repeat_ends
                    removing_ranges[k][1] -= keep_repeat_ends
                removing_ranges[0][1] -= keep_repeat_ends
                if removing_ranges[0][1] < 1:
                    removing_ranges[-1][1] += removing_ranges[0][1]
                    del removing_ranges[0]
                removing_ranges[-1][0] += keep_repeat_ends
                if removing_ranges[-1][0] > seq_len:
                    removing_ranges[0][0] += removing_ranges[-1][0] - seq_len
                    del removing_ranges[-1]
        else:
            if len(removing_ranges) == 1:
                if removing_ranges[0][0] != 1:
                    removing_ranges[0][0] += keep_repeat_ends
                if removing_ranges[0][1] != seq_len:
                    removing_ranges[0][1] -= keep_repeat_ends
            else:
                if removing_ranges[0][0] != 1:
                    removing_ranges[0][0] += keep_repeat_ends
                removing_ranges[0][1] -= keep_repeat_ends
                if removing_ranges[-1][1] != seq_len:
                    removing_ranges[-1][1] -= keep_repeat_ends
                removing_ranges[-1][0] += keep_repeat_ends
                for k in range(1, len(removing_ranges)-1):
                    removing_ranges[k][0] += keep_repeat_ends
                    removing_ranges[k][1] -= keep_repeat_ends

    # create new seq
    out_seq = seq_to_modify[:removing_ranges[0][0]-1]
    for k in range(1, len(removing_ranges)):
        out_seq += "N"*add_gap
        out_seq += seq_to_modify[removing_ranges[k-1][1]:removing_ranges[k][0]-1]
    out_seq += "N"*add_gap
    if circular:
        out_seq = seq_to_modify[removing_ranges[-1][1]:] + out_seq
    else:
        out_seq += seq_to_modify[removing_ranges[-1][1]:]
    # log removing
    removing_report = 'Removing:'
    for removing_item in removing_ranges:
        removing_report += '\n'+32*' '+str(removing_item[1]-removing_item[0]+1)+': '+str(tuple(removing_item))
    log.info(removing_report)
    log.info("Removing repeats finished.")
    return out_seq


def check_db(reference_fa_base, min_repeat, word_size, keep_repeat_ends, del_randomly, circular, add_gap, log):
    log.info('Making BLAST db ... ')
    if os.path.isfile(reference_fa_base):
        ref_fasta = read_fasta_gb_head(reference_fa_base)
        if len(ref_fasta[0]) > 1:
            """Removing repeats from reference"""
            if min_repeat:
                repeats = detect_repeats(ref_fasta[1][0], min_repeat, circular, log, word_size=word_size)
                if repeats[0]:
                    this_seq = remove_repeats(ref_fasta[1][0], repeats, keep_repeat_ends, del_randomly, circular,
                                              add_gap, log)
                else:
                    this_seq = ref_fasta[1][0]
            reference_fa_base += '.modified'
            write_fasta(out_dir=reference_fa_base,
                        matrix=[[ref_fasta[0][0]], [this_seq], ref_fasta[2]], overwrite=True)
            log.warning('multi-seqs in reference file, only use the 1st sequence.')
        elif len(ref_fasta[0]) == 0:
            log.error('illegal reference file!')
            exit()
        else:
            """Removing repeats from reference"""
            if min_repeat:
                repeats = detect_repeats(ref_fasta[1][0], min_repeat, circular, log, word_size=word_size)
                if repeats[0]:
                    this_seq = remove_repeats(ref_fasta[1][0], repeats, keep_repeat_ends, del_randomly, circular,
                                              add_gap, log)
                    reference_fa_base += '.modified'
                    write_fasta(out_dir=reference_fa_base,
                                matrix=[[ref_fasta[0][0]], [this_seq], ref_fasta[2]], overwrite=True)

        try:
            # python2
            makedb_result = subprocess.getstatusoutput('makeblastdb -dbtype nucl -in '+reference_fa_base +
                                                       ' -out '+reference_fa_base+'.index')
        except AttributeError:
            # python3
            makedb_result = commands.getstatusoutput('makeblastdb -dbtype nucl -in ' + reference_fa_base +
                                                     ' -out ' + reference_fa_base + '.index')
        if 'Error' in str(makedb_result[1]) or 'error' in str(makedb_result[1]) or '不是内部或外部命令' in str(makedb_result[1]):
            if not os.path.exists(reference_fa_base+'.index.nhr'):
                log.error('Blast terminated with following info:\n'+str(makedb_result[1]))
                exit()
        in_index = reference_fa_base+'.index'
        log.info('Making BLAST db finished.\n')
    else:
        log.error('No illegal reference input!')
        exit()
    return in_index


def execute_blast(query, blast_db, output, outfmt, b_word_size, b_evalue, resume, log):
    log.info("Executing BLAST ...")
    if resume and os.path.exists(output):
        log.info("BLAST result existed. Skipped.")
    else:
        this_command = 'blastn -num_threads 4 -query ' + query + ' -db ' + blast_db + ' -out ' + output + ' -outfmt ' \
                       + str(outfmt) + ' -word_size ' + str(b_word_size) + ' -evalue ' + str(b_evalue)
        make_blast = subprocess.Popen(this_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        blast_out, blast_error = make_blast.communicate()
        if 'error' in str(blast_out) or "Error" in str(blast_out) or "command not found" in str(blast_out):
            log.error("BLAST error")
            exit()
        log.info("Executing BLAST finished.")


def parse_blast_xml_result(xml_file, log):
    log.info("Parsing BLAST result ...")
    tree = ET.parse(xml_file)
    root = tree.getroot()
    hsp_dict = {}
    for grandchild in root.find('BlastOutput_iterations'):
        # hit one reference
        hits = grandchild.find('Iteration_hits')
        if len(hits) == 1:
            this_name = grandchild.find('Iteration_query-def').text
            temp_dict = {}
            for hsp in hits.find('Hit').find('Hit_hsps'):
                temp_sub_dict = {}
                for hsp_character in hsp:
                    temp_sub_dict[hsp_character.tag] = hsp_character.text
                temp_dict[temp_sub_dict['Hsp_num']] = temp_sub_dict
            re_order_list = sorted(list(temp_dict), key=lambda x: (int(temp_dict[x]['Hsp_query-from']),
                                                                   int(temp_dict[x]['Hsp_query-to'])))
            for re_ordered_num in range(len(re_order_list)):
                hsp_dict[this_name+'--Hsp_num_'+str(re_ordered_num+1)] = temp_dict[re_order_list[re_ordered_num]]

    log.info("Parsing BLAST result finished.")
    return hsp_dict
""" hsp_dict = {query1: {'Hsp_num': '5',
                         'Hsp_bit-score': '364.91',
                         'Hsp_score': '197',
                         'Hsp_evalue': '1.70367e-100',
                         'Hsp_query-from': '10141',
                         'Hsp_query-to': '10621',
                         'Hsp_hit-from': '398089'
                         'Hsp_hit-to': '397606',
                         'Hsp_query-frame': '1'
                         'Hsp_hit-frame': '-1',
                         'Hsp_identity': '392',
                         'Hsp_positive': '392',
                         'Hsp_align-len': '487',
                         'Hsp_gaps': '9',
                         'Hsp_qseq': 'TCACCTTCTCTCCCCCCCCCGTTCCCACACTTTTCCCAATTCCCACCTCTTTATTTTTTTAACTTTTTTTAAAA',
                         'Hsp_hseq': 'TCACCTTCTCCCCCCCCACTGTTCCCACTGTCTTCCTAATTCCCACCTCTTGATTTTTTTGACTTTCTTTAAAA',
                         'Hsp_midline': '|||||||||| |||||| | ||||||||  | |||| |||||||||||||| |||||||| ||||| |||||||'},
                 query2: {...} ...}"""


def initialize_site_dict(dict_length):
    site_dicts = {site: {} for site in range(1, dict_length + 1)}
    for i in range(1, dict_length):
        site_dicts[(i, i + 1)] = []
    site_dicts[(i + 1, 1)] = []
    return site_dicts


def hsp_hits_to_hit_site_dicts(custom_hsp_dict, total_len, log):
    log.info("Hsp information to site dicts ... ")
    # initialization
    query_range_dicts = {}
    hit_site_dicts = initialize_site_dict(total_len)

    for query_name, hsp_hit in custom_hsp_dict.items():
        query_seq = hsp_hit['Hsp_qseq']
        query_from = int(hsp_hit['Hsp_query-from'])
        query_to = int(hsp_hit['Hsp_query-to'])
        hit_seq = hsp_hit['Hsp_hseq']
        hit_from = int(hsp_hit['Hsp_hit-from'])
        hit_to = int(hsp_hit['Hsp_hit-to'])
        hit_direction = hsp_hit['Hsp_hit-frame'] == '1'
        hsp_score = int(hsp_hit['Hsp_score'])
        if not hit_direction:
            hit_from, hit_to = hit_to, hit_from
            query_from, query_to = query_to, query_from
            query_seq = complementary_seq(query_seq)
            hit_seq = complementary_seq(hit_seq)
        query_seq = list(query_seq)
        hit_seq = list(hit_seq)
        query_range_dicts[query_name] = {'h': [hit_from, hit_to], 'q': [query_from, query_to],
                                         'dire': hit_direction, 'score': hsp_score}
        count_site = int(hit_from)
        query_base = query_seq.pop(0)
        hit_base = hit_seq.pop(0)
        try:
            while True:
                if hit_base != '-':
                    hit_site_dicts[count_site][query_name] = [query_base, int(((query_base == hit_base)-0.5)*2)]
                    count_site += 1
                    query_base = query_seq.pop(0)
                    hit_base = hit_seq.pop(0)
                else:
                    count_gap = 0
                    gap_position = (count_site - 1, count_site)
                    while query_seq and hit_base == '-':
                        if count_gap < len(hit_site_dicts[gap_position]):
                            hit_site_dicts[gap_position][count_gap][query_name] = [query_base, -1]
                        else:
                            hit_site_dicts[gap_position].append({query_name: [query_base, -1]})
                        query_base = query_seq.pop(0)
                        hit_base = hit_seq.pop(0)
                        count_gap += 1
        except IndexError:
            pass
    log.info("Hsp information to site dicts finished.")
    return hit_site_dicts, query_range_dicts
""" hit_site_dicts = {1: {name1: ['A', 1], name2:['A', 1]}],
                     (1, 2): [{name1:['A', -1]}, {name1: ['A', -1]}],
                      2: {name1: ['T', -1], name2:['A', 1]}],
                      ...} """


def delete_in_middle(hit_site_dicts, q_range_dicts, name_to_del_in_middle,
                     overlap_start, overlap_hit_len, overlap_query_len, direction_trans):
    count_previous_q_len = 0
    for this_hit_site in range(q_range_dicts[name_to_del_in_middle]['h'][0], overlap_start):
        count_previous_q_len += hit_site_dicts[this_hit_site][name_to_del_in_middle][0] != '-'
        # no need to worry about circular here, since this_hit_site < go_to <= total_len
        for here_gap_site_dict in hit_site_dicts[(this_hit_site, this_hit_site + 1)]:
            count_previous_q_len += name_to_del_in_middle in here_gap_site_dict
    # make a new contig
    new_name = name_to_del_in_middle + "-new"
    q_range_dicts[new_name] = {
        'h': [overlap_start + overlap_hit_len, q_range_dicts[name_to_del_in_middle]['h'][1]],
        'q': [q_range_dicts[name_to_del_in_middle]['q'][0] + direction_trans * (
            count_previous_q_len + overlap_query_len), q_range_dicts[name_to_del_in_middle]['q'][1]],
        'dire': q_range_dicts[name_to_del_in_middle]['dire'],
        'score': q_range_dicts[name_to_del_in_middle]['score']}
    # shorten the previous
    q_range_dicts[name_to_del_in_middle]['q'][1] = \
        q_range_dicts[name_to_del_in_middle]['q'][0] + direction_trans * (count_previous_q_len - 1)
    q_range_dicts[name_to_del_in_middle]['h'][1] = overlap_start - 1
    # replace the hit dicts
    for this_hit_site in range(overlap_start + overlap_hit_len, q_range_dicts[new_name]['h'][1] + 1):
        if this_hit_site > 1:
            here_gap_site = (this_hit_site - 1, this_hit_site)
            for gap_dict_index in range(len(hit_site_dicts[here_gap_site])):
                if name_to_del_in_middle in hit_site_dicts[here_gap_site][gap_dict_index]:
                    hit_site_dicts[here_gap_site][gap_dict_index][new_name] = \
                        hit_site_dicts[here_gap_site][gap_dict_index][name_to_del_in_middle]
                    del hit_site_dicts[here_gap_site][gap_dict_index][name_to_del_in_middle]
        hit_site_dicts[this_hit_site][new_name] = hit_site_dicts[this_hit_site][name_to_del_in_middle]
        del hit_site_dicts[this_hit_site][name_to_del_in_middle]


def del_gap_of_name(hit_site_dicts, name_to_del, gap_site):
    count_gap = 0
    while count_gap < len(hit_site_dicts[gap_site]):
        if name_to_del in hit_site_dicts[gap_site][count_gap]:
            del hit_site_dicts[gap_site][count_gap][name_to_del]
            if not hit_site_dicts[gap_site][count_gap]:
                del hit_site_dicts[gap_site][count_gap]
            else:
                count_gap += 1
        else:
            count_gap += 1


def delete_whole_hit_of_name(hit_site_dicts, q_range_dics, name_cluster, short_name, num_of_name):
    name_to_del = name_cluster[short_name][num_of_name]
    start_h, end_h = q_range_dics[name_to_del]['h']
    previous_gap = (start_h - 1, start_h)
    if previous_gap in hit_site_dicts:
        del_gap_of_name(hit_site_dicts, name_to_del, previous_gap)
    for this_site in range(start_h, end_h + 1):
        if name_to_del in hit_site_dicts[this_site]:
            del hit_site_dicts[this_site][name_to_del]
        this_gap = (this_site, this_site + 1)
        if this_gap in hit_site_dicts:
            del_gap_of_name(hit_site_dicts, name_to_del, this_gap)
    del q_range_dics[name_to_del]
    del name_cluster[short_name][num_of_name]


def strip_new(old_string):
    while old_string.endswith('-new'):
        old_string = old_string[:-4]
    return old_string


def detect_continuity(this_name, previous_name_id, next_name_id, this_name_cluster,
                      q_range_dicts, hit_site_dicts, total_len, circular_h,
                      extra_prev_name=set(), extra_next_name=set()):
    continuity = True
    this_name_info = q_range_dicts[this_name]
    this_direction = this_name_info['dire']
    this_direction_trans = int((this_direction - 0.5) * 2)
    # part one, test previous query
    if previous_name_id >= 0:
        supposed_prev_name = strip_new(this_name_cluster[previous_name_id])
        if circular_h:
            go_to_h_for_prev = (this_name_info['h'][int(not this_direction)]
                                - this_direction_trans - 1) % total_len + 1
            while not hit_site_dicts[go_to_h_for_prev]:
                go_to_h_for_prev = (go_to_h_for_prev - 1 - this_direction_trans) % total_len + 1
        else:
            go_to_h_for_prev = this_name_info['h'][int(not this_direction)] - this_direction_trans
            while go_to_h_for_prev in hit_site_dicts and not hit_site_dicts[go_to_h_for_prev]:
                go_to_h_for_prev -= this_direction_trans
        if go_to_h_for_prev in hit_site_dicts:
            detected_prev_names = set([strip_new(pre_name)
                                       for pre_name in hit_site_dicts[go_to_h_for_prev]])
        for extra_name in extra_prev_name:
            detected_prev_names.add(strip_new(extra_name))
        if supposed_prev_name not in detected_prev_names:
            continuity = False
    # part two, test next query
    if next_name_id < len(this_name_cluster) and continuity:
        supposed_next_name = strip_new(this_name_cluster[next_name_id])
        if circular_h:
            go_to_h_for_next = (this_name_info['h'][int(this_direction)]
                                + this_direction_trans - 1) % total_len + 1
        else:
            go_to_h_for_next = this_name_info['h'][int(this_direction)] + this_direction_trans
            while go_to_h_for_next in hit_site_dicts and not hit_site_dicts[go_to_h_for_next]:
                go_to_h_for_next += this_direction_trans
        if go_to_h_for_next in hit_site_dicts:
            detected_next_names = set([strip_new(nex_name)
                                       for nex_name in hit_site_dicts[go_to_h_for_next]])
        for extra_name in extra_next_name:
            detected_next_names.add(strip_new(extra_name))
        if supposed_next_name not in detected_next_names:
            continuity = False
    return continuity


def remove_multiple_hits_per_query(hit_site_dicts, q_range_dics, total_len, circular_h,
                                   cut_off, match_s, mismatch_s, gap_s, name_cluster, log):
    log.info("Removing redundant hits per query ...")
    for this_q_name in q_range_dics:
        this_short_name = this_q_name.split('--')[0]
        if this_short_name in name_cluster:
            name_cluster[this_short_name].append(this_q_name)
        else:
            name_cluster[this_short_name] = [this_q_name]
    for short_name in name_cluster:
        range_modified = True
        while len(name_cluster[short_name]) > 1 and range_modified:
            name_cluster[short_name].sort(key=lambda x: (min(q_range_dics[x]['q']), min(q_range_dics[x]['h'])))
            i = 0
            while i + 1 < len(name_cluster[short_name]):
                this_name = name_cluster[short_name][i]
                next_name = name_cluster[short_name][i+1]
                this_name_info = q_range_dics[this_name]
                next_name_info = q_range_dics[next_name]
                this_direction = this_name_info['dire']
                next_direction = next_name_info['dire']
                next_direction_trans = int((next_direction - 0.5) * 2)
                this_direction_trans = int((this_direction - 0.5) * 2)
                this_1, this_2 = this_name_info['q']
                this_q_len = abs(this_2 - this_1) + 1
                next_1, next_2 = next_name_info['q']
                next_q_len = abs(next_2 - next_1) + 1
                q_overlap = min(max(this_name_info['q']) - min(next_name_info['q']) + 1, next_q_len)
                if q_overlap > 0:
                    """test whether sequential queries has sequential hits"""
                    # continuity test for this
                    continuity_1 = detect_continuity(this_name, i - 1, i + 2, name_cluster[short_name], q_range_dics,
                                                     hit_site_dicts, total_len, circular_h, extra_next_name={next_name})
                    # continuity test for next
                    continuity_2 = detect_continuity(next_name, i - 1, i + 2, name_cluster[short_name], q_range_dics,
                                                     hit_site_dicts, total_len, circular_h, extra_prev_name={this_name})
                    """judge according to continuity"""
                    if continuity_1 > continuity_2 \
                            and (q_overlap/float(next_q_len)) >= cut_off:
                        delete_whole_hit_of_name(hit_site_dicts, q_range_dics, name_cluster, short_name, i+1)
                        # print(this_name, next_name, 'case 1')
                        range_modified = True
                        break
                    elif continuity_2 > continuity_1 \
                            and (q_overlap/float(this_q_len)) >= cut_off:
                        delete_whole_hit_of_name(hit_site_dicts, q_range_dics, name_cluster, short_name, i)
                        # print(this_name, next_name, 'case 2')
                        range_modified = True
                        break
                    else:
                        # print(this_name, next_name, 'case 3')
                        # let the similarity (score) to judge
                        # keep the uniformity between searching best hit for a query and best query for a hit
                        # calculate the score for this
                        count_q_del = 0
                        count_h_del = 0
                        site_to_del_1 = q_range_dics[this_name]['h'][this_direction]
                        score_1 = 0
                        while count_q_del < q_overlap:
                            if hit_site_dicts[site_to_del_1 + count_h_del][this_name][0] != '-':
                                count_q_del += 1
                                score_1 += match_s
                            else:
                                score_1 -= mismatch_s
                            if this_direction:
                                this_gap = ((site_to_del_1 + count_h_del - 2) % total_len + 1,
                                            (site_to_del_1 + count_h_del - 1) % total_len + 1)
                            else:
                                this_gap = ((site_to_del_1 + count_h_del - 1) % total_len + 1,
                                            (site_to_del_1 + count_h_del) % total_len + 1)
                            for this_gap_dict in hit_site_dicts[this_gap]:
                                if this_name in this_gap_dict:
                                    count_q_del += 1
                                    score_1 -= gap_s
                                    if count_q_del == q_overlap:
                                        break
                            count_h_del -= this_direction_trans
                        # calculate score for next
                        count_q_del = 0
                        count_h_del = 0
                        site_to_del_2 = q_range_dics[next_name]['h'][not next_direction]
                        score_2 = 0
                        while count_q_del < q_overlap:
                            if hit_site_dicts[site_to_del_2+count_h_del][next_name][0] != '-':
                                count_q_del += 1
                                score_2 += match_s
                            else:
                                score_2 -= mismatch_s
                            if next_direction:
                                this_gap = ((site_to_del_2 + count_h_del - 1) % total_len + 1,
                                            (site_to_del_2 + count_h_del) % total_len + 1)
                            else:
                                this_gap = ((site_to_del_2 + count_h_del - 2) % total_len + 1,
                                            (site_to_del_2 + count_h_del - 1) % total_len + 1)
                            for this_gap_dict in hit_site_dicts[this_gap]:
                                if next_name in this_gap_dict:
                                    count_q_del += 1
                                    score_2 -= gap_s
                                    if count_q_del == q_overlap:
                                        break
                            count_h_del += next_direction_trans

                        count_q_del = 0
                        count_h_del = 0
                        # if this has better or equal score than next
                        if score_1 > score_2 or (score_1 == score_2 and continuity_1 >= continuity_2):
                            # del next
                            site_to_del_2 = q_range_dics[next_name]['h'][not next_direction]
                            while count_q_del < q_overlap:
                                if hit_site_dicts[site_to_del_2+count_h_del][next_name][0] != '-':
                                    count_q_del += 1
                                del hit_site_dicts[site_to_del_2+count_h_del][next_name]
                                if next_direction:
                                    this_gap = ((site_to_del_2 + count_h_del - 1) % total_len + 1,
                                                (site_to_del_2 + count_h_del) % total_len + 1)
                                    count_gap = 0
                                    while count_gap < len(hit_site_dicts[this_gap]):
                                        if next_name in hit_site_dicts[this_gap][count_gap]:
                                            del hit_site_dicts[this_gap][count_gap][next_name]
                                            count_q_del += 1
                                            if not hit_site_dicts[this_gap][count_gap]:
                                                del hit_site_dicts[this_gap][count_gap]
                                            else:
                                                count_gap += 1
                                        else:
                                            count_gap += 1
                                else:
                                    this_gap = ((site_to_del_2 + count_h_del - 2) % total_len + 1,
                                                (site_to_del_2 + count_h_del - 1) % total_len + 1)
                                    count_gap = len(hit_site_dicts[this_gap]) - 1
                                    while count_gap >= 0:
                                        if next_name in hit_site_dicts[this_gap][count_gap]:
                                            del hit_site_dicts[this_gap][count_gap][next_name]
                                            count_q_del += 1
                                            if not hit_site_dicts[this_gap][count_gap]:
                                                del hit_site_dicts[this_gap][count_gap]
                                        count_gap -= 1
                                count_h_del += next_direction_trans
                            q_range_dics[next_name]['q'][not next_direction] += count_q_del
                            q_range_dics[next_name]['h'][not next_direction] += count_h_del
                            if (q_range_dics[next_name]['q'][0] != q_range_dics[next_name]['q'][1]) and \
                                    (q_range_dics[next_name]['q'][0] > q_range_dics[next_name]['q'][1]) \
                                    == q_range_dics[next_name]['dire']:
                                del q_range_dics[next_name]
                                del name_cluster[short_name][i+1]
                        # if next has better score and next is within this
                        elif this_name_info['q'][this_direction]  > next_name_info['q'][next_direction]:
                            # del this
                            if score_1*this_q_len >= score_2*next_q_len:
                                delete_whole_hit_of_name(hit_site_dicts, q_range_dics, name_cluster, short_name, i + 1)
                            else:
                                delete_whole_hit_of_name(hit_site_dicts, q_range_dics, name_cluster, short_name, i)
                        # if next has better score
                        else:
                            site_to_del_1 = this_name_info['h'][this_direction]
                            while count_q_del < q_overlap:
                                if hit_site_dicts[site_to_del_1 + count_h_del][this_name][0] != '-':
                                    count_q_del += 1
                                del hit_site_dicts[site_to_del_1 + count_h_del][this_name]
                                if this_direction:
                                    this_gap = ((site_to_del_1 + count_h_del - 2) % total_len + 1,
                                                (site_to_del_1 + count_h_del - 1) % total_len + 1)
                                    count_gap = len(hit_site_dicts[this_gap]) - 1
                                    while count_gap >= 0:
                                        if this_name in hit_site_dicts[this_gap][count_gap]:
                                            del hit_site_dicts[this_gap][count_gap][this_name]
                                            count_q_del += 1
                                            if not hit_site_dicts[this_gap][count_gap]:
                                                del hit_site_dicts[this_gap][count_gap]
                                        count_gap -= 1
                                else:
                                    this_gap = ((site_to_del_1 + count_h_del - 1) % total_len + 1,
                                                (site_to_del_1 + count_h_del) % total_len + 1)
                                    count_gap = 0
                                    while count_gap < len(hit_site_dicts[this_gap]):
                                        if this_name in hit_site_dicts[this_gap][count_gap]:
                                            del hit_site_dicts[this_gap][count_gap][this_name]
                                            count_q_del += 1
                                            if not hit_site_dicts[this_gap][count_gap]:
                                                del hit_site_dicts[this_gap][count_gap]
                                            else:
                                                count_gap += 1
                                        else:
                                            count_gap += 1
                                count_h_del -= this_direction_trans
                            # change q_range_dicts
                            q_range_dics[this_name]['q'][this_direction] -= count_q_del
                            q_range_dics[this_name]['h'][this_direction] += count_h_del
                            # check validity of q_range_dicts
                            if (q_range_dics[this_name]['q'][0] != q_range_dics[this_name]['q'][1]) and \
                                    (q_range_dics[this_name]['q'][0] > q_range_dics[this_name]['q'][1]) \
                                    == q_range_dics[this_name]['dire']:
                                del q_range_dics[this_name]
                                del name_cluster[short_name][i]
                        range_modified = True
                        break
                range_modified = False
                i += 1
    log.info("Removing redundant hits per query finished.")


def update_to_cluster(name_cluster):
    to_cluster_info = {}
    for short_name in name_cluster:
        for name_id in range(len(name_cluster[short_name])):
            to_cluster_info[name_cluster[short_name][name_id]] = (short_name, name_id)
    return to_cluster_info


def sort_overlaps(overlap_value, overlap_count, to_cluster_info, name_cluster, q_range_dicts, hit_site_dicts,
                  total_len, circular_h):
    # sort by passing value_in_overlap*coverage+log(length)
    ave_depth = sum([abs(x[1][1])*x[1][3] for x in list(overlap_value.items())])\
                / sum([abs(x[1][1]) for x in list(overlap_value.items())])
    overlaps = sorted(list(overlap_value.items()),
                      key=lambda x: -(x[1][1] * x[1][3] / ave_depth + math.log(x[1][5]) + x[1][6]/overlap_count))
    # if "FC280-short-ndhB--Hsp_num_65" in overlap_value:
    #     print(overlaps)
    if overlaps[0][1][1]*overlaps[0][1][3]/ave_depth + math.log(overlaps[0][1][5])+overlaps[0][1][6]/overlap_count == \
            overlaps[1][1][1]*overlaps[1][1][3]/ave_depth + math.log(overlaps[1][1][5])+overlaps[1][1][6]/overlap_count:
        # sort by passing value_in_overlap*coverage+log(length),
        # continuity, first appear in query (for first IR)==name
        for this_name in overlap_value:
            this_short_name, this_name_id = to_cluster_info[strip_new(this_name)]
            overlap_value[this_name][4] = detect_continuity(this_name, this_name_id - 1, this_name_id + 1,
                                                            name_cluster[this_short_name],
                                                            q_range_dicts, hit_site_dicts,
                                                            total_len, circular_h)
        overlaps = sorted(list(overlap_value.items()),
                          key=lambda x: (-(x[1][1] * x[1][3]/ave_depth + math.log(x[1][5]) + x[1][6]/overlap_count),
                                         -x[1][4], x[0]))
    return overlaps


def initialize_overlap_value(overlap_names, overlap_value, overlap_count, go_to, hit_site_dicts, q_range_dicts,
                             edge_connections, match_s, mismatch_s, gap_s, in_seq_dict, original_len, subset_mode):
    # {name: [(query_start, query_end), value_in_overlap, count_query, coverage, continuity, seq_length, q_length]}
    if not overlap_value:
        for this_name in overlap_names:
            # try:
            this_query_1, this_query_2 = q_range_dicts[this_name]['q']
            # except:
            #     print(overlap_names)
            #     exit()
            if edge_connections:
                this_coverage = float(this_name.split('cov_')[1].split('--')[0].split('\'')[0].split(':')[0])
            else:
                this_coverage = 1
            overlap_value[this_name] = [[this_query_1, this_query_2],
                                        match_s if hit_site_dicts[go_to][this_name][1] else -mismatch_s,
                                        int(hit_site_dicts[go_to][this_name][0] != '-'),
                                        this_coverage,
                                        True,
                                        len(in_seq_dict[this_name.split('--')[0]]),
                                        abs(this_query_2 - this_query_1) + 1]
        for gap in hit_site_dicts[((go_to-2) % original_len+1, go_to)]:
            for name in gap:
                if name in overlap_value:
                    overlap_value[name][1] -= gap_s
                    overlap_value[name][2] += 1
        for gap in hit_site_dicts[(go_to, go_to % original_len + 1)]:
            for name in gap:
                if name in overlap_value:
                    overlap_value[name][1] -= gap_s
                    overlap_value[name][2] += 1
    while go_to + overlap_count <= original_len \
            and (set(hit_site_dicts[go_to + overlap_count].keys()) == overlap_names
                 or (subset_mode and overlap_names.issubset(set(hit_site_dicts[go_to + overlap_count].keys())))):
        next_go_to = go_to + overlap_count
        for name in overlap_names:
            overlap_value[name][1] += match_s if hit_site_dicts[next_go_to][name][1] else -mismatch_s
            overlap_value[name][2] += int(hit_site_dicts[next_go_to][name][0] != '-')
        for gap in hit_site_dicts[(next_go_to, next_go_to % original_len+1)]:
            for name in gap:
                if name in overlap_value:
                    overlap_value[name][1] -= gap_s
                    overlap_value[name][2] += 1
        overlap_count += 1
    return overlap_value, overlap_count


def only_best_name_in_hit(hit_site_dicts, go_to, best_name, overlap_count, total_len):

    def in_clear_gap(this_gap_site):
        j = 0
        while j < len(hit_site_dicts[this_gap_site]):
            if best_name in hit_site_dicts[this_gap_site][j]:
                hit_site_dicts[this_gap_site][j] = {best_name: hit_site_dicts[this_gap_site][j][best_name]}
                j += 1
            else:
                del hit_site_dicts[this_gap_site][j]
    in_clear_gap(((go_to - 2) % total_len + 1, go_to))
    for i in range(0, overlap_count):
        hit_site_dicts[go_to + i] = {best_name: hit_site_dicts[go_to + i][best_name]}
        in_clear_gap((go_to + i, (go_to + i) % total_len + 1))


# If both the best_name and dropping_name have base in the gap of the edge of the overlap, it would be kept.
# remove edge gaps to prevent this risk.
def drop_name_from_hit(hit_site_dicts, go_to, drop_name, overlap_count, total_len):

    def in_drop_gap(this_gap_site):
        j = 0
        while j < len(hit_site_dicts[this_gap_site]):
            if drop_name in hit_site_dicts[this_gap_site][j]:
                del hit_site_dicts[this_gap_site][j][drop_name]
                if len(hit_site_dicts[this_gap_site][j]):
                    j += 1
                else:
                    del hit_site_dicts[this_gap_site][j]
            else:
                j += 1
    in_drop_gap(((go_to - 2) % total_len + 1, go_to))
    for i in range(0, overlap_count):
        # KeyError: "EDGE_13_length_9487_cov_79.84273:EDGE_14_length_11143_cov_288.13731'--Hsp_num_1"
        # for 20161109-Sequences/FC708.fastg
        # skip this time@20161110, 22:25
        if drop_name in hit_site_dicts[go_to + i]:
            del hit_site_dicts[go_to + i][drop_name]
        in_drop_gap((go_to + i, (go_to + i) % total_len + 1))


def remove_multiple_queries_per_hit(hit_site_dicts, q_range_dicts, total_len, circular_h,
                                    sequential_name, edge_connections, q_cut_off, match_s, mismatch_s, gap_s,
                                    name_cluster, to_cluster_info, in_seq_dict, verbose, log):
    log.info("Removing redundant queries per hit ...")
    #
    go_to = 1
    original_len = len(hit_site_dicts)//2
    while go_to <= original_len:
        if len(hit_site_dicts[go_to]) == 0:
            go_to += 1
        elif len(hit_site_dicts[go_to]) == 1:
            this_name = list(hit_site_dicts[go_to])[0]
            if not sequential_name or this_name != sequential_name[-1]:
                sequential_name.append(this_name)
            go_to += 1
        else:
            overlap_value, overlap_count = initialize_overlap_value(set(hit_site_dicts[go_to].keys()), {}, 1, go_to,
                                                                    hit_site_dicts, q_range_dicts,
                                                                    edge_connections, match_s, mismatch_s,
                                                                    gap_s, in_seq_dict, original_len, False)
            # except:
            #     print(hit_site_dicts[go_to], go_to)
            #     exit()
            overlaps = sort_overlaps(overlap_value, overlap_count, to_cluster_info, name_cluster, q_range_dicts,
                                     hit_site_dicts, total_len, circular_h)
            best_name = overlaps[0][0]
            ###
            for other_name, other_value in overlaps[1:]:
                # print('other name =', other_name)
                # if best_name == "EDGE_39_length_201_cov_64.199:EDGE_38_length_6920_cov_57.04986'--Hsp_num_2":
                #     # print(overlaps)
                    # print('best value', overlaps[0][1][1] * overlaps[0][1][3]/ave_depth, math.log(overlaps[0][1][5]), overlaps[0][1][6]/overlap_count)
                    # print('drop value', overlaps[1][1][1] * overlaps[1][1][3]/ave_depth, math.log(overlaps[1][1][5]), overlaps[1][1][6]/overlap_count)
                if other_value[2]/float(other_value[6]) > q_cut_off:
                    # delete the whole hit
                    this_short_name, this_name_id = to_cluster_info[strip_new(other_name)]
                    delete_whole_hit_of_name(hit_site_dicts, q_range_dicts, name_cluster, this_short_name, this_name_id)
                    to_cluster_info = update_to_cluster(name_cluster)
                    # drop_name_from_hit(hit_site_dicts, go_to, o, overlap_count, total_len)
                else:
                    # start point
                    direction_trans = int((q_range_dicts[other_name]['dire'] - 0.5) * 2)
                    if go_to == q_range_dicts[other_name]['h'][0]:
                        q_range_dicts[other_name]['q'][0] += direction_trans * other_value[2]
                        q_range_dicts[other_name]['h'][0] += overlap_count
                        drop_name_from_hit(hit_site_dicts, go_to, other_name, overlap_count, total_len)
                    # end point
                    elif go_to+overlap_count-1 == q_range_dicts[other_name]['h'][1]:
                        q_range_dicts[other_name]['q'][1] -= direction_trans * other_value[2]
                        q_range_dicts[other_name]['h'][1] -= overlap_count
                        drop_name_from_hit(hit_site_dicts, go_to, other_name, overlap_count, total_len)
                    # middle
                    else:
                        pw_value, pw_count = initialize_overlap_value({best_name, other_name},
                                                                      {pair_name: overlap_value[pair_name][:]
                                                                       for pair_name in (best_name, other_name)},
                                                                      overlap_count, go_to, hit_site_dicts,
                                                                      q_range_dicts, edge_connections, match_s,
                                                                      mismatch_s, gap_s, in_seq_dict, original_len,
                                                                      True)
                        pw_overlaps = sort_overlaps(pw_value, pw_count, to_cluster_info, name_cluster, q_range_dicts,
                                                    hit_site_dicts, total_len, circular_h)
                        if verbose and best_name != pw_overlaps[0][0]:
                            log.info("Best name changed:"+32*' '+'old '+best_name+'\n'+32*' '+'new '+pw_overlaps[0][0])
                        best_name = pw_overlaps[0][0]
                        drop_name, drop_values = pw_overlaps[1]
                        direction_trans = int((q_range_dicts[drop_name]['dire'] - 0.5) * 2)
                        if drop_values[2]/float(drop_values[6]) > q_cut_off:
                            this_short_name, this_name_id = to_cluster_info[strip_new(drop_name)]
                            delete_whole_hit_of_name(hit_site_dicts, q_range_dicts, name_cluster, this_short_name,
                                                     this_name_id)
                            to_cluster_info = update_to_cluster(name_cluster)
                        else:
                            if go_to == q_range_dicts[drop_name]['h'][0]:
                                q_range_dicts[drop_name]['q'][0] += direction_trans * drop_values[2]
                                q_range_dicts[drop_name]['h'][0] += pw_count
                                drop_name_from_hit(hit_site_dicts, go_to, drop_name, pw_count, total_len)
                            elif go_to+pw_count-1 == q_range_dicts[other_name]['h'][1]:
                                q_range_dicts[drop_name]['q'][1] -= direction_trans * drop_values[2]
                                q_range_dicts[drop_name]['h'][1] -= pw_count
                                drop_name_from_hit(hit_site_dicts, go_to, drop_name, pw_count, total_len)
                            else:
                                # simple criteria
                                if verbose:
                                    log.info('Deleting anchor: '+drop_name)
                                this_short_name, this_name_id = to_cluster_info[strip_new(drop_name)]
                                delete_whole_hit_of_name(hit_site_dicts, q_range_dicts, name_cluster, this_short_name,
                                                         this_name_id)
                                to_cluster_info = update_to_cluster(name_cluster)
            #####
            if not sequential_name or best_name != sequential_name[-1]:
                sequential_name.append(best_name)
            go_to += overlap_count
    ########################
    for query_name, query_info in list(q_range_dicts.items()):
        if (query_info['q'][0] > query_info['q'][1]) == query_info['dire']:
            del q_range_dicts[query_name]
    count_sequential_name = 0
    while count_sequential_name < len(sequential_name):
        if sequential_name[count_sequential_name] not in q_range_dicts:
            del sequential_name[count_sequential_name]
        else:
            count_sequential_name += 1
    log.info("Removing redundant queries per hit finished")
    # return sequential_name, to_cluster_info


def del_overlap_sites(go_to_site, in_query_gap, in_name_short_1, h_site_dicts, total_len):
    in_count_del = 0
    while in_count_del < abs(in_query_gap):
        # print('del site', hit_site_dicts[count_site])
        if not h_site_dicts[go_to_site] \
                or list(h_site_dicts[go_to_site].values()).pop()[0] != '-':
            in_count_del += 1
        h_site_dicts[go_to_site] = {in_name_short_1: '-'}
        in_this_gap_site = ((go_to_site - 1) % total_len, go_to_site)
        while in_count_del < abs(in_query_gap) and h_site_dicts[in_this_gap_site]:
            # print('del gap', hit_site_dicts[this_gap_site][-1])
            del h_site_dicts[in_this_gap_site][-1]
            in_count_del += 1
        go_to_site = (go_to_site - 1) % total_len


def fill_gaps_with_query_seq(h_site_dicts, total_len, in_h_gap_end, in_h_gap_start, in_seq_to_add,
                             add_seq_name, name_hit_and_back):
    if in_h_gap_end == total_len and in_h_gap_start == 1:
        while in_seq_to_add:
            h_site_dicts[(total_len, 1)].append({add_seq_name: in_seq_to_add.pop(0)})
    elif in_h_gap_end - in_h_gap_start + 1 < 0:
        get_from_start = True
        for site_to_add in generate_from_ends(list(range(in_h_gap_start, total_len + in_h_gap_end + 1))):
            if in_seq_to_add:
                if get_from_start:
                    h_site_dicts[(site_to_add - 1) % total_len + 1] = {add_seq_name: in_seq_to_add.pop(0)}
                    get_from_start = False
                else:
                    h_site_dicts[(site_to_add - 1) % total_len + 1] = {add_seq_name: in_seq_to_add.pop()}
                    get_from_start = True
            else:
                h_site_dicts[(site_to_add - 1) % total_len + 1] = {add_seq_name: '-'}
        middle = in_h_gap_start + (total_len + in_h_gap_end - in_h_gap_start) // 2
        while in_seq_to_add:
            h_site_dicts[((middle - 1) % total_len + 1, middle % total_len + 1)].append(
                {add_seq_name: in_seq_to_add.pop(0)})
    else:
        get_from_start = True
        for site_to_add in generate_from_ends(list(range(in_h_gap_start, in_h_gap_end + 1))):
            if in_seq_to_add:
                if get_from_start:
                    h_site_dicts[site_to_add] = {add_seq_name: in_seq_to_add.pop(0)}
                    get_from_start = False
                else:
                    h_site_dicts[site_to_add] = {add_seq_name: in_seq_to_add.pop()}
                    get_from_start = True
            else:
                h_site_dicts[site_to_add] = {add_seq_name: '-'}
        middle = in_h_gap_start + (in_h_gap_end - in_h_gap_start) // 2
        middle_gap = (middle, middle % total_len + 1)
        # if hit is continuous, but there are extra bases to be add according to the query
        if h_site_dicts[middle_gap] and list(h_site_dicts[middle_gap][0].keys()).pop() == name_hit_and_back:
            # if there are extra bases of name2 that are not deleted
            insert_position = -len(h_site_dicts[middle_gap][0])
            while in_seq_to_add:
                h_site_dicts[middle_gap].insert(insert_position, {add_seq_name: in_seq_to_add.pop(0)})
        else:
            # if there are no extra bases or there are extra bases of name1 that are not deleted
            while in_seq_to_add:
                h_site_dicts[middle_gap].append({add_seq_name: in_seq_to_add.pop(0)})


def merge_hit_site_dicts(h_site_dicts, q_range_dicts, q_range_sets, sequential_name, in_seq_dict, total_len,
                         options, edge_connections, kmer, merged_names, log):
    log.info("Merging site dicts ...")
    min_overlap = options.min_overlap
    for query_name, query_info in q_range_dicts.items():
        short_name = query_name.split('--')[0]
        if short_name in q_range_sets:
            q_range_sets[short_name].append(sorted(query_info['q']))
        else:
            q_range_sets[short_name] = [sorted(query_info['q'])]
    for short_name in q_range_sets:
        q_range_sets[short_name].sort()

    i = 0
    len_seq_names = len(sequential_name)
    while i < len_seq_names:
        name1 = sequential_name[i]
        name2 = sequential_name[(i+1) % len_seq_names]
        name1_info = q_range_dicts[name1]
        name2_info = q_range_dicts[name2]
        direction_1 = name1_info['dire']
        direction_2 = name2_info['dire']
        direction_trans_1 = int(2 * (direction_1 - 0.5))
        direction_trans_2 = int(2 * (direction_2 - 0.5))
        name_short_1 = name1.split('--')[0]
        name_short_2 = name2.split('--')[0]
        h_gap_start = ((name1_info['h'][1] + 1) - 1) % total_len + 1
        h_gap_end = ((name2_info['h'][0] - 1) - 1) % total_len + 1
        hit_gap = (h_gap_end - h_gap_start + 1) % total_len
        this_seq_1 = in_seq_dict[name_short_1]
        this_seq_2 = in_seq_dict[name_short_2]
        merged = False

        """merging the gaps"""
        # just delete the overlap
        # no need to change q_range_dicts if it is no more used
        if name_short_1 == name_short_2:
            """case 1"""
            if direction_1 == direction_2:
                q_gap_start = name1_info['q'][1] + direction_trans_1
                q_gap_end = name2_info['q'][0] - direction_trans_1
                query_gap = direction_trans_1*(q_gap_end - q_gap_start) + 1
                max_gap = max(query_gap, hit_gap)
                if not bool(query_gap) or not bool(max_gap) or max_gap <= options.max_gap\
                        or abs(query_gap - hit_gap) / float(max_gap) <= options.max_dif:
                    # if there is a gap between two hits but smooth connection or overlap between two queries
                    # delete these overlap firstly
                    if query_gap <= 0 < direction_trans_1*(name2_info['q'][0]-name1_info['q'][0]):
                        del_overlap_sites(name1_info['h'][1], query_gap, name_short_1, h_site_dicts, total_len)
                        merged = True
                        merged_names[(name1, direction_1)] = None
                        merged_names[(name2, not direction_2)] = None
                    # fill the gaps in hits with base or '-'
                    if query_gap <= 0 < direction_trans_1*(name2_info['q'][0]-name1_info['q'][0]) or \
                            not contain_other_query(q_gap_start, q_gap_end, direction_1, q_range_sets[name_short_1]):
                        if query_gap <= 0:
                            seq_to_add = []
                        else:
                            if direction_1:
                                seq_to_add = list(this_seq_1[q_gap_start-1:q_gap_end])
                            else:
                                seq_to_add = list(complementary_seq(this_seq_1[q_gap_end-1:q_gap_start]))
                        fill_gaps_with_query_seq(h_site_dicts, total_len, h_gap_end, h_gap_start, seq_to_add,
                                                 name_short_1, name2)
                        merged = True
                        merged_names[(name1, direction_1)] = None
                        merged_names[(name2, not direction_2)] = None
        elif edge_connections:
            """case 2"""
            edge_name_1 = '_'.join(name_short_1.split('_')[1:]).split('_length')[0]
            edge_name_2 = '_'.join(name_short_2.split('_')[1:]).split('_length')[0]
            # Because the reverse contigs of fastg file are removed via del_complementary(),
            # (not name_short_1.split(';')[0].split(':')[0].endswith('\'')) should normally always be True.
            # thus I use name1_info['dire'] rather than
            # (not name_short_1.split(';')[0].split(':')[0].endswith('\'')) == name1_info['dire']
            if (edge_name_2, direction_2) in edge_connections[edge_name_1][direction_1]:
                len_seq_1 = len(this_seq_1)
                len_seq_2 = len(this_seq_2)
                end_1 = name1_info['q'][1]
                start_2 = name2_info['q'][0]
                if direction_1:
                    gap_1 = len_seq_1 - end_1
                else:
                    gap_1 = end_1 - 1
                if direction_2:
                    gap_2 = start_2 - 1
                else:
                    gap_2 = len_seq_2 - start_2
                query_gap = gap_1 + gap_2 - kmer
                max_gap = max(query_gap, hit_gap)
                if not max_gap or max_gap <= options.max_gap \
                        or abs(query_gap - hit_gap) / float(max_gap) <= options.max_dif:
                    if gap_1 == 0:
                        break_1 = False
                    elif direction_1:
                        break_1 = contain_other_query(end_1 + 1, len_seq_1, True, q_range_sets[name_short_1])
                    else:
                        break_1 = contain_other_query(gap_1, 1, False, q_range_sets[name_short_1])
                    if gap_2 == 0:
                        break_2 = False
                    elif direction_2:
                        break_2 = contain_other_query(1, gap_2, True, q_range_sets[name_short_2])
                    else:
                        break_2 = contain_other_query(len_seq_2, start_2 + 1, False, q_range_sets[name_short_2])
                    if (not break_1) and (not break_2):
                        if query_gap <= 0:
                            # if there is a gap between two hits but smooth connection or overlap between two queries
                            # delete these overlap firstly
                            # add the gap between two hits with '-'
                            del_overlap_sites(name1_info['h'][1], query_gap, name_short_1, h_site_dicts, total_len)
                            seq_to_add = []
                        else:
                            # add the gap between two hits with bases or '-'
                            this_seq_1 = in_seq_dict[name_short_1]
                            if direction_1:
                                seq_to_add_1 = this_seq_1[end_1:]
                            else:
                                seq_to_add_1 = complementary_seq(this_seq_1[:end_1 - 1])
                            this_seq_2 = in_seq_dict[name_short_2]
                            if direction_2:
                                seq_to_add_2 = this_seq_2[:start_2 - 1]
                            else:
                                seq_to_add_2 = complementary_seq(this_seq_2[start_2:])
                            if gap_1 >= kmer:
                                seq_to_add = list(seq_to_add_1[:-kmer] + seq_to_add_2)
                            else:
                                seq_to_add = list(seq_to_add_2[kmer - gap_1:])
                        name_to_add = 'EDGE_' + edge_name_1 + '\''*direction_1 + \
                                      '--EDGE_' + edge_name_2 + '\''*direction_2
                        fill_gaps_with_query_seq(h_site_dicts, total_len, h_gap_end, h_gap_start, seq_to_add,
                                                 name_to_add, name2)
                        merged = True
                        merged_names[(name1, direction_1)] = None
                        merged_names[(name2, not direction_2)] = None
        # detect overlap with minimum overlap information
        if not hit_gap and not merged:
            overlap_count = 0
            seq_1 = in_seq_dict[name_short_1]
            seq_2 = in_seq_dict[name_short_2]
            if not direction_1:
                seq_1 = complementary_seq(seq_1)[::-1]
            if not direction_2:
                seq_2 = complementary_seq(seq_2)[::-1]

            go_to_seq_1 = name1_info['q'][1] + direction_trans_1 - 1
            go_to_seq_2 = name2_info['q'][0] - 1
            while 0 <= go_to_seq_1 < len(seq_1) and 0 <= go_to_seq_2 < len(seq_2) \
                    and seq_1[go_to_seq_1] == seq_2[go_to_seq_2]:
                overlap_count += 1
                go_to_seq_1 += direction_trans_1
                go_to_seq_2 += direction_trans_2
                if overlap_count >= min_overlap:
                    merged = True
                    merged_names[(name1, direction_1)] = go_to_seq_1
                    merged_names[(name2, not direction_2)] = name2_info['q'][0] - 1
                    break
            go_to_seq_1 = name1_info['q'][1] - 1
            go_to_seq_2 = name2_info['q'][0]-direction_trans_2 - 1
            while 0 <= go_to_seq_1 < len(seq_1) and 0 <= go_to_seq_2 < len(seq_2) \
                    and seq_1[go_to_seq_1] == seq_2[go_to_seq_2]:
                overlap_count += 1
                go_to_seq_1 -= direction_trans_1
                go_to_seq_2 -= direction_trans_2
                if overlap_count >= min_overlap:
                    merged = True
                    merged_names[(name1, direction_1)] = name1_info['q'][1] - 1
                    merged_names[(name2, not direction_2)] = go_to_seq_2
                    break

        # if not merged, expand to half of the unused bases.
        if not merged:  # and options.expand_percent and options.expand_limit
            if name2_info['h'][0] < name1_info['h'][1]:
                add_h_middle = total_len
            else:
                add_h_middle = (name2_info['h'][0] + name1_info['h'][1]) // 2
            ########
            add_h_from_1 = (name1_info['h'][1] % total_len) + 1
            add_h_to_1 = add_h_middle
            if direction_1:
                for this_q_range in q_range_sets[name_short_1] + [(len(this_seq_1) + 1, len(this_seq_1) + 1)]:
                    if name1_info['q'][1] < this_q_range[0]:
                        add_q_from_1 = name1_info['q'][1]
                        if this_q_range[0] == len(this_seq_1) + 1:
                            add_q_to_1 = add_q_from_1 + min(len(this_seq_1) - add_q_from_1, options.expand_limit)
                        else:
                            add_q_to_1 = int(add_q_from_1 + this_q_range[0]) / 2.
                            add_q_to_1 = add_q_from_1 + min(options.expand_limit,
                                                            int(options.expand_percent * (add_q_to_1 - add_q_from_1)))
                        seq_to_add_here_1 = this_seq_1[add_q_from_1:add_q_to_1]
                        break
            else:
                for this_q_range in q_range_sets[name_short_1][::-1] + [(0, 0)]:
                    if name1_info['q'][1] > this_q_range[1]:
                        add_q_to_1 = name1_info['q'][1] - 1
                        if this_q_range[1] == 0:
                            add_q_from_1 = add_q_to_1 - min(options.expand_limit, add_q_to_1)
                        else:
                            add_q_from_1 = int((name1_info['q'][1] + this_q_range[1]) / 2.)
                            add_q_from_1 = add_q_to_1 - min(options.expand_limit,
                                                            int(options.expand_percent * (add_q_to_1 - add_q_from_1)))
                        seq_to_add_here_1 = complementary_seq(this_seq_1[add_q_from_1:add_q_to_1])
                        break
            ########
            add_h_from_2 = add_h_middle % total_len + 1
            add_h_to_2 = ((name2_info['h'][0] - 2) % total_len) + 1
            if direction_2:
                for this_q_range in q_range_sets[name_short_2][::-1] + [(0, 0)]:
                    if name2_info['q'][0] > this_q_range[1]:
                        add_q_to_2 = name2_info['q'][0] - 1
                        if this_q_range[1] == 0:
                            add_q_from_2 = add_q_to_2 - min(options.expand_limit, add_q_to_2)
                        else:
                            add_q_from_2 = int((name2_info['q'][0] + this_q_range[1]) / 2.)
                            add_q_from_2 = add_q_to_2 - min(options.expand_limit,
                                                            int(options.expand_percent * (add_q_to_2 - add_q_from_2)))
                        seq_to_add_here_2 = this_seq_2[add_q_from_2:add_q_to_2]
                        break
            else:
                for this_q_range in q_range_sets[name_short_2] + [(len(this_seq_2) + 1, len(this_seq_2) + 1)]:
                    if name2_info['q'][1] < this_q_range[0]:
                        add_q_from_2 = name2_info['q'][0]
                        if this_q_range[0] == len(this_seq_2) + 1:
                            add_q_to_2 = add_q_from_2 + min(options.expand_limit, len(this_seq_2) - add_q_from_2)
                        else:
                            add_q_to_2 = (name2_info['q'][0] + this_q_range[0]) / 2.
                            add_q_to_2 = add_q_from_2 + min(options.expand_limit,
                                                            int(options.expand_percent * (add_q_to_2 - add_q_from_2)))
                        seq_to_add_here_2 = complementary_seq(this_seq_2[add_q_from_2:add_q_to_2])
                        break
            ########
            add_len_1 = len(seq_to_add_here_1)
            add_len_2 = len(seq_to_add_here_2)
            if direction_1:
                seq_to_check_1 = this_seq_1[max(add_q_from_1 - min_overlap, 0):add_q_to_1]
            else:
                seq_to_check_1 = complementary_seq(this_seq_1[add_q_from_1:add_q_to_1 + min_overlap])
            if direction_2:
                seq_to_check_2 = this_seq_2[add_q_from_2:add_q_to_2 + min_overlap]
            else:
                seq_to_check_2 = complementary_seq(this_seq_2[max(add_q_from_2 - min_overlap, 0):add_q_to_2])
            check_len_1 = len(seq_to_check_1)
            check_len_2 = len(seq_to_check_2)
            if check_len_1 >= min_overlap and check_len_2 >= min_overlap:
                # risky to use dict here because if there were 2-time repeats both in the two contigs, the result would
                # include a 3-time repeat.
                candidate_overlap_words_1 = {seq_to_check_1[l:l+min_overlap]: l
                                             for l in range(check_len_1-min_overlap, -1, -1)}
                candidate_overlap_words_2 = [seq_to_check_2[l:l+min_overlap] for l in range(check_len_2-min_overlap+1)]
                for word_index_2 in range(len(candidate_overlap_words_2)):
                    if candidate_overlap_words_2[word_index_2] in candidate_overlap_words_1:
                        word_index_1 = candidate_overlap_words_1[candidate_overlap_words_2[word_index_2]]
                        merged = True
                        break
            if merged:
                point_to_seq_1 = word_index_1 - (check_len_1 - add_len_1)
                point_to_seq_2 = word_index_2
                extra = 10
                skip = 3
                if options.verbose:
                    max_print = max(word_index_1, word_index_2)
                    start_differ = abs(word_index_1 - word_index_2)
                    start_cut_print = start_differ - min(start_differ, extra)
                    end_differ = abs((check_len_1 - word_index_1) - (check_len_2 - word_index_2))
                    end_cut_print = end_differ - min(end_differ, extra)
                    if start_cut_print:
                        if word_index_1 > word_index_2:
                            prefix_print_1 = '.' * skip + seq_to_check_1[word_index_1 - extra:word_index_1]
                            prefix_print_2 = ' ' * (skip + extra) + seq_to_check_2[:word_index_2]
                        else:
                            prefix_print_1 = ' ' * skip + seq_to_check_1[:word_index_1]
                            prefix_print_2 = '.' * (skip + extra) + seq_to_check_2[word_index_2 - extra:word_index_2]
                    else:
                        prefix_print_1 = ' ' * (max_print - word_index_1) + seq_to_check_1[:word_index_1]
                        prefix_print_2 = ' ' * (max_print - word_index_2) + seq_to_check_2[:word_index_2]
                    if end_cut_print:
                        if check_len_1 - word_index_1 > check_len_2 - word_index_2:
                            postfix_print_1 = seq_to_check_1[word_index_1 + min_overlap: -end_cut_print] + '.' * skip
                            postfix_print_2 = seq_to_check_2[word_index_2 + min_overlap:]
                        else:
                            postfix_print_1 = seq_to_check_1[word_index_1 + min_overlap:]
                            postfix_print_2 = seq_to_check_2[word_index_2 + min_overlap: - end_cut_print] + '.' * skip
                    else:
                        postfix_print_1 = seq_to_check_1[word_index_1 + min_overlap:]
                        postfix_print_2 = seq_to_check_2[word_index_2 + min_overlap:]
                    seq_to_print_1 = prefix_print_1 + '*' + seq_to_check_1[word_index_1:word_index_1 + min_overlap] +\
                                     '*' + postfix_print_1
                    seq_to_print_2 = prefix_print_2 + '*' + seq_to_check_2[word_index_2:word_index_2 + min_overlap] +\
                                     '*' + postfix_print_2
                    equal_marks = '' + ' '*3*int(bool(start_cut_print))
                    for go_to_base in range(3*int(bool(start_cut_print)), min(len(seq_to_print_1), len(seq_to_print_2))):
                        if seq_to_print_1[go_to_base] == seq_to_print_2[go_to_base] == '*':
                            equal_marks += ' '
                        elif seq_to_print_1[go_to_base] == seq_to_print_2[go_to_base] != ' ':
                            equal_marks += '|'
                        else:
                            equal_marks += ' '
                    log.info("Expanded and merged: " + name1 + "(" +
                             str(add_len_1) + ") " + name2 + "(" + str(len(seq_to_add_here_2)) + ")\n" + ' ' * 32 +
                             seq_to_print_1 + '\n' + ' ' * 32 + equal_marks + '\n' + ' ' * 32 + seq_to_print_2)
                if point_to_seq_1 < 0:
                    del_overlap_sites(name1_info['h'][1], -point_to_seq_1, name_short_1, h_site_dicts, total_len)
                fill_gaps_with_query_seq(h_site_dicts, total_len, add_h_to_2, add_h_from_1,
                                         list(seq_to_add_here_1[:max(point_to_seq_1, 0)]
                                              + seq_to_add_here_2[min(point_to_seq_2, len(seq_to_add_here_2)):]),
                                         name_short_1 + '+' + name_short_2, name2)
                merged_names[(name1, direction_1)] = (add_q_from_1, add_q_to_1)[direction_1] + point_to_seq_1
                merged_names[(name2, not direction_2)] = (add_q_from_2, add_q_to_2)[not direction_2] - word_index_2
        i += 1
    log.info("Merging site dicts finished.")


def extend_unmerged(h_site_dicts, q_range_dicts, q_range_sets, sequential_name, in_seq_dict, name_cluster,
                    to_cluster_info, total_len, options, merged_names, add_gap, log):
    log.info("Extending unmerged site dicts ...")
    i = 0
    len_seq_names = len(sequential_name)
    # for jack in sequential_name:
    #     print(jack)
    while i < len_seq_names:
        name1 = sequential_name[i]
        name2 = sequential_name[(i + 1) % len_seq_names]
        name1_info = q_range_dicts[name1]
        name2_info = q_range_dicts[name2]
        direction_1 = name1_info['dire']
        direction_2 = name2_info['dire']
        direction_trans_1 = int(2 * (direction_1 - 0.5))
        direction_trans_2 = int(2 * (direction_2 - 0.5))
        if (name1, direction_1) in merged_names:
            # if name1 in {"FC202-short--Hsp_num_31", "FC202-short--Hsp_num_32"}:
            #     print("merged !!", name1, name2)
            pass
        else:
            name_short_1 = name1.split('--')[0]
            name_short_2 = name2.split('--')[0]
            this_seq_1 = in_seq_dict[name_short_1]
            this_seq_2 = in_seq_dict[name_short_2]

            if name2_info['h'][0] < name1_info['h'][1]:
                add_h_middle = total_len
            else:
                add_h_middle = (name2_info['h'][0] + name1_info['h'][1]) // 2
            ########
            add_h_from_1 = (name1_info['h'][1] % total_len) + 1
            add_h_to_1 = add_h_middle
            #
            short_name_1, name_1_id = to_cluster_info[strip_new(name1)]
            next_to_name_1_id = name_1_id + direction_trans_1
            modified_q_position_1 = -1
            if 0 <= next_to_name_1_id < len(name_cluster[short_name_1]):
                next_to_name_1 = name_cluster[short_name_1][next_to_name_1_id]
                if (next_to_name_1, not direction_1) in merged_names:
                    try:
                        modified_q_position_1 = merged_names[(next_to_name_1, not direction_1)] - direction_trans_1
                    except TypeError:
                        pass
            if direction_1:
                for this_q_range in q_range_sets[name_short_1] + [(len(this_seq_1) + 1, len(this_seq_1) + 1)]:
                    if name1_info['q'][1] < this_q_range[0]:
                        add_q_from_1 = name1_info['q'][1]
                        if modified_q_position_1 < 0:
                            if this_q_range[0] == len(this_seq_1) + 1:
                                add_q_to_1 = len(this_seq_1)
                            else:
                                add_q_to_1 = int(add_q_from_1 + this_q_range[0]) / 2.
                        else:
                            add_q_to_1 = modified_q_position_1 + 1
                        add_q_to_1 = add_q_from_1 + min(options.expand_limit, int(
                            options.expand_percent * (add_q_to_1 - add_q_from_1)))
                        seq_to_add_here_1 = this_seq_1[add_q_from_1:add_q_to_1]
                        break
            else:
                for this_q_range in q_range_sets[name_short_1][::-1] + [(0, 0)]:
                    if name1_info['q'][1] > this_q_range[1]:
                        add_q_to_1 = name1_info['q'][1] - 1
                        if modified_q_position_1 < 0:
                            if this_q_range[1] == 0:
                                add_q_from_1 = 0
                            else:
                                add_q_from_1 = int((name1_info['q'][1] + this_q_range[1]) / 2.)
                        else:
                            add_q_from_1 = modified_q_position_1
                        add_q_from_1 = add_q_to_1 - min(options.expand_limit, int(
                            options.expand_percent * (add_q_to_1 - add_q_from_1)))
                        seq_to_add_here_1 = complementary_seq(this_seq_1[add_q_from_1:add_q_to_1])
                        break
            ########
            add_h_from_2 = add_h_middle % total_len + 1
            add_h_to_2 = ((name2_info['h'][0] - 2) % total_len) + 1
            short_name_2, name_2_id = to_cluster_info[strip_new(name2)]
            prior_to_name_2_id = name_2_id - direction_trans_2
            modified_q_position_2 = -1
            if 0 <= prior_to_name_2_id < len(name_cluster[short_name_2]):
                prior_to_name_2 = name_cluster[short_name_2][prior_to_name_2_id]
                if (prior_to_name_2, direction_2) in merged_names:
                    try:
                        modified_q_position_2 = merged_names[(prior_to_name_2, direction_2)] + direction_trans_2
                    except TypeError:
                        pass
            if direction_2:
                for this_q_range in q_range_sets[name_short_2][::-1] + [(0, 0)]:
                    if name2_info['q'][0] > this_q_range[1]:
                        add_q_to_2 = name2_info['q'][0] - 1
                        if modified_q_position_2 < 0:
                            if this_q_range[1] == 0:
                                add_q_from_2 = 0
                            else:
                                add_q_from_2 = int((name2_info['q'][0] + this_q_range[1]) / 2.)
                        else:
                            add_q_from_2 = modified_q_position_2
                        add_q_from_2 = add_q_to_2 - min(options.expand_limit, int(
                            options.expand_percent * (add_q_to_2 - add_q_from_2)))
                        seq_to_add_here_2 = this_seq_2[add_q_from_2:add_q_to_2]
                        break
            else:
                for this_q_range in q_range_sets[name_short_2] + [(len(this_seq_2) + 1, len(this_seq_2) + 1)]:
                    if name2_info['q'][1] < this_q_range[0]:
                        add_q_from_2 = name2_info['q'][0]
                        if modified_q_position_2 < 0:
                            if this_q_range[0] == len(this_seq_2) + 1:
                                add_q_to_2 = len(this_seq_2)
                            else:
                                add_q_to_2 = (name2_info['q'][0] + this_q_range[0]) / 2.
                        else:
                            add_q_to_2 = modified_q_position_2 + 1
                        add_q_to_2 = add_q_from_2 + min(options.expand_limit, int(
                            options.expand_percent * (add_q_to_2 - add_q_from_2)))
                        seq_to_add_here_2 = complementary_seq(this_seq_2[add_q_from_2:add_q_to_2])
                        break
            ########
            add_len_1 = len(seq_to_add_here_1)
            add_len_2 = len(seq_to_add_here_2)

            n_gap = "N" * max((add_h_to_2 - add_h_from_1 + 1) % total_len - add_len_1 - add_len_2, add_gap)

            n_gap_1 = n_gap[:len(n_gap) // 2]
            n_gap_2 = n_gap[len(n_gap) // 2:]
            fill_gaps_with_query_seq(h_site_dicts, total_len, add_h_to_1, add_h_from_1,
                                     list(seq_to_add_here_1 + n_gap_1), name_short_1, name2)
            fill_gaps_with_query_seq(h_site_dicts, total_len, add_h_to_2, add_h_from_2,
                                     list(n_gap_2 + seq_to_add_here_2), name_short_2, name2)
        i += 1
    log.info("Extending unmerged site dicts finished.")


def contain_other_query(t_s, t_e, direction, sorted_discrete_ranges):
    if not direction:
        t_s, t_e = t_e, t_s
    if t_e < t_s:
        if t_s < sorted_discrete_ranges[-1][1] or t_e > sorted_discrete_ranges[0][0]:
            return True
    else:
        for r_1, r_2 in sorted_discrete_ranges:
            if (t_s < r_1 < t_e) or (t_s < r_2 < t_e):
                return True
    return False


def generate_from_ends(seq_to_chop):
    while seq_to_chop:
        yield seq_to_chop.pop(0)
        if seq_to_chop:
            yield seq_to_chop.pop()


def parse_fastg(query_matrix, log):
    # ----------------------------------------
    # find start and end points of query
    # initialize candidates: fastg topologies and sequences
    log.info("Parsing fastg file ...")
    len_fastg = len(query_matrix[0])
    edge_infos = {}
    short_names = []
    for i in range(len_fastg):
        full_name = query_matrix[0][i]
        short_name = '_'.join(full_name.split()[0].split('_')[1:]).split('_length')[0]
        coverage = float(full_name.split('cov_')[1].split(';')[0].split('\'')[0].split(':')[0])
        edge_infos[short_name] = {False: set(), True: set(), 'coverage': coverage}
        short_names.append(short_name)
    for i in range(len_fastg):
        full_name = query_matrix[0][i]
        short_name = short_names[i]
        connected_edges = set()
        if ':' in full_name:
            for edge in full_name.rstrip(';').split(':')[1].split(','):
                edge_short_name = '_'.join(edge.split('_')[1:]).split('_length')[0]
                if edge_short_name in edge_infos:
                    if edge.endswith('\''):
                        connected_edges.add((edge_short_name, False))
                    else:
                        connected_edges.add((edge_short_name, True))
        if full_name.split(';')[0].split(':')[0].endswith('\''):
            sequence = query_matrix[1][i]
            new_items = {('index', False): i,
                         ('seq', False): sequence,
                         ('seq', True): complementary_seq(sequence),
                         False: connected_edges}
            edge_infos[short_name].update(new_items)
        else:
            sequence = query_matrix[1][i]
            len_seq = len(sequence)
            new_items = {'identity': [0]*len_seq,
                         'start_block': {'q': (len_seq, len_seq), 'r': []},
                         'end_block': {'q': (0, 0), 'r': []},
                         ('index', True): i,
                         ('seq', True): sequence,
                         ('seq', False): complementary_seq(sequence),
                         'len_seq': len_seq,
                         True: connected_edges}
            edge_infos[short_name].update(new_items)
    # -----------------------------------
    # detect k-mer
    k_mer = 0
    try:
        for short_name in edge_infos:
            for direction in [True, False]:
                for next_edge_info in edge_infos[short_name][direction]:
                    if k_mer:
                        if edge_infos[short_name][('seq', direction)][-k_mer:] != \
                                edge_infos[next_edge_info[0]][('seq', next_edge_info[1])][:k_mer]:
                            raise ValueError
                    else:
                        for k_mer in range(127, 19, -2):
                            if edge_infos[short_name][('seq', direction)][-k_mer:] == \
                                    edge_infos[next_edge_info[0]][('seq', next_edge_info[1])][:k_mer]:
                                break
                        else:
                            raise ValueError
    except ValueError:
        k_mer = 0
        pass
    # ------------------------------------
    if not k_mer:
        edge_infos = {}
        log.warning("Parsing fastg file failed.")
    else:
        log.info("Detect k-mer = "+str(k_mer))
        log.info("Parsing fastg file finished.")
    return edge_infos, k_mer


def combine_site_dict(add_item, total_site_dict, log):
    q_name, site_dicts = add_item
    log.info("Adding "+str(q_name)+" to final result ...")
    site_set = set(range(1, len(site_dicts)+1))
    for site in site_dicts:
        if site_dicts[site]:
            if site in site_set:
                total_site_dict[site][q_name] = list(site_dicts[site].values())[0]
            else:
                gap_sites = site_dicts[site]
                exist_len = len(total_site_dict[site])
                for i in range(len(gap_sites)):
                    if i < exist_len:
                        total_site_dict[site][i][q_name] = list(gap_sites[i].values())[0]
                        # if list(len(gap_sites[i].values())) != 1: print("Warning", gap_sites[i])
                    else:
                        total_site_dict[site].append({q_name: list(gap_sites[i].values())[0]})
        else:
            pass
    log.info("Adding "+str(q_name)+" to final result finished.")


def hit_site_dicts_to_sequence(hit_matrix, hit_site_dicts, filling, log):
    log.info("Hit sites to sequence ...")
    site_values = list(hit_site_dicts[1].values())
    if site_values:
        sequence = site_values[0][0]
    else:
        sequence = filling
    for go_to_site in range(2, len(hit_matrix[1][0]) + 1):
        for gap_dict in hit_site_dicts[(go_to_site - 1, go_to_site)]:
            site_values = list(gap_dict.values())
            if site_values:
                sequence += site_values[0][0]
        site_values = list(hit_site_dicts[go_to_site].values())
        if site_values:
            sequence += site_values[0][0]
        else:
            sequence += filling
    for gap_dict in hit_site_dicts[(go_to_site, 1)]:
        site_values = list(gap_dict.values())
        if site_values:
            sequence += site_values[0][0]
    log.info("Hit sites to sequence finished. ")
    return sequence


def hit_site_dicts_to_sequence_mark_conservative(hit_matrix, hit_site_dicts, filling, cont_cons_sites, log):
    log.info("Hit sites to sequence ...")
    site_values = list(hit_site_dicts[1].values())
    if site_values:
        sequence = site_values[0][0]
    else:
        sequence = filling
    seq_count = 1
    if 1 in cont_cons_sites:
        cont_cons_sites[1].append(seq_count)
    for go_to_site in range(2, len(hit_matrix[1][0]) + 1):
        for gap_dict in hit_site_dicts[(go_to_site - 1, go_to_site)]:
            site_values = list(gap_dict.values())
            if site_values:
                sequence += site_values[0][0]
                seq_count += 1
        site_values = list(hit_site_dicts[go_to_site].values())
        if site_values:
            sequence += site_values[0][0]
        else:
            sequence += filling
        seq_count += 1
        if go_to_site in cont_cons_sites:
            cont_cons_sites[go_to_site].append(seq_count)
    for gap_dict in hit_site_dicts[(go_to_site, 1)]:
        site_values = list(gap_dict.values())
        if site_values:
            sequence += site_values[0][0]
    log.info("Hit sites to sequence finished. ")
    return sequence


def get_range(continuous_conservative_sites, groups, seq_length, seq_number):
    split_points = []
    for group in groups:
        split_points.append((group[0] + group[-1]) // 2)
    this_range = [[0]]
    for split_p in split_points:
        corresponding = continuous_conservative_sites[split_p][seq_number]
        this_range[-1].append(corresponding)
        this_range.append([corresponding])
    this_range[-1].append(seq_length)
    return this_range


def get_groups(continuous_conservative_sites, insertions_in_conservative):
    groups = [[]]
    for candidate_site in sorted(list(continuous_conservative_sites)):
        groups[-1].append(candidate_site)
        if candidate_site + 1 not in continuous_conservative_sites \
                or (candidate_site, candidate_site + 1) in insertions_in_conservative:
            groups.append([])
    if not groups[-1]:
        del groups[-1]
    return groups


def check_conservative_continuous(continuous_conservative_sites, h_site_dicts, insertions, conservative_len, log):
    log.info("Checking conservative continuous sites ...")
    for candidate_site in sorted(list(continuous_conservative_sites)):
        if not h_site_dicts[candidate_site] or list(h_site_dicts[candidate_site].values())[0][1] == -1:
            del continuous_conservative_sites[candidate_site]
    #
    groups = get_groups(continuous_conservative_sites, insertions)
    #
    count_group = 0
    while count_group < len(groups):
        if len(groups[count_group]) < conservative_len:
            for candidate_site in groups[count_group]:
                del continuous_conservative_sites[candidate_site]
            del groups[count_group]
        else:
            count_group += 1
    #
    for count_group in range(len(groups)):
        count_member = 0
        while count_member < len(groups[count_group]) - 1:
            insertion_site = (groups[count_group][count_member], groups[count_group][count_member] + 1)
            if h_site_dicts[insertion_site]:
                groups.append(groups[count_group][:count_member+1])
                groups[count_group] = groups[count_group][count_member+1:]
                insertions.add(insertion_site)
                count_member = 0
            else:
                count_member += 1
    #
    count_group = 0
    while count_group < len(groups):
        if len(groups[count_group]) < conservative_len:
            for candidate_site in groups[count_group]:
                del continuous_conservative_sites[candidate_site]
            insertions.discard((groups[count_group][0]-1, groups[count_group][0]))
            insertions.discard((groups[count_group][-1], groups[count_group][-1]+1))
            del groups[count_group]
        else:
            count_group += 1
    # groups.sort(key=lambda x: min(x))
    log.info("Conservative continuous groups: "+str(len(groups)))
    log.info("Checking conservative continuous sites finished!")


def add_info_to_cc_sites(continuous_conservative_sites, sequence_with_gaps):
    count_seq = 0
    for base_count in range(len(sequence_with_gaps)):
        if sequence_with_gaps[base_count] != '-':
            count_seq += 1
            if count_seq in continuous_conservative_sites:
                continuous_conservative_sites[count_seq].append(base_count)


def alignment_multiple_with_hit_site_dicts(hit_matrix, hit_site_dicts, filling, log):
    log.info("Making alignments with site dicts ...")
    new_matrix = [[], [], hit_matrix[2]]
    new_matrix[0].append(hit_matrix[0][0])
    new_matrix[1].append([hit_matrix[1][0][0]])
    sys.stdout.write("Framing (1/2) - ")
    for go_to_site in range(2, len(hit_matrix[1][0])+1):
        new_matrix[1][0].append('-'*len(hit_site_dicts[(go_to_site - 1, go_to_site)]))
        new_matrix[1][0].append(hit_matrix[1][0][go_to_site-1])
        sys.stdout.write(str(go_to_site) + '\b' * len(str(go_to_site)))
        sys.stdout.flush()
    new_matrix[1][0] = ''.join(new_matrix[1][0])
    #
    sys.stdout.write('\b' * len("Framing (1/2) - "))
    sys.stdout.flush()
    sys.stdout.write("Framing (2/2) - ")
    names = set()
    for name in hit_site_dicts[1]:
        names.add(name)
    for go_to_site in range(2, len(hit_matrix[1][0])+1):
        for gap_dict in hit_site_dicts[(go_to_site-1, go_to_site)]:
            for name in gap_dict:
                names.add(name)
        for name in hit_site_dicts[go_to_site]:
            names.add(name)
        sys.stdout.write(str(go_to_site) + '\b' * len(str(go_to_site)))
        sys.stdout.flush()
    for gap_dict in hit_site_dicts[(go_to_site, 1)]:
        for name in gap_dict:
            names.add(name)
    #
    sys.stdout.write('\b' * len("(2/2) - "))
    sys.stdout.flush()
    sys.stdout.write("- Filling - ")
    query_seqs = {name: [hit_site_dicts[1].get(name, filling)[0]] for name in names}
    for go_to_site in range(2, len(hit_matrix[1][0])+1):
        for gap_dict in hit_site_dicts[(go_to_site-1, go_to_site)]:
            for name in names:
                query_seqs[name].append(gap_dict.get(name, '-')[0])
        for name in names:
            query_seqs[name].append(hit_site_dicts[go_to_site].get(name, filling)[0])
        sys.stdout.write(str(go_to_site)+'\b'*len(str(go_to_site)))
        sys.stdout.flush()
    #
    last_gaps = hit_site_dicts[(go_to_site, 1)]
    if last_gaps:
        new_matrix[1][0] += '-'*len(last_gaps)
        for gap_dict in last_gaps:
            for name in names:
                query_seqs[name].append(gap_dict.get(name, '-'))
    #
    sys.stdout.write("Joining")
    for name in query_seqs:
        new_matrix[0].append(name)
        new_matrix[1].append(''.join(query_seqs[name]))
    sys.stdout.write("\b" * len("Framing - Filling - Joining"))
    sys.stdout.flush()
    log.info("Making alignments with site dicts finished.")
    return new_matrix


def del_complementary(fastg_file):
    temp_matrix = read_fasta_gb_head(fasta_dir=fastg_file)
    i = 0
    while i < len(temp_matrix[0]):
        if temp_matrix[0][i].rstrip(';').split(':')[0].endswith('\''):
            del temp_matrix[0][i]
            del temp_matrix[1][i]
        else:
            i += 1
    write_fasta(out_dir=fastg_file + '.fasta', matrix=temp_matrix, overwrite=True)
    return fastg_file + '.fasta'


def read_fasta_gb_head(fasta_dir):
    fasta_file = open(fasta_dir, 'rU')
    names = []
    seqs = []
    this_line = fasta_file.readline()
    interleaved = 0
    while this_line:
        if this_line.startswith('>'):
            names.append(this_line[1:].rstrip('\n').rstrip('\r').rstrip(';').rstrip('-').rstrip('.').strip())
            this_seq = ''
            this_line = fasta_file.readline()
            seq_line_count = 0
            while this_line and not this_line.startswith('>'):
                if seq_line_count == 1:
                    interleaved = len(this_seq)
                this_seq += this_line.strip()
                this_line = fasta_file.readline()
                seq_line_count += 1
            seqs.append(this_seq.replace('?', "N"))
        else:
            this_line = fasta_file.readline()
    fasta_file.close()
    return [names, seqs, interleaved]


def write_fasta(out_dir, matrix, overwrite, remove_gap_and_n=False):
    global keep_n_len
    if not overwrite:
        while os.path.exists(out_dir):
            out_dir = '.'.join(out_dir.split('.')[:-1])+'_.'+out_dir.split('.')[-1]
    fasta_file = open(out_dir, 'w')
    if remove_gap_and_n:
        new_matrix = [matrix[0][:], [], matrix[2]]
        for sequence in matrix[1]:
            sequence = sequence.replace('-', '')
            new_sequence = []
            count_n = 0
            for base in sequence:
                if base == "N":
                    if count_n < keep_n_len:
                        new_sequence.append(base)
                    count_n += 1
                else:
                    new_sequence.append(base)
                    count_n = 0
            new_sequence = ''.join(new_sequence)
            new_matrix[1].append(new_sequence)
    else:
        new_matrix = matrix
    if new_matrix[2]:
        for i in range(len(new_matrix[0])):
            fasta_file.write('>'+new_matrix[0][i]+'\n')
            j = new_matrix[2]
            while j < len(new_matrix[1][i]):
                fasta_file.write(new_matrix[1][i][(j-new_matrix[2]):j]+'\n')
                j += new_matrix[2]
            fasta_file.write(new_matrix[1][i][(j-new_matrix[2]):j]+'\n')
    else:
        for i in range(len(new_matrix[0])):
            fasta_file.write('>'+new_matrix[0][i] + '\n')
            fasta_file.write(new_matrix[1][i] + '\n')
    fasta_file.close()


def constant_len_number(this_number, number_length):
    str_number = str(this_number)
    return (number_length - len(str_number))*"0"+str_number


def main():
    time0 = time.time()
    print_title = "\nAlign contigs to reference" \
                  "\n" \
                  "\nThis scripts would produce both mapping-like alignments" \
                  " and consensus alignments for phylogenetic analysis." \
                  "\nCopyright Aug 24 2016 Jianjun Jin (jinjianjun@mail.kib.ac.cn)" \
                  "\n"
    options, args, log = require_options(print_title)
    try:
        """Make blastdb"""
        del_randomly = False
        out_dir = options.output_directory
        database = check_db(options.reference_fasta, options.min_repeat, options.word_size, options.keep_repeat_ends,
                            del_randomly, options.circular_refer, options.add_gap_repeat, log)
        reference_matrix = read_fasta_gb_head(database[:-6])
        total_len = len(reference_matrix[1][0])
        if options.separate_result:
            log.info("Initializing continuous conservative sites ...\n")
            continuous_conservative_sites = {candidate_site: [] for candidate_site in range(1, total_len+1)}
            insertions_in_conservative = set()
        if options.no_final_aligned:
            final_out_alignment = [[], [], reference_matrix[2]]
        else:
            total_site_dict = initialize_site_dict(total_len)
        for count_fasta in range(len(args)):
            log.info(str(count_fasta+1)+'/'+str(len(args)))
            log.info("Analysing " + args[count_fasta] + ' ...')

            """Initialization and execute blast"""
            raw_seq_file = str(args[count_fasta])
            is_fastg = raw_seq_file.endswith('.fastg')
            if is_fastg:
                b_word_size = options.blast_word_size_fg
                b_evalue = options.blast_evalue_fg
                log.info("Pre-treating fastg ...")
                blast_input = del_complementary(raw_seq_file)
                if options.rm_low_coverage:
                    self_blast = subprocess.Popen("rm_low_coverage_duplicated_contigs.py --blur "+blast_input,
                                                  stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
                    output, error = self_blast.communicate()
                    if 'error' in str(output) or "Error" in str(output) or "command not found" in str(output):
                        log.warning(str(output).strip())
                    elif os.path.isfile(blast_input + '.purified.fastg'):
                        if not options.keep_temp:
                            os.remove(blast_input)
                        os.rename(blast_input + '.purified.fastg', blast_input + '.purified')
                        blast_input += '.purified'
                log.info("Pre-treating fastg finished.")
            else:
                b_word_size = options.blast_word_size
                b_evalue = options.blast_evalue
                if options.min_repeat:
                    modified = False
                    seqs_to_remove_repeats = read_fasta_gb_head(raw_seq_file)
                    for seq_id in range(len(seqs_to_remove_repeats[0])):
                        this_seq = seqs_to_remove_repeats[1][seq_id]
                        this_repeats = detect_repeats(this_seq, options.min_repeat, options.circular_query, log,
                                                      word_size=options.word_size)
                        if this_repeats and this_repeats[0]:
                            modified = True
                            seqs_to_remove_repeats[1][seq_id] = remove_repeats(this_seq, this_repeats,
                                                                               options.keep_repeat_ends,
                                                                               del_randomly,
                                                                               options.circular_query,
                                                                               options.add_gap_repeat, log)
                    if modified:
                        raw_seq_file += '.repeat_removed'
                        write_fasta(out_dir=raw_seq_file, matrix=seqs_to_remove_repeats, overwrite=True)
                    del seqs_to_remove_repeats
                blast_input = raw_seq_file

            this_base_name = re.sub(".fasta$", '', re.sub(".fastg$", '', os.path.basename(raw_seq_file)))
            blast_result_file = os.path.join(out_dir, 'Blast', str(count_fasta+1)
                                             + '.' + this_base_name + '_blast_result')
            execute_blast(blast_input, database, blast_result_file, 5, b_word_size, b_evalue, options.resume, log)
            if is_fastg and not options.keep_temp:
                os.remove(blast_input)

            """Parse blast result and produce the immediate result"""
            hsp_infos = parse_blast_xml_result(blast_result_file, log)
            h_site_dicts, q_range_dict = hsp_hits_to_hit_site_dicts(hsp_infos, total_len, log)
            #
            if options.mapping_like:
                out_alignment = alignment_multiple_with_hit_site_dicts(reference_matrix, h_site_dicts, "-", log)
                ref_name, ref_seq = out_alignment[0].pop(0), out_alignment[1].pop(0)
                sorted_alignment = sorted(
                    [(out_alignment[0][x], out_alignment[1][x]) for x in range(len(out_alignment[0]))],
                    key=lambda x: x[1], reverse=True)
                out_alignment = [[x[0] for x in sorted_alignment], [x[1] for x in sorted_alignment], out_alignment[2]]
                out_alignment[0].insert(0, ref_name)
                out_alignment[1].insert(0, ref_seq)
                write_fasta(os.path.join(out_dir, 'Raw', str(count_fasta+1)
                                         + '.' + this_base_name+'_mapping_like.fasta'), out_alignment, True)
            """Reorder queries according to blast positions (core)"""
            if is_fastg:
                edge_connections, kmer = parse_fastg(read_fasta_gb_head(args[count_fasta]), log)
            else:
                edge_connections, kmer = {}, 0
            sequential_name = []
            raw_query_seqs = read_fasta_gb_head(raw_seq_file)
            in_seq_dict = {raw_query_seqs[0][k]: raw_query_seqs[1][k] for k in range(len(raw_query_seqs[0]))}
            # test_name = "EDGE_86_length_31013_cov_41.73211:EDGE_86_length_31013_cov_41.73211"
            # if test_name in in_seq_dict:
            #     print(in_seq_dict[test_name])
            name_cluster = {}
            remove_multiple_hits_per_query(h_site_dicts, q_range_dict, total_len, options.circular_refer,
                                           options.hit_cut_off, options.match_score, options.mismatch_score,
                                           options.gap_score, name_cluster, log)
            to_cluster_info = update_to_cluster(name_cluster)
            remove_multiple_queries_per_hit(h_site_dicts, q_range_dict, total_len, options.circular_refer,
                                            sequential_name, edge_connections, options.query_cut_off,
                                            options.match_score, options.mismatch_score, options.gap_score,
                                            name_cluster, to_cluster_info, in_seq_dict, options.verbose, log)
            to_cluster_info = update_to_cluster(name_cluster)
            if options.separate_result and continuous_conservative_sites:
                check_conservative_continuous(continuous_conservative_sites, h_site_dicts, insertions_in_conservative,
                                              options.conserve_len, log)
            merged_names = {}
            q_range_sets = {}
            merge_hit_site_dicts(h_site_dicts, q_range_dict, q_range_sets, sequential_name, in_seq_dict, total_len,
                                 options, edge_connections, kmer, merged_names, log)
            extend_unmerged(h_site_dicts, q_range_dict, q_range_sets, sequential_name, in_seq_dict, name_cluster,
                            to_cluster_info, total_len, options, merged_names, options.add_gap_disconnect, log)
            if options.mapping_like:
                out_alignment = alignment_multiple_with_hit_site_dicts(reference_matrix, h_site_dicts, "-", log)
                ref_name, ref_seq = out_alignment[0].pop(0), out_alignment[1].pop(0)
                sorted_alignment = sorted(
                    [(out_alignment[0][x], out_alignment[1][x]) for x in range(len(out_alignment[0]))],
                    key=lambda x: x[1], reverse=True)
                out_alignment = [[x[0] for x in sorted_alignment], [x[1] for x in sorted_alignment], out_alignment[2]]
                out_alignment[0].insert(0, ref_name)
                out_alignment[1].insert(0, ref_seq)
                write_fasta(os.path.join(out_dir, 'Raw', str(count_fasta + 1)
                                         + '.' + this_base_name + '_mapping_like.modified.fasta'), out_alignment, True)
            """Add the result of the sequence to total alignment"""
            if options.no_final_aligned:
                final_out_alignment[0].append(re.sub(".fasta$", '', re.sub(".fastg$", '',
                                                                           os.path.basename(args[count_fasta]))))
                if options.separate_result and continuous_conservative_sites:
                    final_out_alignment[1].append(
                        hit_site_dicts_to_sequence_mark_conservative(reference_matrix, h_site_dicts, "N",
                                                                     continuous_conservative_sites, log))
                else:
                    final_out_alignment[1].append(hit_site_dicts_to_sequence(reference_matrix, h_site_dicts, "N", log))
            else:
                combine_site_dict((re.sub(".fasta$", '', re.sub(".fastg$", '', os.path.basename(args[count_fasta]))),
                                   h_site_dicts), total_site_dict, log)
            del q_range_dict, h_site_dicts
            log.info("Analysing " + args[count_fasta] + ' finished.\n')
            if not is_fastg and options.min_repeat and modified and not options.keep_temp:
                os.remove(raw_seq_file)

        log.info("Generating final result ...")
        mafft_commands = []
        fst_aligned = []
        if not options.no_final_aligned:
            final_out_alignment = alignment_multiple_with_hit_site_dicts(reference_matrix, total_site_dict, "N", log)
            # insertions_in_conservative
            if options.separate_result and continuous_conservative_sites:
                add_info_to_cc_sites(continuous_conservative_sites, final_out_alignment[1][0])
                groups = get_groups(continuous_conservative_sites, insertions_in_conservative)
                range_of_index = get_range(continuous_conservative_sites, groups, len(final_out_alignment[1][0]), 0)
                code_length = len(str(len(range_of_index)))
                for file_num in range(len(range_of_index)):
                    this_alignment = [final_out_alignment[0], [], final_out_alignment[2]]
                    for seq_count in range(len(this_alignment[0])):
                        from_id, to_id = range_of_index[file_num]
                        this_alignment[1].append(final_out_alignment[1][seq_count][from_id: to_id])
                    fst = 'Alignment_' + constant_len_number(file_num+1, code_length) + '.fasta'
                    fst_aligned.append(fst[:-5] + "mafft.fasta")
                    write_fasta(os.path.join(out_dir, fst), this_alignment, True)
                    mafft_commands.append("mafft --auto --maxiterate 1000 --thread 4 " + fst + " > " + fst_aligned[-1])
            else:
                write_fasta(os.path.join(out_dir, 'Alignment.fasta'), final_out_alignment, True)
                mafft_commands.append("mafft --auto --maxiterate 1000 --thread 8 "
                                      "Alignment.fasta >  Alignment.mafft.fasta")
        else:
            if options.separate_result and continuous_conservative_sites:
                groups = get_groups(continuous_conservative_sites, insertions_in_conservative)
                ranges_of_index = []
                for seq_count in range(len(final_out_alignment[0])):
                    ranges_of_index.append(get_range(continuous_conservative_sites, groups,
                                                     len(final_out_alignment[1][seq_count]), seq_count))
                code_length = len(str(len(ranges_of_index[0])))
                for file_num in range(len(ranges_of_index[0])):
                    this_alignment = [final_out_alignment[0], [], final_out_alignment[2]]
                    for seq_count in range(len(this_alignment[0])):
                        from_id, to_id = ranges_of_index[seq_count][file_num]
                        this_alignment[1].append(final_out_alignment[1][seq_count][from_id: to_id])
                    fst = 'Sequences_' + constant_len_number(file_num+1, code_length) + '.fasta'
                    fst_aligned.append(fst[:-5] + "mafft.fasta")
                    write_fasta(os.path.join(out_dir, fst), this_alignment, True, True)
                    mafft_commands.append("mafft --auto --maxiterate 1000 --thread 4 " + fst + " > " + fst_aligned[-1])
            else:
                write_fasta(os.path.join(out_dir, 'Sequences.fasta'), final_out_alignment, True, True)
                mafft_commands.append("mafft --auto --maxiterate 1000 --thread 8 "
                                      "Sequences.fasta > Alignment.mafft.fasta")
        if options.separate_result and continuous_conservative_sites:
            mafft_commands.append("concatenate_fasta.py "+" ".join(fst_aligned)+" -o Alignment.mafft.fasta")
        open(os.path.join(out_dir, "mafft.sh"), 'w').writelines([this_line + '\n' for this_line in mafft_commands])
        os.system("chmod 770 "+os.path.join(out_dir, "mafft.sh"))
    except:
        log.exception("")
        log = simple_log(log, out_dir)
        log.info("\nTotal cost " + str(time.time() - time0))
        log.info("If you find bugs, "
                 "Please run the file with \"--raw\" again, and then send me the Raw folder and Map.log.txt via email.")
    else:
        log = simple_log(log, out_dir)
        log.info("\nTotal cost " + str(time.time() - time0))
        log.info("Thanks you!")
    logging.shutdown()


if __name__ == "__main__":
    main()


"""Copyright 2018 Jianjun Jin

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License."""