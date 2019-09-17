#!/usr/bin/env python
import dendropy
import sys
import os
from optparse import OptionParser
from copy import deepcopy
from itertools import combinations


example_tree_1_str = \
    "(Trema_tomentosa, " \
    "(((Trema_angustifolia, Trema_nitida)95, (Trema_aspera, Trema_discolor)87)100, "\
    "(Parasponia_rugosa, Parasponia_parviflora)65)65, " \
    "(Cannabis_sativa, Humulus_lupulus)100);\n"
example_tree_2_str = "((Parasponia_parviflora, Parasponia_rugosa)69, " \
    "(((Trema_angustifolia, Trema_nitida)95, (Trema_aspera, Trema_discolor)90)100, "\
    "Trema_tomentosa)33, " \
    "(Cannabis_sativa, Humulus_lupulus)100);\n"
criteria_tab_str = """\
label\tgenus\tbig_clade
Trema_tomentosa\tTrema\tTrema-Parasponia clade
Trema_angustifolia\tTrema\tTrema-Parasponia clade
Trema_nitida\tTrema\tTrema-Parasponia clade
Trema_aspera\tTrema\tTrema-Parasponia clade
Trema_discolor\tTrema\tTrema-Parasponia clade
Parasponia_rugosa\tParasponia\tTrema-Parasponia clade
Parasponia_parviflora\tParasponia\tTrema-Parasponia clade
Cannabis_sativa\tCannabis\tCannabis-Humulus clade
Humulus_lupulus\tHumulus\tCannabis-Humulus clade
"""


def get_options():
    parser = OptionParser(usage="check_monophyletic.py -c criteria_tab a_list_of_tree_files")
    parser.add_option("-c", dest="criteria_tab",
                      help="criteria tab file with head: label\\tfamily\\tclade ...")
    parser.add_option("-r", dest="rooted", default=False, action="store_true",
                      help="input tree(s) is(are) rooted. [default: %default]")
    parser.add_option("--support-cutoff", dest="support_cutoff", type="float", default=None,
                      help="Cutoff for collapsing low supported nodes. [default: %default]")
    parser.add_option("--clarify", dest="clarify", default=False, action="store_true",
                      help="distinguish unresolved from paraphyletic. [default: %default]")
    parser.add_option("--off-verbose", dest="verbose", default=True, action="store_false",
                      help="turn off verbose. [default: on]")
    parser.add_option("--example", dest="generate_example", default=False, action="store_true",
                      help="Generate local example. [default: off]")
    options, argv = parser.parse_args()
    if options.generate_example:
        tree1 = "check_monophyletic--tree1.tre"
        tree2 = "check_monophyletic--tree2.tre"
        criteria = "check_monophyletic--criteria.tab"
        open(tree1, "w").write(example_tree_1_str)
        open(tree2, "w").write(example_tree_2_str)
        open(criteria, "w").write(criteria_tab_str)
        sys.stdout.write("\nExample files generated: \n" +
                         os.path.realpath(tree1) + "\n" +
                         os.path.realpath(tree2) + "\n" +
                         os.path.realpath(criteria) + "\n" +
                         "\n"
                         "# Please run simple checking by:\n"
                         "\x1b[1;36m"
                         "check_monophyletic.py -c check_monophyletic--criteria.tab check_monophyletic--tree*tre"
                         "\x1b[0m\n"
                         "\n"
                         "# If you want to export the screen-printed result to a tab file:\n"
                         "\x1b[1;36m"
                         "check_monophyletic.py -c check_monophyletic--criteria.tab check_monophyletic--tree*tre"
                         "\x1b[0m"
                         " > result.tab\n"
                         "\n"
                         "# If you want to extract tree(s) with \x1b[1;31mmonophyletic Trema genus\x1b[0m, "
                         "please run filtering by:\n"
                         "mkdir \x1b[1;31mTrema-monophyly\x1b[0m\n"
                         "cp `"
                         "\x1b[1;36m"
                         "check_monophyletic.py -c check_monophyletic--criteria.tab check_monophyletic--tree*tre"
                         "\x1b[0m"
                         "| grep \"\x1b[1;31mgenus\\tTrema\\tmonophyletic\x1b[0m\" "
                         "| cut -f 2` \x1b[1;31mTrema-monophyly\x1b[0m\n\n")
                         # no color version
                         # "# If you want to extract tree(s) with monophyletic Trema genus, please run filtering by:\n"
                         # "mkdir Trema-monophyly\n"
                         # "cp `
                         # check_monophyletic.py -c check_monophyletic--criteria.tab check_monophyletic--tree*tre "
                         # "|grep \"genus\\tTrema\\tmonophyletic\" | cut -f 2` Trema-monophyly\n\n")
        sys.exit()
    elif not (options.criteria_tab and len(argv)):
        parser.print_help()
        sys.stdout.write("Insufficient arguments!\n")
        sys.exit()
    return options, argv


def get_tree(tree_file_or_str_or_obj, name_space=None, schema="newick", edge_label=False, **kwargs):
    if not name_space:
        name_space = dendropy.TaxonNamespace()
    if type(tree_file_or_str_or_obj) == dendropy.Tree:
        tree_a = deepcopy(tree_file_or_str_or_obj)
        tree_a.taxon_namespace = name_space
    elif type(tree_file_or_str_or_obj) == str:
        if os.path.isfile(tree_file_or_str_or_obj):
            tree_a = dendropy.Tree.get(path=tree_file_or_str_or_obj, schema=schema, taxon_namespace=name_space,
                                       is_assign_internal_labels_to_edges=edge_label,
                                       preserve_underscores=True, **kwargs)
        else:
            try:
                tree_a = dendropy.Tree.get_from_string(tree_file_or_str_or_obj, schema=schema,
                                                       taxon_namespace=name_space,
                                                       is_assign_internal_labels_to_edges=edge_label,
                                                       preserve_underscores=True, **kwargs)
            except ValueError:
                raise ValueError("Error: " + str(tree_file_or_str_or_obj) +
                                 " is neither a file or a " + schema + "string\n")
    else:
        sys.stdout.write("Error: parsing tree " + str(tree_file_or_str_or_obj))
        sys.exit()
    return tree_a


def collapse_nodes(in_tree, support_cutoff):
    tree = deepcopy(in_tree)
    for node in tree.postorder_node_iter():
        if not node.is_leaf() and node != tree.seed_node:
            try:
                this_label = float(node.edge.label)
            except (ValueError, TypeError):
                sys.stdout.write("\tWarning: supporting information not found in node! Skip!\n")
                return in_tree
            else:
                if this_label < support_cutoff:
                    node.edge.collapse()
    return tree


def read_criteria(criteria_f):
    tab_context = [line.strip("\n").split("\t") for line in open(criteria_f) if line.strip()]
    head = tab_context.pop(0)
    if "label" not in head:
        sys.stdout.write("Error: \"label\" not found in criteria head!\n")
        sys.exit()
    context_dict = {}
    for record in tab_context:
        context_dict[record[0]] = {}
        for j, key_c in enumerate(head[1:]):
            context_dict[record[0]][key_c] = record[j + 1]
    return context_dict, head[1:]


def make_criteria_set(original_dict, keys, limit_set, verbose=True):
    criteria_set = {k: {} for k in keys}
    for label in sorted(original_dict):
        if label not in limit_set:
            if verbose:
                sys.stdout.write("\tWarning: " + label + " in the criteria not found in the tree!\n")
        else:
            for k in keys:
                value = original_dict[label][k]
                if value not in criteria_set[k]:
                    criteria_set[k][value] = {label}
                else:
                    criteria_set[k][value].add(label)
    for k in keys:
        for value in criteria_set[k]:
            criteria_set[k][value] = tuple(sorted(criteria_set[k][value]))
    return criteria_set


def slim_tree(tree, taxon_set, limit_set, verbose=True):
    tips_to_drop = []
    for label in sorted(taxon_set):
        if label not in limit_set:
            tips_to_drop.append(label)
            if verbose:
                sys.stdout.write("\tWarning: " + label + " in the tree not found in the criteria!\n")
    if len(tips_to_drop) == len(limit_set):
        raise Exception("All tips removed!")
    else:
        tree.prune_taxa_with_labels(tips_to_drop)


def main():
    options, argv = get_options()
    criteria, keys = read_criteria(criteria_f=options.criteria_tab)
    sys.stdout.write("TreeId\tTreeName\tGroupLevel\tGroupName\tGroupStatus\n")
    for go_to, this_tre_f in enumerate(argv):
        sys.stdout.write("\t\t\t\t\n")
        tree_id = str(go_to + 1)
        this_tre = get_tree(this_tre_f, edge_label=True)
        if options.support_cutoff is not None:
            this_tre = collapse_nodes(this_tre, options.support_cutoff)
        tree_taxon_set = set([t.label for t in this_tre.taxon_namespace])
        criteria_set = make_criteria_set(criteria, keys, tree_taxon_set, options.verbose)
        slim_tree(this_tre, tree_taxon_set, criteria, options.verbose)

        # make searching table; straightforward, but not efficient;
        # mrca bipartition comparison would be faster
        mono_partition = set()
        unresolved_partition = set()
        for node in this_tre.postorder_node_iter():
            if not node.is_leaf():
                directed = options.rooted or node != this_tre.seed_node
                children = node.child_nodes()
                for child in children:
                    this_partition = set([leaf.taxon.label for leaf in child.leaf_iter()])
                    mono_partition.add(tuple(sorted(this_partition)))
                    if not options.rooted:
                        mono_partition.add(tuple(sorted(tree_taxon_set - this_partition)))
                if options.clarify and len(children) > 2:
                    for comb_num in range(2, int(len(children) / (2 - directed)) + (not directed)):
                        for assumptive_clade in combinations(children, comb_num):
                            this_partition = set()
                            for child in assumptive_clade:
                                for leaf in child.leaf_iter():
                                    this_partition.add(leaf.taxon.label)
                            unresolved_partition.add(tuple(sorted(this_partition)))
                            if not options.rooted:
                                unresolved_partition.add(tuple(sorted(tree_taxon_set - this_partition)))

        for k in sorted(criteria_set):
            for value in sorted(criteria_set[k]):
                if criteria_set[k][value] not in mono_partition:
                    if options.clarify:
                        if criteria_set[k][value] in unresolved_partition:
                            sys.stdout.write(tree_id + "\t" + this_tre_f + "\t" + k + "\t" + value + "\t" + "unresolved\n")
                        else:
                            sys.stdout.write(tree_id + "\t" + this_tre_f + "\t" + k + "\t" + value + "\t" + "para/poly-phyletic\n")
                    else:
                        sys.stdout.write(tree_id + "\t" + this_tre_f + "\t" + k + "\t" + value + "\t" + "para/poly-phyletic/unresolved\n")
                elif options.verbose:
                    sys.stdout.write(tree_id + "\t" + this_tre_f + "\t" + k + "\t" + value + "\t" + "monophyletic\n")
        # if go_to < len(argv) - 1:
        #     sys.stdout.write("\t\t\t\n")


if __name__ == '__main__':
    main()
