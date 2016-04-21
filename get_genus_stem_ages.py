__author__ = 'Kinggerm'
# python2
import dendropy
import os


def main():
    tree_f = raw_input('Input mcc nexus tree:').strip()
    threshold = float(raw_input('Input threshold to cut:'))
    if not (0 < threshold < 1):
        print 'Error: threshold has to be set in (0, 1)!'
        os._exit(0)
    tree = dendropy.Tree.get(path=tree_f, schema='nexus')
    tree.calc_node_ages()
    terminals = tree.leaf_nodes()
    # tree.taxon_namespace
    stems = {}
    for terminal in terminals:
        this_age = terminal.parent_node.age
        stems[str(terminal.taxon).strip("'")] = this_age
    # print list
    print ''
    taxa = stems.keys()
    taxa.sort(key=lambda x:stems[x])
    ages = stems.values()
    max_t_l = max([len(str(x)) for x in taxa])+4
    max_a_l = max([len(str(x)) for x in ages])+1
    for taxon in taxa:
        print repr(str(taxon)).ljust(max_t_l).replace("'", ''), repr(stems[taxon]).ljust(max_a_l)

    # print cut
    # ages = list(set(ages))
    ages.sort()
    print '\nCut age at:', ages[int(len(ages)*threshold)]


if __name__ == '__main__':
    main()

