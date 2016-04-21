__author__ = 'Kinggerm'
# python2
import dendropy
import os
import sys
# import numpy as np
# from Tkinter import *
import matplotlib.pyplot as plt
# from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
# from matplotlib.figure import Figure


# def drawCurve_interactively():
def draw_curves(x_dots, y_dots, replicates, out_f):
    fig = plt.figure(figsize=(16, 8))
    curves = fig.add_subplot(111)
    i = 0
    for x_rep, y_rep in replicates:
        i += 1
        curves.plot(x_rep, y_rep, 'k-', linewidth=2, color='#E57C5E')
        sys.stdout.write(str(i)+'\b'*len(str(i)))
        sys.stdout.flush()
    curves.plot(x_dots, y_dots, 'k-', linewidth=2, color='#DC4D38')
    curves.set_title('Sliding window analysis')
    curves.set_xlabel('Ma')
    curves.set_xlim(-90, 0)
    curves.set_ylim(0, 1.5)
    fig.savefig(out_f, bbox_inches='tight')


def cal_times_rates(tree, window_size, step):
    origin_t = tree.max_distance_from_root()
    this_time = window_size
    times = []
    rates = []
    while this_time <= origin_t:
        start_l = tree.num_lineages_at(this_time-window_size)
        if start_l:
            times.append(this_time-origin_t)
            end_l = tree.num_lineages_at(this_time)
            rates.append((end_l-start_l)/float(start_l))
        this_time += step
    times = [x-window_size*0.5 for x in times]
    return times, rates


def write_csv_result(out_f, replicates):
    out_csv = open(out_f, 'wb')
    for replicate in replicates:
        out_csv.write('\t'.join([str(x) for x in replicate[0]])+'\n')
        out_csv.write('\t'.join([str(x) for x in replicate[1]])+'\n')
    out_csv.close()


def main():
    tree_f = raw_input('Input mcc nexus tree:').strip()
    trees_f = raw_input('Input nexus trees (skip):').strip()
    step = float(raw_input('Input step:'))
    window_size = int(raw_input('Input window size:'))
    tree = dendropy.Tree.get(path=tree_f, schema='nexus')
    times, rates = cal_times_rates(tree, window_size, step)
    replicates = []
    i = 1
    if trees_f.strip():
        namespace = dendropy.TaxonNamespace()
        tree_yielder = dendropy.Tree.yield_from_files(files=[trees_f], schema='nexus', taxon_namespace=namespace, store_tree_weights=True, preserve_underscores=True, ignore_unrecognized_keyword_arguments=True)
        i = 0
        for tree in tree_yielder:
            i += 1
            this_time, this_rate = cal_times_rates(tree, window_size, step)
            replicates.append((this_time, this_rate))
            sys.stdout.write(str(i)+'\b'*len(str(i)))
            sys.stdout.flush()
        write_csv_result(trees_f+'.s'+str(step)+'.w'+str(window_size)+'.t'+str(i)+'.csv', replicates)
    write_csv_result(tree_f+'.s'+str(step)+'.w'+str(window_size)+'.t'+str(1)+'.csv', [(times, rates)])
    draw_curves(times, rates, replicates, tree_f + '.sliding.s'+str(step)+'.w'+str(window_size)+'.t'+str(i)+'.pdf')

if __name__ == '__main__':
    main()
