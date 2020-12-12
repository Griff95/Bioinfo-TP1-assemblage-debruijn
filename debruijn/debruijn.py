#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import statistics
import itertools
import random
random.seed(9001)
from random import randint
from operator import itemgetter

import matplotlib
import tkinter
import matplotlib.pyplot as plt
import networkx as nx



__author__ = "Pierre Côte"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Pierre Côte"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Pierre Côte"
__email__ = "cotepierre@eisti.eu"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()


def read_fastq(fastq_file):
    with open(fastq_file) as file:
        for line in file:
            seq = next(file)
            # delete newline char
            if seq[-1] == '\n':
                seq = seq[:-1]
            next(file)
            next(file)
            yield seq


def cut_kmer(read, kmer_size):
    for i in range(len(read) - kmer_size +1):
        yield read[i:i+kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    build = dict()
    for seq in read_fastq(fastq_file):
        for kmer in cut_kmer(seq, kmer_size):
            if kmer in build:
                build[kmer] += 1
            else:
                build[kmer] = 1
    return build


def build_graph(kmer_dict):
    G = nx.DiGraph()
    for key in kmer_dict.keys():
        prefix = key[:-1]
        suffix = key[1:]
        # print(prefix, suffix)
        G.add_edge(prefix, suffix, weight=kmer_dict[key])
    return G


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    for path in path_list:
        for n1, n2 in itertools.zip_longest(path[:-1], path[1:]):
            graph.remove_edge(n1, n2)
    if delete_entry_node:
        for entry in map(lambda x: x[0], path_list):
            graph.remove_node(entry)
    if delete_sink_node:
        for sink in map(lambda x: x[-1], path_list):
            graph.remove_node(sink)
    remove = [node for node, degree in dict(graph.degree()).items() if degree == 0]
    graph.remove_nodes_from(remove)
    return graph

def std(data):
    return statistics.stdev(data)


def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
    path_data = list(zip(weight_avg_list, path_length, range(len(path_list))))
    ranked = sorted(path_data, key = lambda x: (-x[0], -x[1]))
    best = ranked[0]
    candidates = [x for x in ranked if x[0] == best[0] and x[1] == best[1]]
    r = randint(0, len(candidates)-1)
    best_path_index = candidates[r][2]
    to_remove = path_list[:best_path_index] + path_list[best_path_index+1:]
    return remove_paths(graph, to_remove, delete_entry_node, delete_sink_node)

def path_average_weight(graph, path):
    weight_sum = 0
    count = 0
    for n1, n2 in itertools.zip_longest(path[:-1], path[1:]):
        weight_sum += graph[n1][n2]['weight']
        count += 1
    return weight_sum/count

def solve_bubble(graph, ancestor_node, descendant_node):
    path_list = []
    path_length = []
    path_avg_weight = []
    for path in nx.all_simple_paths(graph, source=ancestor_node, target=descendant_node):
        path_list.append(path)
        path_length.append(len(path))
        path_avg_weight.append(path_average_weight(graph, path))
    return select_best_path(graph, path_list, path_length, path_avg_weight)


def simplify_bubbles(graph):
    graph_copy = graph.copy()
    for node in graph.nodes():
        if graph.in_degree(node) > 1:
            predecessors = list(graph.predecessors(node))
            pairs = list(itertools.zip_longest(predecessors[:-1], predecessors[1:]))
            ancestors = nx.all_pairs_lowest_common_ancestor(graph, pairs)
            for ((n1, n2), lca) in ancestors:
                solve_bubble(graph_copy, lca, node)
    return graph_copy

def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass



def get_starting_nodes(graph):
    starting_nodes = []
    for node in graph.nodes:
        if len(list(graph.predecessors(node))) == 0:
            starting_nodes.append(node)
    return starting_nodes

def get_sink_nodes(graph):
    ending_nodes = []
    for node in graph.nodes:
        if len(list(graph.successors(node))) == 0:
            ending_nodes.append(node)
    return ending_nodes

def get_contigs(graph, starting_nodes, ending_nodes):
    contigs = []
    for start in starting_nodes:
        for end in ending_nodes:
            for path in nx.all_simple_paths(graph, source=start, target=end):
                contig = path[0]
                for n in path[1:]:
                    contig += n[-1]
                contig_size = len(contig)
                contigs.append((contig, contig_size))
    return contigs


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def save_contigs(contigs_list, output_file):
    with open(output_file, 'w') as f:
        for i in range(len(contigs_list)):
            seq, size = contigs_list[i]
            f.write(">"+str(i+1) + " len=" + str(size)+'\n')
            f.write(fill(seq))
            f.write('\n')

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    dict = build_kmer_dict(args.fastq_file, 21)
    # print(dict)
    graph = build_graph(dict)
    # plt.figure()
    # nx.draw(graph, with_labels=True, font_weight='bold')
    # plt.show()
    # print(graph.nodes)
    # print(graph.edges)
    # print(graph.degree)
    starts = get_starting_nodes(graph)
    ends = get_sink_nodes(graph)
    print(starts)
    print(ends)
    contigs = get_contigs(graph, starts, ends)
    print(*contigs, sep='\n')
    save_contigs(contigs, "test")
    print(list(cut_kmer("TCAGA", 3)))


if __name__ == '__main__':
    main()
