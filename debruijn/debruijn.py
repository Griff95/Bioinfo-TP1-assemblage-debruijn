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
import itertools
import statistics
import random
from random import randint
from collections import defaultdict
import networkx as nx
random.seed(9001)

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
    """Generator that yield each line sequence of fastq file"""
    with open(fastq_file) as file:
        lines = file.read().splitlines()
    desired_lines = lines[1::4]
    for line in desired_lines:
        yield line


def cut_kmer(read, kmer_size):
    """Generator of kmer of size kmer_size from sequence read"""
    for i in range(len(read) - kmer_size +1):
        yield read[i:i+kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    """Return a dictionary of kmer from fastq_file of size kmer_size"""
    build = dict()
    for seq in read_fastq(fastq_file):
        for kmer in cut_kmer(seq, kmer_size):
            if kmer in build:
                build[kmer] += 1
            else:
                build[kmer] = 1
    return build


def build_graph(kmer_dict):
    "Return networkx weighted DiGraph based on kmer dict"
    graph = nx.DiGraph()
    for key in kmer_dict.keys():
        prefix = key[:-1]
        suffix = key[1:]
        # print(prefix, suffix)
        graph.add_edge(prefix, suffix, weight=kmer_dict[key])
    return graph


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """Return graph with removed path from path_list"""
    for path in path_list:
        for node_1, node_2 in itertools.zip_longest(path[:-1], path[1:]):
            graph.remove_edge(node_1, node_2)
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
    """Return standard deviation of data list"""
    return statistics.stdev(data)


def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
    """Return graph with selected best path of path_list, remove others"""
    path_data = list(zip(weight_avg_list, path_length, range(len(path_list))))
    ranked = sorted(path_data, key = lambda x: (-x[0], -x[1]))
    best = ranked[0]
    candidates = [x for x in ranked if x[0] == best[0] and x[1] == best[1]]
    r_int = randint(0, len(candidates)-1)
    best_path_index = candidates[r_int][2]
    to_remove = path_list[:best_path_index] + path_list[best_path_index+1:]
    return remove_paths(graph, to_remove, delete_entry_node, delete_sink_node)


def path_average_weight(graph, path):
    """Return average weight of all edge weight in path"""
    weight_sum = 0
    count = 0
    for node_1, node_2 in itertools.zip_longest(path[:-1], path[1:]):
        weight_sum += graph[node_1][node_2]['weight']
        count += 1
    return weight_sum/count


def solve_bubble(graph, ancestor_node, descendant_node):
    """Return graph without any bubbles between ancestor and descendant nodes"""
    path_list = []
    path_length = []
    path_avg_weight = []
    for path in nx.all_simple_paths(graph, ancestor_node, descendant_node):
        path_list.append(path)
        path_length.append(len(path))
        path_avg_weight.append(path_average_weight(graph, path))
    return select_best_path(graph, path_list, path_length, path_avg_weight)


def simplify_bubbles(graph):
    """Return a graph copy without any bubbles"""
    graph_copy = graph.copy()
    for node in graph.nodes():
        if graph.in_degree(node) > 1:
            predecessors = list(graph.predecessors(node))
            pairs = list(itertools.combinations(predecessors, 2))
            ancestors = nx.all_pairs_lowest_common_ancestor(graph, pairs)
            for (_, lca) in ancestors:
                solve_bubble(graph_copy, lca, node)
    return graph_copy


def solve_entry_tips(graph, starting_nodes):
    """Return graph without tips of a starting node"""
    intersections_dict = defaultdict(list)
    # get first met intersection points for each starting node
    # find successor node with 2 or more predecessors
    for starting in starting_nodes:
        intersection = starting
        while graph.out_degree(intersection)==1:
            intersection = list(graph.successors(intersection))[0]
            if graph.in_degree(intersection) > 1:
                break
        # check that we get an intersection point
        # save starting point in the list for key 'intersection' in dict
        if graph.in_degree(intersection) > 1:
            intersections_dict[intersection].append(starting)

    # for each intersection point (as key in dict)
    # get all candidates (= starting points that share the same intersection)
    # select best path for all combinations of candidates with the intersection point
    for inter in intersections_dict.keys():
        candidates = intersections_dict[inter]
        if len(candidates) > 1:
            for (candidate_1, candidate_2) in itertools.combinations(candidates, 2):
                path1 = nx.shortest_path(graph, candidate_1, inter)
                path2 = nx.shortest_path(graph, candidate_2, inter)
                path_list = [path1, path2]
                path_weight = [path_average_weight(graph, p) for p in path_list]
                path_length = [len(p) for p in path_list]
                graph = select_best_path(graph, path_list, path_length, path_weight)
    return graph


def solve_out_tips(graph, ending_nodes):
    """Return graph without tips of an ending node"""
    intersections_dict = defaultdict(list)
    # get first met intersection points for each ending node
    # find predecessors node with 2 or more successors
    for ending in ending_nodes:
        intersection = ending
        while graph.in_degree(intersection)==1:
            intersection = list(graph.predecessors(intersection))[0]
            if graph.out_degree(intersection) > 1:
                break
        # check that we get an intersection point
        # save ending point in the list for key 'intersection' in dict
        if graph.out_degree(intersection) > 1:
            intersections_dict[intersection].append(ending)

    # for each intersection point (as key in dict)
    # get all candidates (= ending points that share the same intersection)
    # select best path for all combinations of candidates with the intersection point
    for inter in intersections_dict.keys():
        candidates = intersections_dict[inter]
        if len(candidates) > 1:
            for (candidate_1, candidate_2) in itertools.combinations(candidates, 2):
                path1 = nx.shortest_path(graph, inter, candidate_1)
                path2 = nx.shortest_path(graph, inter, candidate_2)
                path_list = [path1, path2]
                path_weight = [path_average_weight(graph, p) for p in path_list]
                path_length = [len(p) for p in path_list]
                graph = select_best_path(graph, path_list, path_length, path_weight)
    return graph



def get_starting_nodes(graph):
    """Return list of starting nodes (no predecessors) from graph"""
    starting_nodes = []
    for node in graph.nodes:
        if len(list(graph.predecessors(node))) == 0:
            starting_nodes.append(node)
    return starting_nodes


def get_sink_nodes(graph):
    """Return list of ending nodes (no successor) from graph"""
    ending_nodes = []
    for node in graph.nodes:
        if len(list(graph.successors(node))) == 0:
            ending_nodes.append(node)
    return ending_nodes


def get_contigs(graph, starting_nodes, ending_nodes):
    """Return list of contig in graph
    from a node of starting_nodes list to a node of ending_nodes list"""
    contigs = []
    for start in starting_nodes:
        for end in ending_nodes:
            for path in nx.all_simple_paths(graph, source=start, target=end):
                contig = path[0]
                for node in path[1:]:
                    contig += node[-1]
                contig_size = len(contig)
                contigs.append((contig, contig_size))
    return contigs


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def save_contigs(contigs_list, output_file):
    """Write contigs_list in output_file as fasta format"""
    with open(output_file, 'w') as file:
        for i, contig in enumerate(contigs_list):
            seq, size = contig
            file.write(">contig_"+str(i) + " len=" + str(size)+'\n')
            file.write(fill(seq))
            file.write('\n')




# def draw_graph(graph, graphimg_file):
#     """Draw the graph
#     """
#     fig, ax = plt.subplots()
#     elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
#     #print(elarge)
#     esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
#     #print(elarge)
#     # Draw the graph with networkx
#     #pos=nx.spring_layout(graph)
#     pos = nx.random_layout(graph)
#     nx.draw_networkx_nodes(graph, pos, node_size=6)
#     nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
#     nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5,
#                            edge_color='b', style='dashed')
#     #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
#     # save image
#     plt.savefig(graphimg_file)

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    # build graph
    kmer_dict = build_kmer_dict(args.fastq_file, 21)
    graph = build_graph(kmer_dict)

    # solve bubble
    graph = simplify_bubbles(graph)

    # solve tips
    starts = get_starting_nodes(graph)
    graph = solve_entry_tips(graph, starts)
    ends = get_sink_nodes(graph)
    graph = solve_out_tips(graph, ends)

    # save contigs
    starts = get_starting_nodes(graph)
    ends = get_sink_nodes(graph)
    contigs = get_contigs(graph, starts, ends)
    save_contigs(contigs, "contig")
    # draw_graph(graph, 'graph.png')


if __name__ == '__main__':
    main()
