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
import networkx as nx
import matplotlib
from operator import itemgetter
import random

from networkx.algorithms.simple_paths import all_simple_paths
from networkx.algorithms.tree.mst import ALGORITHMS
random.seed(9001)
from random import randint
import statistics

__author__ = "Mamadou DANSOKHO"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
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
                        default=22, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    with open(fastq_file,"r+") as file:
        for line in file :
            if "@" not in line and "+" not in line and "J" not in line:
                yield line.rstrip()





def cut_kmer(read, kmer_size):
    for i in range(len(read)-kmer_size +1):
        kmer = read[i:i+kmer_size]
        yield kmer




def build_kmer_dict(fastq_file, kmer_size):
    kmer_dict = {}
    for seq in read_fastq(fastq_file):
        for kmer in cut_kmer(seq, kmer_size):
            if kmer in kmer_dict:
                kmer_dict[kmer]+=1
            else:
                kmer_dict.update({kmer:1})
    return kmer_dict
            


def build_graph(kmer_dict):
    graph  = nx.DiGraph()
    for kmer in kmer_dict:
        node1 = kmer[:-1]
        node2 = kmer[1:]
        graph.add_edge(node1,node2, weight = kmer_dict[kmer])
    return graph 


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """nodes_list = list(graph.nodes)
    path_list = list(graph.edge)
    for nodes in nodes_list :"""

def std(data):
    pass


def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    pass

def path_average_weight(graph, path):
    pass

def solve_bubble(graph, ancestor_node, descendant_node):
    pass

def simplify_bubbles(graph):
    pass

def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    nodes_list = []
    for node in graph:
        pred_node = [node for node in graph.predecessors(node)]  
        if len(pred_node)==0:
            nodes_list.append(node)
    return nodes_list

def get_sink_nodes(graph):
    node_list = []
    for node in graph:
        succ_node = [node for node in graph.successors(node)]
        if len(succ_node) == 0:
            node_list.append(node)
    return node_list


def get_contig_from_path(path):
    contig = ""
    contig+=path[0]
    for node in path[1:]:
        contig+=node[-1]
    return contig

def get_contigs(graph, starting_nodes, ending_nodes):
    contigs  = []

    for starting_node in starting_nodes:
        for ending_node in ending_nodes:
            tmp_contigs = [get_contig_from_path(path) for path in nx.all_simple_paths(graph, starting_node, ending_node)]
            contigs+=[(contig, len(contig)) for contig in tmp_contigs]

    
    return contigs


def save_contigs(contigs_list, output_file):

    pass


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def draw_graph(graph, graphimg_file):
    """Draw the graph
    """                                    
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


def save_graph(graph, graph_file):
    """Save the graph with pickle
    """
    with open(graph_file, "wt") as save:
            pickle.dump(graph, save)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)
    # Save the graph in file
    # if args.graph_file:
    #     save_graph(graph, args.graph_file)





if __name__ == '__main__':
    main()