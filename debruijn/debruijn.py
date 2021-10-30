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
__credits__ = ["Mamadou DANSOKHO"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Mamadou DANSOKHO"
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
    for path in path_list:
        if delete_entry_node is True and delete_sink_node is True:
            graph.remove_nodes_from(path)
        elif delete_entry_node is True and delete_sink_node is False:
            graph.remove_nodes_from(path[:-1])
        elif delete_entry_node is False and delete_sink_node is True:
            graph.remove_nodes_from(path[1:])
        else:
            graph.remove_nodes_from(path[1:-1])
    return graph


def std(data):
    if len(data) == 1:
        return 0
    return statistics.stdev(data)


def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
    if std(weight_avg_list) > 0:
        best_index = weight_avg_list.index(max(weight_avg_list))
        del path_list[best_index]
        return remove_paths(graph,path_list,delete_entry_node,delete_sink_node)
    elif std(path_length) > 0:
        best_index = path_length.index(max(path_length))
        del path_list[best_index]
        return remove_paths(graph,path_list,delete_entry_node,delete_sink_node)
    else:
        best_index = randint(0,len(path_list)-1)
        del path_list[best_index]
        return remove_paths(graph,path_list,delete_entry_node,delete_sink_node)

def path_average_weight(graph, path):

    return statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])

def solve_bubble(graph, ancestor_node, descendant_node):
    path_length = []
    path_list= [p for p in nx.all_simple_paths(graph,ancestor_node,descendant_node)]
    weight_avg_list = []
    for path in path_list:
        path_length.append(len(path))
        weight_avg_list.append(path_average_weight(graph,path))
    return select_best_path(graph, path_list, path_length, weight_avg_list)
    

def simplify_bubbles(graph):
    bubble = False
    for node in graph.nodes():
        predecessors_list = [n for n in graph.predecessors(node)]
        if len(predecessors_list) > 1:
            for i in range(0,len(predecessors_list)-1):
                ancestor_node = nx.lowest_common_ancestor(graph,predecessors_list[i],predecessors_list[i+1])
                if ancestor_node != None:
                    bubble = True
                    break
            if bubble == True:
                break
    if bubble == True:
        graph = simplify_bubbles(solve_bubble(graph,ancestor_node,node))
    return graph


def solve_entry_tips(graph, starting_nodes):
    for node in graph.nodes():
        path_list = []
        path_length = []
        weight_avg_list = []
        tip_exist = False
        predecessors_nodes_list = [n for n in graph.predecessors(node)]
        if len(predecessors_nodes_list) > 1:
            for starting_node in starting_nodes:
                if nx.has_path(graph, starting_node, node):
                    path = [n for n in nx.all_simple_paths(graph, starting_node, node)]
                    path_list.append(path[0])
                    path_length.append(len(path[0]))
                    weight_avg_list.append(path_average_weight(graph,path[0]))
                    tip_exist = True
        if tip_exist is True:
            break
    if tip_exist is True:
        graph = select_best_path(graph, path_list, path_length, weight_avg_list, delete_entry_node=True)
        solve_entry_tips(graph,starting_nodes)
    return graph

def solve_out_tips(graph, ending_nodes):
    for node in graph.nodes():
        path_list = []
        path_length = []
        weight_avg_list = []
        tip_exist = False
        successors_nodes_list = [n for n in graph.successors(node)]
        if len(successors_nodes_list) > 1:
            for ending_node in ending_nodes:
                if nx.has_path(graph, node, ending_node):
                    path = [n for n in nx.all_simple_paths(graph, node, ending_node)]
                    path_list.append(path[0])
                    path_length.append(len(path[0]))
                    weight_avg_list.append(path_average_weight(graph,path[0]))
                    tip_exist = True
        if tip_exist is True:
            break
    if tip_exist is True:
        graph = select_best_path(graph, path_list, path_length, weight_avg_list, delete_sink_node=True)
        solve_out_tips(graph,ending_nodes)
    return graph

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

    with open(output_file,"w+", newline="\n") as file_out:
        for index, contigs in enumerate(contigs_list):
            contig = contigs[0]
            file_out.write(">contig_{} len={}\n".format(index,len(contig)))
            file_out.write("{}\n".format(fill(contig,80)))

    return file_out


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
    """"
    # Get arguments
    args = get_arguments()
    #****************************1) reading file and building graph
    kmer_dict =build_kmer_dict(args.fastq_file,args.kmer_size)
    graph = build_graph(kmer_dict)
    #******************************2) Solving bulles
    graph = simplify_bubbles(graph)
    #*****************************3)solving entry and output tips
    graph = solve_entry_tips(graph, get_starting_nodes(graph))
    graph = solve_out_tips(graph, get_sink_nodes(graph))
    #******************************4) writing contigs
    save_contigs(get_contigs(graph, get_starting_nodes(graph),get_sink_nodes(graph)), args.saved_contigs_output_file)
    """



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
