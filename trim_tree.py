"""
Trim a tree at a given phylogenetic leaf e.g. phylum and write a new leaf list
"""

import os
import sys
import argparse
from ete3 import Tree

global tag 

def load_jplacer(jpf):
    """
    lo   ad the jplacer file and return the tree
    :param jpf: The jplacer file
    :return: the data structure of the tree
    """

    with open(jpf, 'r') as f:
        data = json.load(f)
    tree = Tree(data['tree'], quoted_node_names=True, format=1)

    return tree

def is_leaf_node(node):
    global tag
    if tag in node.name:
        return True
    else:
        return False


def write_leaves(tree, outputf):
    """
    Write a list of all the leaves, one line per leaf.
    :param tree: the tree
    :param outputf: the file to write
    :return:
    """

    with open(outputf, 'w') as out:
        for n in tree.get_leaves():
            out.write("{}\n".format(n.name))


def write_tree(tree, outputf):
    """
    Write the tree to a file.
    :param tree: The tree to write
    :param outputf: The output filename
    :return:
    """

    tree.write(outfile=outputf, format=1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="trim tree")
    parser.add_argument('-t', help='tree file to trim', required=True)
    parser.add_argument('-p', help='phylogenetic level to trim at', required=True)
    parser.add_argument('-o', help='output tree to write')
    parser.add_argument('-l', help='list of leaves to write')
    parser.add_argument('-j', help='tree is a jplacer file. Default is to assume tree will be a newick file', action='store_true', default=False)
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    if args.j:
        tree = load_jplacer(args.t)
    else:
        tree = Tree(args.t, quoted_node_names=True, format=1)

    tag = "r_{}".format(args.p)

    trimmed = Tree( tree.write(is_leaf_fn=is_leaf_node) )


    if args.o:
        write_tree(trimmed, args.o)
    if args.l:
        write_leaves(trimmed, args.l)
