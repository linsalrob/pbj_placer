"""
Just print the children for a node. This is to help explore the tree.
"""

import os
import sys
import argparse
from ete3 import Tree
from ete3.parser.newick import NewickError
import re

def read_tree(treefile):
    """
    Read the tree file and return the tree
    :param treefile: The tree file to read
    :return: The ete3 tree object
    """

    return Tree(treefile, quoted_node_names=True, format=1)

def print_children(tree, node):
    """
    Print the children of 'node'
    :param tree: The tree
    :param node: The name node to find the children of
    :return:
    """

    for n in tree.traverse('preorder'):
        if n.name == node:
            children = n.get_children()
            print("Children for {}:\n\t{}".format(n.name, "\n\t".join([x.name for x in children])))

def regexp_children(tree, node):
    """
    Find the node to print the children of 'node' by regexp
    :param tree: The tree
    :param node: The name node to find the children of
    :return:
    """

    for n in tree.traverse('preorder'):
        if re.search(node, n.name):
            children = n.get_children()
            print("Children for {}:\n\t{}".format(n.name, "\n\t".join([x.name for x in children])))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Print the children of a tree")
    parser.add_argument('-t', help='tree file', required=True)
    parser.add_argument('-n', help='full name of the node')
    parser.add_argument('-r', help='partial name of the node to be matched with regexp')
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    tree = read_tree(args.t)
    if args.n:
        print_children(tree, args.n)
    elif args.r:
        regexp_children(tree, args.r)
    else:
        sys.stderr.write("Sorry, one of -r or -n must be specified\n")
        sys.exit(-1)
