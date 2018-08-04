"""
Rewrite the leaves in a jplacer file with taxonomic level (One of 'superkingdom', 'phylum', 'class', 'order', 'family',
'genus', 'species', 'subspecies') and write that tree to a newick file.
"""

import os
import sys
import argparse
import json
import re
from ete3 import Tree
from ete3.parser.newick import NewickError


from taxon import get_taxonomy_db, get_taxonomy


def load_jplacer(jpf):
    """
    load the jplacer file and return the tree
    :param jpf: The jplacer file
    :return: the data structure of the tree
    """

    with open(jpf, 'r') as f:
        data = json.load(f)

    return data

def parse_jplacer_tree(data):
    """
    Extract the tree from the jplacer data structure and make it into an ete3 object
    :param data: the jplacer data structure
    :return:
    """

    try:
        tree = Tree(data['tree'], quoted_node_names=True, format=1)
    except NewickError as n:
        tt = re.sub(r'(\:[\d\.]+){\d+}', r'\1', data['tree'])
        tt = re.sub(r'{\d+};$', ';', tt)
        tree = Tree(tt, quoted_node_names=True, format=1)

    return tree

def rename_nodes_ncbi(tree, verbose=False):
    """
    Rename the nodes based on everything below me, but also give each node a unique branch number.
    The format of this number is _b\d+

    :param tree: the tree to rename
    :param verbose: more output
    :return: the renamed tree
    """

    # connect to the SQL dataabase
    c = get_taxonomy_db()

    wanted_levels = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'subspecies']
    wanted_levels.reverse()  # too lazy to write in reverse :)
    taxonomy = {}
    # first get all the leaves and their parents. This is just to speed things up ... maybe
    for l in tree.get_leaves():
        m = re.search('\[(\d+)\]', l.name)
        if not m:
            if verbose:
                sys.stderr.write("No taxid in {}\n".format(l.name))
            continue
        tid = m.groups()[0]
        taxonomy[l.name] = {}
        t, n = get_taxonomy(tid, c)
        if not t:
            continue
        while t.parent != 1 and t.taxid != 1:
            if t.rank in wanted_levels:
                taxonomy[l.name][t.rank] = n.scientific_name
            t, n = get_taxonomy(t.parent, c)

    # now traverse every node that is not a leaf and see if we can some up with a
    # unique name for the node!
    if verbose:
        sys.stderr.write("Traversing the tree to rename the nodes\n")
    branchnum = 0
    for n in tree.traverse("postorder"):
        if n.is_leaf():
            continue
        ## if both our children have the same name, we acquire that name and reset their names
        ## otherwise we figure out what our name should be based on the last common ancestor
        ## of the leaves.

        children = n.get_children()
        names = set([re.sub('b_\d+', '', x.name) for x in children])
        if len(names) == 1:
            n.name = "{} b_{}".format(names.pop(), branchnum)
            if verbose:
                sys.stderr.write("Reset name to {} because both children are the same\n".format(n.name))
            for c in children:
                oldname = c.name
                c.name = re.sub('r_\w+', '', c.name)
                if verbose:
                    sys.stderr.write("\tAs both children the same set {} to {}\n".format(oldname, c.name))
        else:
            ## We have to figure out what our unique name should be
            taxs = {w: set() for w in wanted_levels}
            for l in n.get_leaves():
                if l.name not in taxonomy:
                    continue
                for w in wanted_levels:
                    if w in taxonomy[l.name]:
                        taxs[w].add(taxonomy[l.name][w])
            # which is the LOWEST level with a single taxonomy
            for w in wanted_levels:
                if len(taxs[w]) == 1:
                    newname = "{} r_{} b_{}".format(taxs[w].pop(), w, branchnum)
                    if verbose:
                        sys.stderr.write("Changing name from: {} to {}\n".format(n.name, newname))
                    n.name = newname
                    break
        branchnum += 1
    return tree

def reroot_tree(tree, verbose=False):
    """
    Reroot the tree between bacteria and archaea.

    This will only work after renaming the leaves on the tree.

    :param tree: the tree
    """

    didreroot = False

    if verbose:
        sys.stderr.write("rerooting the tree\n")
    for n in tree.traverse("preorder"):
        childs = n.get_children()
        if verbose:
            cname = ""
            for c in childs:
                cname += "| {} |".format(c.name)
            sys.stderr.write("{}\t{}\t{}\n".format(len(childs), n.name, cname))
        if len(childs) == 2:
            if ("Archaea r_superkingdom" in childs[0].name and "Eukaryota r_superkingdom" in childs[1].name) or (
                    "Archaea r_superkingdom" in childs[1].name and "Eukaryota r_superkingdom" in childs[0].name):
                tree.set_outgroup(n)
                if verbose:
                    sys.stderr.write("Rerooted on {}\n".format(n.name))
                didreroot = True
                break
            if "Bacteria r_superkingdom" in childs[0].name and "Archaea r_superkingdom" in childs[1].name:
                tree.set_outgroup(childs[0])
                if verbose:
                    sys.stderr.write("Rerooted on {}\n".format(childs[0].name))
                didreroot = True
                break
            if "Bacteria r_superkingdom" in childs[1].name and "Archaea r_superkingdom" in childs[0].name:
                tree.set_outgroup(childs[1])
                if verbose:
                    sys.stderr.write("Rerooted on {}\n".format(childs[1].name))
                didreroot = True
                break

    if not didreroot:
        for n in tree.traverse("preorder"):
            if "Bacteria r_superkingdom" in n.name:
                tree.set_outgroup(n)
                if verbose:
                    sys.stderr.write("Rerooted on {} because it is bacteria\n".format(n.name))
                break

    return tree

def clean_newick_id(name):
    """
    Return a version of name suitable for placement in a newick file
    :param name: The name to clean up
    :return: a name with no colons, spaces, etc
    """
    name = name.replace(' ', '_')
    name = name.replace(':', '_')
    name = name.replace('[', '_')
    name = name.replace(']', '_')

    return name

def get_placements(data):
    """
    Get the placements and return a dict with the keys being the edge numbers where to do
    the insertions and the values being a set of nodes to insert at that point.

    TODO: we have not implemented the approach for a single insertion as we don't have an example of that (yet!)

    :param data: the parsed jplacer tree
    :return: a dict of placement edge_numbers and sets of ids to add
    """

    # first make sure the tree fields are in the correct order!
    posn = data['fields'].index('edge_num')
    placements = {}

    for pl in data['placements']:
        addhere = set()
        if 'n' in pl:
            sys.stderr.write("Crap, not sure what to do because I've never seen an example. You should be able to figure out from what I did with nm\n")
            sys.exit(-1)
        if 'nm' in pl:
            for i in pl['nm']:
                addhere.add(clean_newick_id(i[0]))
        for p in pl['p']:
            edge_num = p[posn]
            if edge_num not in placements:
                placements[edge_num] = set()
            placements[edge_num].update(addhere)

    return placements

def write_placement_tuples(pl, tree, tpoutfile, verbose=False):
    """
    Write a file with the tuples of new edge node (from metagenome) and existing node where it would be inserted
    :param pl: The placements from get_placements
    :param tree: The tree (either with or without rewriting)
    :param tpoutfile: The file to write the tuples to
    :param verbose: more information
    :return:
    """

    with open(tpoutfile, 'w') as out:
        for t in tree.traverse("postorder"):
            m = re.search('{(\d+)}', t.name)
            if not m:
                continue
            thisid = int(m.groups()[0])

            if thisid in pl:
                for p in pl[thisid]:
                    out.write("{}\t{}\n".format(clean_newick_id(t.name), p))

def write_tree(tree, outputf):
    """
    Write the tree to a file.
    :param tree: The tree to write
    :param outputf: The output filename
    :return:
    """

    tree.write(outfile=outputf, format=1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse a jplacer file and rewrite the leaves with the taxonomy')
    parser.add_argument('-j', help='jplacer file', required=True)
    parser.add_argument('-o', help='output file to write the tree to', required=True)
    parser.add_argument('-m', help='write a mapping file from the mapped nodes to nodes on the tree.' +
                                   ' Essentially converts the jplacer placements into tsv for later processing')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    data = load_jplacer(args.j)
    tree = parse_jplacer_tree(data)
    tree = rename_nodes_ncbi(tree, args.v)
    tree = reroot_tree(tree, args.v)
    write_tree(tree, args.o)

    if args.m:
        pl = get_placements(data)
        write_placement_tuples(pl, tree, args.m, args.v)


