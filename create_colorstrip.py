"""
Create a color strip file from a tsv file that has the sequence ids.
"""

import os
import sys
import argparse

def read_labels(lf, col, verbose=False):
    """
    Read the labels file and return a dict with tree labels and values
    :param lf: labels file
    :param col: the column to use
    :param verbose: extra output
    :return: a dict
    """

    ret = {}
    with open(lf, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            if len(p) < col:
                continue
            if not p[col]:
                continue
            ret[p[0]] = p[col]
    return ret

def read_mapping(mapf, verbose=False):
    """
    Read the mapping from metagenomes to nodes in the tree
    :param mapf: the mapping file
    :param verbose: more output
    :return:
    """

    mapping = {}
    with open(mapf, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            if p[1] not in mapping:
                mapping[p[1]] = set()
            mapping[p[1]].add(p[0])
    return mapping

def remap(data, mapping, verbose=False):
    """
    Map from the mapping file to the data file. In this step we figure out where on the tree
    we should place the color strip
    :param data: the data dictionary where the metagenome read id is the key and the label is the value
    :param mapping: the mapping from metagenome read id to node in the tree
    :param verbose:
    :return:
    """

    ndata = {}

    for r in mapping:
        for t in mapping[r]:
            if t in ndata and ndata[t] != data[r]:
                sys.stderr.write("WARNING: We originally had {} for {} and now we have {}\n".format(
                    ndata[t], t, data[r]
                ))
            ndata[t] = data[r]
    return ndata

def write_output(data, colors, label, lshape, outputfile, verbose):
    """
    Write the colorstrip file
    :param data: the data dict of leaves and valus
    :param colors: the array of colors
    :param label: the label for the color strip
    :param lshape: the label shape
    :param outputfile: the file to write
    :param verbose: more output
    :return:
    """

    vals = list(set(data.values()))
    if len(vals) > len(colors):
        sys.stderr.write("WARNING: NOT ENOUGH COLORS! We have {} values and {} colors\n".format(len(vals), len(colors)))
        sys.exit(-1)

    valcols = {v:colors[vals.index(v)] for v in vals}

    with open(outputfile, 'w') as out:
        out.write("DATASET_COLORSTRIP\n")
        out.write("SEPARATOR COMMA\n")
        out.write(f"DATASET_LABEL,{label}\n")
        out.write("COLOR,#ff0000\n")
        out.write(f"LEGEND_TITLE,{label}\n")
        out.write("LEGEND_COLORS,{}\n".format(",".join(valcols.values())))
        out.write("LEGEND_SHAPES,{}\n".format(",".join([lshape for v in valcols.values()])))
        out.write("LEGEND_LABELS,{}\n".format(",".join(vals)))
        out.write("STRIP_WIDTH,25\n")
        out.write("COLOR_BRANCHES,1\n")
        out.write("DATA\n")
        for d in data:
            out.write("{},{},{}\n".format(d, valcols[data[d]], data[d]))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Create a color strip file")
    parser.add_argument('-f', help='The labeled leaves file from fastq2ids.py', required=True)
    parser.add_argument('-m', help='The mapping file from rename_tree_leaves.py', required=True)
    parser.add_argument('-n', help='Column in the labeled leaves file to use. 0 indexed', required=True, type=int)
    parser.add_argument('-l', help='color strip legend (e.g. Kingdom, Fish, Species', required=True)
    parser.add_argument('-o', help='Output file', required=True)
    parser.add_argument('-s', help='Legend shape (a number). Default = 1', default="1", type=str)
    parser.add_argument('-c', help='Colors to use. These will be prepended to our default list', action='append')
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf', '#999999']
    if args.c:
        colors = args.c + colors

    data = read_labels(args.f, args.n, args.v)
    mapping = read_mapping(args.m, args.v)
    mapdata = remap(data, mapping, args.v)
    write_output(mapdata, colors, args.l, args.s, args.o, args.v)