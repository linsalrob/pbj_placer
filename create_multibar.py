"""
Create a multibar file from our labeled file and our tree. We will do this at different taxonomic levels
and different label levels.

For each label level we have, we'll create a single file.

"""

import os
import sys
import argparse
from ete3 import Tree


# TODO:
# There appears to be great redundancy in the file sharks_stingray.leaves.labels
# so need to figure out where that comes from as the counts in read_labels were wrong

def read_labels(lf, col, verbose=False):
    """
    Read the labels file and return a dict with tree labels and values
    :param lf: labels file
    :param col: the column to use
    :param verbose: extra output
    :return: a dict of the leaves and their labels and a dict of the labels and their counts
    """

    ret = {}
    mreads = {}
    with open(lf, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            if len(p) <= col:
                continue
            if not p[col]:
                continue

            ret[p[0]] = p[col]
            if p[col] not in mreads:
                mreads[p[col]] = set()
            mreads[p[col]].add(p[0])

    counts = {x:len(mreads[x]) for x in mreads}

    if verbose:
        sys.stderr.write("After read_labels: data has {} keys and counts has {} keys\n".format(len(ret.keys()), len(counts.keys())))

    return ret, counts

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

    s=set()
    for m in mapping:
        s.update(mapping[m])
    if verbose:
        sys.stderr.write("After read_mapping: mapping has {} keys and {} values\n".format(len(mapping.keys()), len(s)))

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

    if verbose:
        sys.stderr.write("There are {} keys in mapping\n".format(len(mapping.keys())))
    s = set()
    for k in mapping:
        s.update(mapping[k])


    if verbose:
        sys.stderr.write("There are {} vals in mapping\n".format(len(s)))
    for mgid in mapping:
        for posn_in_tree in mapping[mgid]:
            if posn_in_tree not in ndata:
                ndata[posn_in_tree] = {}
            ndata[posn_in_tree][data[mgid]] = ndata[posn_in_tree].get(data[mgid], 0) + 1
    if verbose:
        sys.stderr.write("There are  {} Keys in ndata\n".format(len(ndata.keys())))

    return ndata

def multibar_counts(treefile, data, taxa, proportions, verbose=False):
    """
    Calculate the counts that will be added to the multibar and return a mutlidimensional
    dict of shark type, tree name, and count.
    :param treefile: The tree file in newick format
    :param data: The mapped data we need to read
    :param taxa: The taxonomic level we desire
    :param proportions: whether to use counts or proportions
    :param verbose: more output
    :return: a dict of dicts.
    """

    total = {} ## the total number of times we see a node
    val = {} ## how many times we see the children of this node

    if verbose:
        sys.stderr.write("Reading tree\n")

    tree = Tree(treefile, quoted_node_names=True, format=1)

    for k in counts:
        val[k] = {}
        for n in tree.traverse("preorder"):
            if taxa in n.name:
                for l in n.iter_descendants():
                    if l.name in data and k in data[l.name]:
                        val[k][n.name] = val[k].get(n.name, 0) + data[l.name][k]
                        total[n.name] = total.get(n.name, 0) + data[l.name][k]
                        if verbose:
                            sys.stderr.write("Value for {} and {} is now {}\n".format(k, n.name, val[k][n.name]))
                    elif l.name in data and verbose:
                        sys.stderr.write("Skipped {} as it is a {} and we're a {}\n".format(l.name, data[l.name], k))

    if proportions:
        # how many times did we see each thing:
        for k in val:
            sums = sum(val[k].values())
            for n in val[k]:
                val[k][n] /= sums
    if verbose:
        sys.stderr.write("We found a total of {} metagenomes\n".format(sum(total.values())))


    return val


def write_directory(counts, outputdir, colors, proportions, verbose=False):
    """
    Write a directory with one multibar file per type
    :param counts: the dict of dicts. The first key is the shark type, the second the genus/species
    :param outputdir: the directory to create
    :param colors: the array of colors to choose from
    :param proportions: whether we are using proportions or not
    :param verbose: more output
    :return:
    """

    # assign colors the keys
    allkeys = list(counts.keys())
    if len(allkeys) > len(colors):
        sys.stderr.write("ERROR: Not enough colors. We have {}  keys and {} colors\n".format(len(allkeys), len(colors)))
        sys.exit(-1)
    keycolors = {x: colors[allkeys.index(x)] for x in allkeys}

    # what is our maxvalue
    maxval = 50
    for k in counts:
        m = max(counts[k].values())
        if m > maxval:
            maxval = m


    if not os.path.exists(outputdir):
        try:
            os.mkdir(outputdir)
        except Exception as e:
            sys.stderr.write("Cannot make directory: {}\n".format(outputdir))
            sys.stderr.write("{}\n".format(e))
            sys.exit(-1)

    if verbose:
        sys.stderr.write(f"Creating output files in {outputdir}\n")

    for k in counts:
        fnme = k.replace(' ', '_')
        outputf = os.path.join(outputdir, fnme + ".multibar.txt")

        dsl = "DATASET_LABEL,Count of {} reads\n".format(k)
        if proportions:
            dsl = "DATASET_LABEL,Proportion of {} reads\n".format(k)

        with open(outputf, 'w') as out:
            out.write("DATASET_MULTIBAR\nSEPARATOR COMMA\n")
            out.write("{}\n".format(dsl))
            out.write("FIELD_COLORS,{}\n".format(keycolors[k]))
            out.write("FIELD_LABELS,{}\n".format(k))
            out.write("DATASET_SCALE,0-{}-{}\n".format(k, keycolors[k]))
            out.write("BORDER_WIDTH,5\n")
            out.write("BORDER_COLOR,{}\n".format(keycolors[k]))
            out.write("WIDTH,{}\n".format(maxval))
            out.write("HEIGHT_FACTOR,50\n")
            out.write("SHOW_INTERNAL,1\n")
            out.write("ALIGN_FIELDS,1\n")
            out.write("COLOR,{}\n".format(keycolors[k]))
            out.write("DATA\n")
            for n in counts[k]:
                if counts[k][n] > 0:
                    out.write("{},{}\n".format(n, counts[k][n]))


def write_tsv(counts, taxa, outputfile, verbose=False):
    """
    Write the counts as a tsv file for some stats.
    :param counts: The counts dict that has keys as sharks and then keys as taxa and values as counts
    :param taxa: the taxonomic level we're choosing
    :param outputfile: the file to write
    :param verbose: more output
    :return: nothing
    """

    allkeys = counts.keys()
    allvals = set()
    for k in allkeys:
        allvals.update(set(counts[k].keys()))

    with open(outputfile, 'w') as out:
        out.write("{}\t{}\n".format(taxa, "\t".join(allkeys)))
        for v in allvals:
            out.write(v)
            for k in allkeys:
                out.write("\t{}".format(counts[k].get(v, 0)))
            out.write("\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-f', help='The labeled leaves file from fastq2ids.py', required=True)
    parser.add_argument('-m', help='Mapping file from rename_tree_leaves.py', required=True)
    parser.add_argument('-t', help='Newick tree file', required=True)
    parser.add_argument('-d', help='Output directory where to write the files', required=True)
    parser.add_argument('-n', help='Column in the labeled leaves file to use. 0 indexed', required=True, type=int)
    parser.add_argument('-x', help='taxa to use for the labels', required=True)
    parser.add_argument('-p', help='Display proportion of counts not counts', action='store_true')
    parser.add_argument('-c', help='Colors to use. These will be prepended to our default list', action='append')
    parser.add_argument('-o', help='tsv file to write with all the data')
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf', '#999999']
    if args.c:
        colors = args.c + colors

    data, counts = read_labels(args.f, args.n, args.v)
    mapping = read_mapping(args.m, args.v)
    mapdata = remap(data, mapping, args.v)

    allowed_taxa = ['r_superkingdom', 'r_phylum', 'r_class', 'r_order', 'r_family', 'r_genus', 'r_species', 'r_subspecies']

    taxa = args.x
    if not taxa.startswith('r_'):
        taxa = "r_{}".format(taxa)

    if taxa not in allowed_taxa:
        sys.stderr.write("Sorry: {} is not an allowed taxa. Your choices are\n{}\n".format(taxa, " ".join(allowed_taxa)))
        sys.exit(-1)

    mbcounts = multibar_counts(args.t, mapdata, taxa, args.p, args.v)

    write_directory(mbcounts, args.d, colors, args.p, args.v)

    if args.o:
        write_tsv(mbcounts, args.x, args.o, args.v)