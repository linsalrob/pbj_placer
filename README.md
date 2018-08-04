# PB_JPLACER

Rewrite jplacer files produced by pplacer and PhyloSift to make trees that are compatible with ITOL and other tree 
viewing software.

## Why should you use this code?


[pplacer](http://matsen.github.io/pplacer/) and related software like [PhyloSift](https://phylosift.wordpress.com/) don't create trees. They map reads
onto trees and output a probabilistic description of where a read should appear on a tree. For more information why, take a look at the 
[pplacer paper](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0031009) that describes the format. The problem with this format
is that it is not compatible with a lot of common phylogenetic software. We have written pbj_placer to reanalyze the tree, add some new nodal information,
modify the format, and to create additional files that can be used to decorate those phylogenetic trees. 
This allows complex data sets to be visualized in, for example, [ITOL](https://itol.embl.de).

## What does pbj_placer do?

We parse out the sequences that should be placed onto the tree, and write them and the tree to separate files.
We then read a user-provided config file that explains what each sample is. By matching reads -> fastq files -> samples
we generate input files that you put into ITOL to make beautifully decorated trees for your publication.

## How do I cite pbj_placer?

First, make sure you cite [phylosift](https://www.ncbi.nlm.nih.gov/pubmed/24482762), [pplacer](https://www.ncbi.nlm.nih.gov/pubmed/21034504), and any other software you use. 
If you really want to cite us, you can cite pbj_placer as *Doane, MP and Edwards, RA. PBJ_PLACER: reformatting probabilistic trees. DOI: 10.5281/zenodo.1312078* 
[![DOI](https://zenodo.org/badge/140918593.svg)](https://zenodo.org/badge/latestdoi/140918593)

## Who do I thank for pbj_placer?

This is experimental code written by Rob Edwards and Mike Doane. You should thank Rob and criticize Mike!

### What we need

Before we begin, we need a few things ready to go:

- a directory of fastq files for your sequences that were used as the input to PhyloSift.
- a classification file that has the fastq file name and then an arbitrary number of classifications. Each classification should be tab-separated. See below for an example.
- the jplacer output file from [phylosift](https://github.com/gjospin/PhyloSift)
- you will need to install the [NCBI Taxonomy SQLite3 database](https://github.com/linsalrob/EdwardsLab/tree/master/taxon)
- download pbj_placer by cloning this repository. `git clone https://github.com/linsalrob/pbj_placer.git` should create a new directory with the code for you.

### What we will produce

We will output several files that you can import into [itol](https://itol.embl.de)

- a newick file with just the reference tree used by phylosift. Note that the tree is rerooted
so that the root is between Archaea and Bacteria or Archaea and Eukarya (depending exactly on the topology of your tree
and where the shared nodes are)
- multibar files:
  - We create a directory with multibars that can either be percentage plots or raw counts.
  - There is one file per sample in the directory, and you can specify the taxonomic levels at which to apply the data.



## Step one, separate the metagenomes and the tree

In this step, we read the jplacer file and separate out the metagenomes from the tree, rename some of the internal nodes in the tree
based on the classification of the leaves associated with the nodes.

Use the command:

```
python3 rename_tree_leaves.py -j sharks_stingray.jplace -o sharks_stingray.nwk -m sharks_stingray.placements
```

to parse the jplace file and create (a) the tree for itol (sharks_stingray.nwk), and (b) a list of all the metagenome reads
and the positions those mapto the tree (sharks_stingray.placements) that we will use in subsequent commands.

This step requires access to the [SQLite3 taxnomy database](https://github.com/linsalrob/EdwardsLab/tree/master/taxon)
that is an interface to NCBI taxonomy. We use that database to figure out our taxonomic level.


## Step two, create our classifications

We need a directory with the fastq files, and then a list of classifications associated with those files.
The list can be of arbitrary length but must be separated with tabs. Missing values should be empty.

You should probably organise this from the highest to lowest classifications (i.e. the left most 
column should have the least number of unique entries). For example, our classification looks like this:

| fastq file | class 1 | class 2 |
| --- | --- | --- | 
| ts.fastq | sharks | thresher shark | 
| ws.fastq | sharks | whale shark |
| bl.fastq | fish | blennie |
| fl.fastq | fish | flounder |

This step also requires access to the [SQLite3 taxnomy database](https://github.com/linsalrob/EdwardsLab/tree/master/taxon)

We generate a file that has all the leaves found in the tree, including the reference sequences, and whether they are 
bacteria, archaea, eukarya, or from your metagenomes. We also append all the classification information you provide.

```
python3 fastq2ids.py -l sharks_stingray.leaves -p -c ../fastq_classification.tsv -d ../fastq -o sharks_stingray.leaves.labels
```

This creates a new file that has several columns (depending on exactly how many metadata classes you provide in your
fastq classification file):

- The original name of the fastq read
- The modified name as it appears in the tree (some special characters like : and . are automatically replaced)
- Whether the read is from a metagenome. This is assumed if it is in a fastq file.
- Which fastq file the read is found in
- The classification(s) provided in the classification file for that fastq file.

Basically, we're mapping leaves in the tree to reads, and reads to fastq files, and fastq files to their taxonomy.


## Step three, count the metagenomes at different levels and create multibars

We create a directory of output files that you can upload to ITOL. Each of the output files can be dropped onto the tree
visualization to decorate the tree. There are a couple of labels that are used in ITOL, e.g. for the legend, and 
we also create the multibar based on the phylogenetic annoations in the tree and the column in the 
labels file that you are keen to display.

The multibar counts the occurrences of metagenomes in different samples, and either uses raw counts or normalized
counts for the display.

For example, to plot the occurrence of different samples at the class level on the tree, and to separate sharks/stingray/fish
we use 

```angular2html
python3 ~redwards/GitHubs/pbj_placer/create_multibar.py -f  sharks_sting_fish.leaves.labels -m sharks_sting_fish.placements -t sharks_sting_fish.nwk -n 4 -x class -d multibar.class.sharkfish -o class.sharkfish.counts.tsv
```

If we change the column in the labels file to read (the -n parameter) we can switch from just looking at sharks, stingray,
and fish to the different species of fish. 

````angular2html
do python3 ~redwards/GitHubs/pbj_placer/create_multibar.py -f  sharks_sting_fish.leaves.labels -m sharks_sting_fish.placements -t sharks_sting_fish.nwk -n 5 -x class -d  multibar.class.counts -o class.counts.tsv
````

Note that this command also creates a file (the .tsv files) that summarizes the data if you want to bring it into
a stats program for other analysis.

