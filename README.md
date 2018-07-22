# PB_JPLACER

Rewrite jplacer files produced by pplacer and PhyloSift to make trees that are compatible with ITOL and other tree viewing software.

## Why should you use this code?


[pplacer](http://matsen.github.io/pplacer/) and related software like [PhyloSift](https://phylosift.wordpress.com/) don't create trees. They map reads
onto trees and output a probabilistic description of where a read should appear on a tree. For more information why, take a look at the 
[pplacer paper](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0031009) that describes the format. The problem with this format
is that it is not compatible with a lot of common phylogenetic software. We have written pbj_placer to rewrite the format, and to create additional
files that can be used to decorate those phylogenetic trees. This allows complex data sets to be visualized in, for example, [ITOL](https://itol.embl.de).

## What does pbj_placer do?

We parse out the sequences that should be placed onto the tree, and write them and the tree to separate files.
We then read a user-provided config file that explains what each sample is. By matching reads -> fastq files -> samples,
we generate input files that you put into ITOL.

## How do I cite pbj_placer?

First, make sure you cite [phylosift](https://www.ncbi.nlm.nih.gov/pubmed/24482762), [pplacer](https://www.ncbi.nlm.nih.gov/pubmed/21034504), and any other software you use. 
If you really want to cite us, you can cite pbj_placer as *Doane, MP and Edwards, RA. PBJ_PLACER: reformatting probabilistic trees. DOI: 10.5281/zenodo.1312078* 
[![DOI](https://zenodo.org/badge/140918593.svg)](https://zenodo.org/badge/latestdoi/140918593)

## Who do I thank for pbj_placer?

This is experimental code written by Rob Edwards and Mike Doane. You should thank Rob and criticie Mike!

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
- colorstrip files: These files make the strips around the outside of the images.
  - Kingdom: Whether the sample is from Archaea, Metagenome, Eukaryota, Bacteria, Unknown
  - Species: this file has all the species in the samples.
- multibar files:
  - We create a directory with multibars that can either be percentage plots or raw counts.
  - There is one file per sample in the directory, and you can specify the taxonomic levels at which to apply the data.



## Step one, separate the metagenomes and the tree

In this step, we read the jplacer file and separate out the metagenomes from the tree. 

Use the command:

```
python3 rename_tree_leaves.py -j sharks_stingray.jplace -o sharks_stingray.nwk -m sharks_stingray.placements
```

to parse the jplace file and create (a) the tree for itol (sharks_stingray.nwk), and (b) a list of all the metagenome reads
and the positions those mapto the tree (sharks_stingray.placements) that we will use in subsequent commands.

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

## Step three, create color strips for different taxonomic levels

We can create colorstrips for e.g. Bacteria, Archaea, Metagenomes or for Sharks and Fish based on the data in our 
leaves.labels file that we have just created. e.g. to make a file for the Kingdom, use:

```
python3 create_colorstrip.py -f sharks_stingray.leaves.labels -n 2 -l Kingdom -s 1 -o sharks_stingray.kingdom.colorstrip
```

and to create a similar file for shark or fish, use:

```
python3 create_colorstrip.py -f sharks_stingray.leaves.labels -n 4 -l "Fish/Shark" -s 2 -o sharks_stingray.fish_shark.colorstrip
```

and to create one for each of the fish/shark species we use:

```
python3 create_colorstrip.py -f sharks_stingray.leaves.labels -n 5 -l Species -s 3 -o sharks_stingray.species.colorstrip
```

or, just do them all at once and output everything to a directory:
```angular2html
mkdir colorstrip/
python3 create_colorstrip.py -f sharks_stingray.leaves.labels -n 2 -l Kingdom -s 1 -o colorstrip/sharks_stingray.kingdom.colorstrip
python3 create_colorstrip.py -f sharks_stingray.leaves.labels -n 4 -l "Fish/Shark" -s 2 -o colorstrip/sharks_stingray.fish_shark.colorstrip
python3 create_colorstrip.py -f sharks_stingray.leaves.labels -n 5 -l Species -s 3 -o colorstrip/sharks_stingray.species.colorstrip
```


## Step four, count the metagenomes at different levels and create multibars

This creates a directory with a file for each of the bars we need to plot. The first command normalizes the data to percentages. The second command gives you raw counts.

```
python3 create_multibar.py -f sharks_stingray.leaves.labels -t sharks_stingray.nwk -d multibar -n 5 -l Species -x family -s 3 -p
python3 create_multibar.py -f sharks_stingray.leaves.labels -t sharks_stingray.nwk -d multibarnp -n 5 -l Species -x family -s 3
```



