# BiotaPhy analyses

[![Build Status](https://travis-ci.org/biotaphy/analyses.svg?branch=master)](https://travis-ci.org/biotaphy/analyses) [![Coverage Status](https://coveralls.io/repos/github/biotaphy/analyses/badge.svg?branch=master)](https://coveralls.io/github/biotaphy/analyses?branch=master)

Generate ancestral state and distribution reconstructions.

The ancestral distribution code uses a novel approach developed by S.A. Smith and B. O'Meara. Given a set of histograms for species, representing occupancy of environmental space in terms of common bins (i.e., a PNO or predicted niche occupancy profile), this approach reconstructs ancestral histograms of occupancy of climate space. 

This approach is different from those used previously, based on either sampling statistically from present day environmental space or summary statistics (mean, median, maximum, 95th percentile, etc.). Instead of sampling environmental space, probabilities of climate occupancy per bin are explicitly reconstructed. A key advantage of this approach is the ability to reconstruct multimodal ancestral distributions, whereas sampling-based approaches tend to result in normally distributed ancestral reconstructions regardless of extant species distributions.

## Install
```
$ python setup.py install
```

## Usage

### ancestral_distribution.py
```
usage: ancestral_distribution.py [-h] [-l ANNOTATE_LABELS] [-p PLOT_DIRECTORY]
                                 [-c OUT_CSV_FILENAME]
                                 in_tree_filename {newick,nexml,nexus}
                                 data_filename {csv,json,phylip,table}
                                 out_tree_filename {newick,nexml,nexus}

Generates ancestral distribution estimations based on the environmental
distributions at the tips of the tree

positional arguments:
  in_tree_filename      Path to the tree file
  {newick,nexml,nexus}  The format of the tree
  data_filename         Path to file with character state data
  {csv,json,phylip,table}
                        The format of the character data
  out_tree_filename     Path to write the resulting annotated tree
  {newick,nexml,nexus}  The format to use when writing the tree

optional arguments:
  -h, --help            show this help message and exit
  -l ANNOTATE_LABELS, --annotate_labels ANNOTATE_LABELS
                        If provided, annotate the tree labels with this data
                        column
  -p PLOT_DIRECTORY, --plot_directory PLOT_DIRECTORY
                        If provided, write distribution plots to this
                        directory
  -c OUT_CSV_FILENAME, --out_csv_filename OUT_CSV_FILENAME
                        If provided, write the output character matrix CSV to
                        this file location
```

Notes:
  * Use the `-l` option with either a column name or column index to use the reconstructed values for 
    that column as the labels of your output tree.
  * The `-p` option will tell the tool to write out plots for the distributions
  * The `-c` option will write out the reconstruction matrix as a CSV file for processing elsewhere


## Formats

### CSV Alignment file

CSV alignment files can have a header row, where each entry is a label for one of the variables.

The first column of each data row should be a tree tip label followed by the values corresponding with variables in each column.

Each variable can be for a single value, such as mean value, or it can be for a histogram bin.
```
Tip label, Var A, Var B, Var C, Var D - bin 1, Var D - bin 2, Var D - bin 3, Var D - bin 4, Var D - bin 5
Tip A, 1.0, 1.3, 0.5, 3.1, 3.7, 4.1, 3.3, 2.7
Tip B, 0.3, 0.3, 3.0, 0.3, 0.6, 0.2, 0.1, 0.1
Tip C, 0.9, 0.1, 2.5, 1.8, 0.3, 2.1, 2.4, 2.9
Tip D, 1.4, 2.0, 4.3, 3.2, 2.1, 2.3, 1.9, 1.8
...
```
