# ancestral_state

[![Build Status](https://travis-ci.org/biotaphy/ancestral_state.svg?branch=master)](https://travis-ci.org/biotaphy/ancestral_state) [![Coverage Status](https://coveralls.io/repos/github/biotaphy/ancestral_state/badge.svg?branch=master)](https://coveralls.io/github/biotaphy/ancestral_state?branch=master)

Generate ancestral state and distribution reconstructions.

## Install
```
$ python setup.py install
```

## Usage

### ancestral_state.py
```
usage: ancestral_state.py [-h] [-ch COLUMN_HEADERS]
                          in_tree_filename {newick,nexml,nexus} data_filename
                          {csv,json,phylip} out_tree_filename
                          {newick,nexml,nexus}

positional arguments:
  in_tree_filename          Path to the tree file
  {newick,nexml,nexus}      The format of the tree
  data_filename             Path to file with character state data
  {csv,json,phylip}         The format of the character data
  out_tree_filename         Path for new output tree
  {newick,nexml,nexus}      The format of the tree
  out_characters_filename   A file location to write the character data as CSV

optional arguments:
  -h, --help            show this help message and exit
  -ch COLUMN_HEADERS, --column_headers COLUMN_HEADERS
                        Path to file containing column headers, one per line
```

### ancestral_distribution.py
```
usage: ancestral_distribution.py [-h] [-o OUTPUT_DIRECTORY]
                                 in_tree_filename {newick,nexml,nexus}
                                 data_filename {csv,json,phylip,table}

Generates ancestral distribution estimations based on the environmental
distributions at the tips of the tree

positional arguments:
  in_tree_filename      Path to the tree file
  {newick,nexml,nexus}  The format of the tree
  data_filename         Path to file with character state data
  {csv,json,phylip,table}
                        The format of the character data

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        If provided, write the rates and plots here
```

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
