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
  in_tree_filename      Path to the tree file
  {newick,nexml,nexus}  The format of the tree
  data_filename         Path to file with character state data
  {csv,json,phylip}     The format of the character data
  out_tree_filename     Path for new output tree
  {newick,nexml,nexus}  The format of the tree

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


