# pairsnp

A set of scripts for very quickly obtaining pairwise SNP distance matrices from multiple sequence alignments using sparse matrix libraries to improve performance.

For larger alignments such as the Maela pneumoccocal dataset (3e5 x 3e3) the c++ version is approximately an order of magnitude faster than approaches based on pairwise comaprison of every site such as [snp-dists](https://github.com/tseemann/snp-dists).

In order to be most useful implementations in R, python and c++ are available.

## Installation

### R
```
#install.packages("devtools")
devtools::install_github("gtonkinhill/pairsnp", subdir = "R/pairsnp")
```

### python

The python version can be installed using `pip` or by downloading the repository and running `setup.py`.

At the moment it is only available in python2 but I'm planning on converting it to python3.

```
python -m pip install pairsnp
```

or alternatively download the repository and run

```
cd ./pairsnp/python/pairsnp
python ./setup.py install
```

### c++

The c++ code relies on a recent version of Armadillo (currently tested on v8.6) and can be built by running

```
cd ./pairsnp/cpp/
./configure
make
make install
```

The majority of time is spend doing sparse matrix multiplications so linking to a parallelised library for this is likely to improve preformance further.

At the moment there you may need to run `touch ./cpp/*` before compiling to avoid some issues with time stamps.

## Quick Start

### R

```
library(pairsnp)
fasta.file.name <- system.file("extdata", "seqs.fa", package = "pairsnp")
sparse.data <- import_fasta_sparse(fasta.file.name)
d <- snp_dist(sparse.data)
```

### python

The python version can be run from the python interpreter as 

```
from pairsnp import calculate_snp_matrix, calculate_distance_matrix

sparse_matrix, consensus, seq_names = calculate_snp_matrix(fasta.file.name)
d = calculate_distance_matrix(sparse_matrix, consensus, "dist", False)
```

alternatively if installed using pip it can be used at the command line as


```
pairsnp -f /path/to/msa.fasta -o /path/to/output.csv
```
additional options include

```
Program to calculate pairwise SNP distance and similarity matrices.

optional arguments:
  -h, --help            show this help message and exit
  -t {sim,dist}, --type {sim,dist}
                        either sim (similarity) or dist (distance) (default).
  -n, --inc_n           flag to indicate differences to gaps should be
                        counted.
  -f FILENAME, --file FILENAME
                        location of a multiple sequence alignment. Currently
                        only DNA alignments are supported.
  -o OUTPUT, --out OUTPUT
                        location of output file.
```

### c++

The c++ version can be run from the command line as

```
pairsnp -c msa.fasta > output.csv
```

additional options include

```
SYNOPSIS
  Pairwise SNP similarity and distance matrices using fast matrix algerbra libraries
USAGE
  pairsnp [options] alignment.fasta[.gz] > matrix.csv
OPTIONS
  -h	Show this help
  -v	Print version and exit
  -s	Find the similarity matrix
  -c	Output CSV instead of TSV
  -n	Count comparisons with Ns (off by default)
  -b	Blank top left corner cell instead of 'pairsnp 0.0.1'
```

