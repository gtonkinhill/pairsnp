# pairsnp

A set of scripts for very quickly obtaining pairwise SNP distance matrices from multiple sequence alignments using sparse matrix libraries or bitset operations to improve performance.

For larger alignments such as the Maela pneumococcal data set (3e5 x 3e3) the c++ version is approximately an order of magnitude faster than approaches based on pairwise comparison of every site such as [snp-dists](https://github.com/tseemann/snp-dists) from which the skeleton code for the c++ version was taken.

In order to be most useful implementations in c++, python (now implemented in Tracs) and R are available.

| Implementation        | Travis           |
| ------------- |:-------------:|
| [c++](https://github.com/gtonkinhill/pairsnp-cpp) | [![pairsnp-CI](https://github.com/gtonkinhill/pairsnp-cpp/actions/workflows/pairsnp_test.yml/badge.svg)](https://github.com/gtonkinhill/pairsnp-cpp/actions/workflows/pairsnp_test.yml)   |
| [python](https://github.com/gtonkinhill/tracs)      | [![tracs-CI](https://github.com/gtonkinhill/tracs/actions/workflows/tracs_test.yml/badge.svg)](https://github.com/gtonkinhill/tracs/actions/workflows/tracs_test.yml)      |
| [R](https://github.com/gtonkinhill/pairsnp-r)     | No longer supported |


## Installation

### c++

The c++ version can be installed manually, by downloading the binaries in this repository, or with conda as

```
conda install -c bioconda pairsnp
```

The c++ code relies on a recent version of Armadillo (currently tested on v8.6) and after downloading the [repository](https://github.com/gtonkinhill/pairsnp-cpp) can be built by running

```
cd ./pairsnp-cpp/
./configure
make
make install
```

### python

The python version is now included in the Tracs pipeline and can be installed using conda or `pip`.

```
conda install bioconda::tracs
```

```
pip3 install git+https://github.com/gtonkinhill/tracs
```

### R

The R version is not longer supported but should be able to be installed using devtools or downloaded from its [repository](https://github.com/gtonkinhill/pairsnp-r)

```
#install.packages("devtools")
devtools::install_github("gtonkinhill/pairsnp-r")
```


## Quick Start

### c++

The c++ version can be run from the command line as

```
pairsnp -c msa.fasta > output.csv
```

additional options include

```
SYNOPSIS
  Pairwise SNP distance matrices using fast matrix algerbra libraries
USAGE
  pairsnp [options] alignment.fasta[.gz] > matrix.csv
OPTIONS
  -h	Show this help
  -v	Print version and exit
  -s	Output in sparse matrix form (i,j,distance).
  -d	Distance threshold for sparse output. Only distances <= d will be returned.
  -k	Will on return the k nearest neighbours for each sample in sparse output.
  -c	Output CSV instead of TSV
  -n	Count comparisons with Ns (off by default)
  -t	Number of threads to use (default=1)
  -b	Blank top left corner cell instead of 'pairsnp 0.1.0'
```

### Python

See the [Tracs documentation](https://gthlab.au/tracs/#/)

### R

```
library(pairsnp)
fasta.file.name <- system.file("extdata", "seqs.fa", package = "pairsnp")
sparse.data <- import_fasta_sparse(fasta.file.name)
d <- snp_dist(sparse.data)
```




