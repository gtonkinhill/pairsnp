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

At the moment the python library is not set up to work as a package but I'm planning on doing this at some point along with transferring it to python3.

It relies on scipy and numpy

```
pip install scipy
pip install numpy
```

After these are installed the script can simply be called by

```
python ../pairsnp.py -h
```

### c++

The c++ code relies on a recent version of Armadillo and can be built by running

```
cd ./pairsnp/cpp/
./configure
make
make install
```

