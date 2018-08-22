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

The majority of time is spend doing sparse matrix multiplications so linking to a parallelised library for this is likely to improve preformance further.

At the moment there you may need to run `touch ./cpp/*` before compiling to avoid some issues with time stamps.

```
cd ./pairsnp/cpp/
./configure
make
make install
```

