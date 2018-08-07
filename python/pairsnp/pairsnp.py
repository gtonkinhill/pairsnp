import sys, os
import numpy as np
import argparse
from scipy import sparse

INITIALISATION_LENGTH = 1e5

# This function was taken from https://stackoverflow.com/questions/7654971/parsing-a-fasta-file-using-a-generator-python
def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line[1:], []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


def calculate_consensus(fastafile):
    seq_names = []

    with open(fastafile) as fasta:
        # process first sequence and set up count matrix
        fasta_gen = read_fasta(fasta)
        h,s = fasta_gen.next()
        seq_names.append(h)
        s = np.fromstring(s.lower(), dtype=np.int8)
        align_length = len(s)
        count_matrix = np.zeros((4,align_length))

        count_matrix[0,] = count_matrix[0,] + (s==97)
        count_matrix[1,] = count_matrix[1,] + (s==99)
        count_matrix[2,] = count_matrix[2,] + (s==103)
        count_matrix[3,] = count_matrix[3,] + (s==116)

        for h,s in fasta_gen:
            seq_names.append(h)
            s = np.fromstring(s.lower(), dtype=np.int8)
            count_matrix[0,] = count_matrix[0,] + (s==97)
            count_matrix[1,] = count_matrix[1,] + (s==99)
            count_matrix[2,] = count_matrix[2,] + (s==103)
            count_matrix[3,] = count_matrix[3,] + (s==116)
            

        consensus =  np.array([97,99,103,116])[np.argmax(count_matrix, axis=0)]

    return seq_names, align_length, consensus


def calculate_snp_matrix(align_length, consensus, fastafile, nseqs):
    print "Generating SNP matrix... "
    row = np.empty(INITIALISATION_LENGTH)
    col = np.empty(INITIALISATION_LENGTH, dtype=np.int64)
    val = np.empty(INITIALISATION_LENGTH, dtype=np.int8)
    print "Allocated memory"
    r = 0
    n_snps = 0
    current_length = INITIALISATION_LENGTH
    with open(fastafile) as fasta:
        for h,s in read_fasta(fasta):
            s = np.fromstring(s.lower(), dtype=np.int8)
            snps = consensus!=s
            right = n_snps + np.sum(snps)

            if right >= (current_length/2):
                current_length = current_length + INITIALISATION_LENGTH
                row.resize(current_length)
                col.resize(current_length)
                val.resize(current_length)

            row[n_snps:right] = r
            col[n_snps:right] = np.flatnonzero(snps)
            val[n_snps:right] = s[snps]
            r += 1
            n_snps = right

    row = row[0:right] 
    col = col[0:right]
    val = val[0:right]

    sparse_snps = sparse.csc_matrix((val, (row, col)), shape=(nseqs, align_length))

    return sparse_snps



def main():

    parser = argparse.ArgumentParser(description='Program to calculate pairwise SNP distance and similarity matrices.')

    parser.add_argument('-t', '--type', dest='type', type=str, choices=["sim", "dist"],
                       help='either sim (similarity) or dist (distance).')

    parser.add_argument('-f', '--file', dest='filename', required=True,
                       type=str,
                       help='location of a multiple sequence alignment. Currently only DNA alignments are supported.')

    parser.add_argument('-o', '--out', dest='output', required=True,
                       type=str,
                       help='location of output file.')

    args = parser.parse_args()

    seq_names, align_length, consensus = calculate_consensus(args.filename)

    sparse_matrix = calculate_snp_matrix(align_length, consensus, args.filename, len(seq_names))

    print(np.sum(1*(sparse_matrix[2,].todense()>0)))

    d = (1*(sparse_matrix==97)) * (sparse_matrix.transpose()==97)
    d = d + (1*(sparse_matrix==99) * (sparse_matrix.transpose()==99))
    d = d + (1*(sparse_matrix==103) * (sparse_matrix.transpose()==103))
    d = d + (1*(sparse_matrix==116) * (sparse_matrix.transpose()==116))

    d = d.todense()

    if(args.type!="sim"):
        temp_total = np.zeros((len(seq_names), len(seq_names)))
        seq_sum = (1*(sparse_matrix>0)).sum(1)
        temp_total[:] = seq_sum
        total_differences_shared = (1*(sparse_matrix>0)) * (sparse_matrix.transpose()>0)
        d = temp_total + np.transpose(temp_total) - total_differences_shared.todense() - d
    
    with open(args.output, 'w') as outfile:
        np.savetxt(outfile, d, fmt="%d", delimiter=",", 
            header=",".join(seq_names))

    return



if __name__ == '__main__':
    main()