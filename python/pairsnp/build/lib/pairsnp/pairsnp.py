import sys, os
import numpy as np
import argparse
from scipy import sparse

INITIALISATION_LENGTH = 100000

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
        s[(s!=97) & (s!=99) & (s!=103) & (s!=116)] = 110
        align_length = len(s)
        count_matrix = np.zeros((5,align_length))

        count_matrix[0,] = count_matrix[0,] + (s==97)
        count_matrix[1,] = count_matrix[1,] + (s==99)
        count_matrix[2,] = count_matrix[2,] + (s==103)
        count_matrix[3,] = count_matrix[3,] + (s==116)
        count_matrix[4,] = count_matrix[4,] + (s==110)

        for h,s in fasta_gen:
            seq_names.append(h)
            s = np.fromstring(s.lower(), dtype=np.int8)
            s[(s!=97) & (s!=99) & (s!=103) & (s!=116)] = 110
            count_matrix[0,] = count_matrix[0,] + (s==97)
            count_matrix[1,] = count_matrix[1,] + (s==99)
            count_matrix[2,] = count_matrix[2,] + (s==103)
            count_matrix[3,] = count_matrix[3,] + (s==116)
            count_matrix[4,] = count_matrix[4,] + (s==110)
            

        consensus =  np.array([97,99,103,116,110])[np.argmax(count_matrix, axis=0)]

    return seq_names, align_length, consensus


def calculate_snp_matrix(fastafile):
    seq_names, align_length, consensus = calculate_consensus(fastafile)
    nseqs = len(seq_names)

    row = np.empty(INITIALISATION_LENGTH)
    col = np.empty(INITIALISATION_LENGTH, dtype=np.int64)
    val = np.empty(INITIALISATION_LENGTH, dtype=np.int8)

    r = 0
    n_snps = 0
    current_length = INITIALISATION_LENGTH
    with open(fastafile) as fasta:
        for h,s in read_fasta(fasta):
            s = np.fromstring(s.lower(), dtype=np.int8)
            s[(s!=97) & (s!=99) & (s!=103) & (s!=116)] = 110
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

    return sparse_snps, consensus, seq_names

def calculate_distance_matrix(sparse_matrix, consensus, type, inc_n):

    n_seqs = sparse_matrix.shape[0]

    d = (1*(sparse_matrix==97)) * (sparse_matrix.transpose()==97)
    d = d + (1*(sparse_matrix==99) * (sparse_matrix.transpose()==99))
    d = d + (1*(sparse_matrix==103) * (sparse_matrix.transpose()==103))
    d = d + (1*(sparse_matrix==116) * (sparse_matrix.transpose()==116))

    d = d.todense()

    if(type=="dist"):
        temp_sparse_n = sparse_matrix==110
        n_comp = (1*temp_sparse_n * (temp_sparse_n.transpose())).todense()
        d = d + n_comp
        temp_total = np.zeros((n_seqs, n_seqs))
        seq_sum = (1*(sparse_matrix>0)).sum(1)
        temp_total[:] = seq_sum
        total_differences_shared = (1*(sparse_matrix>0)) * (sparse_matrix.transpose()>0)

        if inc_n:
            d = temp_total + np.transpose(temp_total) - total_differences_shared.todense() - d
        else:
            n_total = np.zeros((n_seqs, n_seqs))
            n_sum = (1*temp_sparse_n).sum(1)
            n_total[:] = n_sum
            
            if(sum(consensus==110)<=0):
                tot_cons_snps_N = cons_snps_N = 0
            else:
                matrix_n_cols = 1*(sparse_matrix>0)[:,consensus==110]
                cons_snps_N = matrix_n_cols * np.transpose(matrix_n_cols)

                tot_cons_snps_N = np.zeros((n_seqs, n_seqs))
                tot_cons_snps_N_sum = matrix_n_cols.sum(1)
                tot_cons_snps_N[:] = tot_cons_snps_N_sum
            
            diff_n = n_total + np.transpose(n_total) - 2*n_comp + tot_cons_snps_N + np.transpose(tot_cons_snps_N) - 2*cons_snps_N
            d = temp_total + np.transpose(temp_total) - total_differences_shared.todense() - d - diff_n

    return d

def main():

    parser = argparse.ArgumentParser(description='Program to calculate pairwise SNP distance and similarity matrices.')

    parser.add_argument('-t', '--type', dest='type', type=str, choices=["sim", "dist"], default="dist",
                       help='either sim (similarity) or dist (distance) (default).')

    parser.add_argument('-n', '--inc_n', dest='inc_n', action='store_true',
                       help='flag to indicate differences to gaps should be counted.')

    parser.add_argument('-f', '--file', dest='filename', required=True,
                       type=str,
                       help='location of a multiple sequence alignment. Currently only DNA alignments are supported.')

    parser.add_argument('-o', '--out', dest='output', required=True,
                       type=str,
                       help='location of output file.')

    args = parser.parse_args()

    sparse_matrix, consensus, seq_names = calculate_snp_matrix(args.filename)
 
    d = calculate_distance_matrix(sparse_matrix, consensus, args.type, args.inc_n)

    with open(args.output, 'w') as outfile:
        np.savetxt(outfile, d, fmt="%d", delimiter=",", 
            header=",".join(seq_names))

    return



if __name__ == '__main__':
    main()