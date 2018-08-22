#include <stdio.h>
#include <getopt.h>
#include <iostream>
#include <fcntl.h>
#include <vector>
#include <armadillo>
#include "kseq.h"


#include <time.h>
#include <ctime>

#define VERSION "0.0.1"
#define EXENAME "pairsnp"

using namespace arma;

void show_help(int retcode)
{
  FILE* out = (retcode == EXIT_SUCCESS ? stdout : stderr);
  fprintf(out, "SYNOPSIS\n  Pairwise SNP similarity and distance matrices using fast matrix algerbra libraries\n");
  fprintf(out, "USAGE\n  %s [options] alignment.fasta[.gz] > matrix.csv\n", EXENAME);
  fprintf(out, "OPTIONS\n");
  fprintf(out, "  -h\tShow this help\n");
  fprintf(out, "  -v\tPrint version and exit\n");
  fprintf(out, "  -s\tFind the similarity matrix\n");
  fprintf(out, "  -c\tOutput CSV instead of TSV\n");
  fprintf(out, "  -n\tCount comparisons with Ns (off by default)\n");
  fprintf(out, "  -b\tBlank top left corner cell instead of '%s %s'\n", EXENAME, VERSION);
  exit(retcode);
}

int main(int argc, char *argv[])
{
  // clock_t begin_time = clock();

  // parse command line parameters
  int opt, quiet=0, csv=0, corner=1, allchars=0, keepcase=0, dist=1, count_n=0;
  while ((opt = getopt(argc, argv, "hqsncv")) != -1) {
    switch (opt) {
      case 'h': show_help(EXIT_SUCCESS); break;
      case 'q': quiet=1; break;
      case 'c': csv=1; break;
      case 's': dist=0; break;
      case 'n': count_n=1; break;
      case 'v': printf("%s %s\n", EXENAME, VERSION); exit(EXIT_SUCCESS);
      default : show_help(EXIT_FAILURE);
    }
  }

  char sep = csv ? ',' : '\t';

  // require a filename argument
  if (optind >= argc) {
    show_help(EXIT_FAILURE);
    return 0;
  }
  const char* fasta = argv[optind];

  // say hello
  if (!quiet) fprintf(stderr, "This is %s %s\n", EXENAME, VERSION);

  // open filename via libz
  int fp = open(fasta, O_RDONLY);
  if (! fp) {
    fprintf(stderr, "ERROR: Could not open filename '%s'\n", fasta);
    exit(EXIT_FAILURE);
  }

  
  // cout << "warmup. " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
  // begin_time = clock();

  //Initially run through fasta to get consensus sequence and dimensions of matrix
  int n = 0;
  int l = 0;
  kseq seq;
  FunctorRead r;
  kstream<int, FunctorRead> ks(fp, r);

  // Use first sequence to set size of arrays
  l = ks.read(seq);
  int seq_length = seq.seq.length();
  int allele_counts[5][seq_length];
  memset(allele_counts, 0, 5*seq_length*sizeof(int));

  while((l = ks.read(seq)) >= 0) {
    for(int j=0; j<seq_length; j++){
      if((seq.seq[j]=='a') || (seq.seq[j]=='A')){
        allele_counts[0][j] += 1;
      } else if((seq.seq[j]=='c') || (seq.seq[j]=='C')){
        allele_counts[1][j] += 1;
      } else if((seq.seq[j]=='g') || (seq.seq[j]=='G')){
        allele_counts[2][j] += 1;
      } else if((seq.seq[j]=='t') || (seq.seq[j]=='T')){
        allele_counts[3][j] += 1;
      } else {
        allele_counts[4][j] += 1;
      }
    }
    n++;
  }
  close(fp);

  
  // cout << "initial read.. " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
  // begin_time = clock();

  // Now calculate the consensus sequence
  uvec consensus(seq_length);
  int n_seqs = n+1;
  std::vector<std::string> seq_names;

  int max_allele = -1;

  for(int j=0; j<seq_length; j++){
    for(int i=0; i<5; i++){
      if(allele_counts[i][j]>max_allele){
        max_allele = allele_counts[i][j];
        consensus(j) = i;
      }
    }
    max_allele = -1;
  }

  
  // cout << "consensus calc.. " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
  // begin_time = clock();

  // create matrix to store snp locations
  std::vector<uint> m_i;
  std::vector<uint> m_j;
  std::vector<uint> m_x;
  m_i.reserve(50000);
  m_j.reserve(50000);
  m_x.reserve(50000);

  // open a new stream to the fasta file TODO: I think there's a cleaner way of doing this.
  int fp2 = open(fasta, O_RDONLY);
  kstream<int, FunctorRead> ks2(fp2, r);
  char temp_char;
  int n_snps = 0;
  n=0;
  std::string nt ("AaCcGgTt");

  while((l = ks2.read(seq)) >= 0) {
    // Record sequence names
    auto it = std::find(seq.name.begin(), seq.name.end(), '>');
    if (it != seq.name.end()){
      seq.name.erase(it);
    }

    seq_names.push_back(seq.name);

    if ((m_i.capacity() - 2*n_snps)<100){
      // Reserve memory
      m_i.reserve( 2*m_i.capacity() );
      m_j.reserve( 2*m_j.capacity() );
      m_x.reserve( 2*m_x.capacity() );
    }
    for(int j=0; j<seq_length; j++){
      temp_char = seq.seq[j];
      if((((temp_char=='A') || (temp_char=='a')) && (consensus(j)!=0))){
        m_i.push_back(n);
        m_j.push_back(j);
        m_x.push_back(1);
        n_snps += 1;
      } else if((((temp_char=='C') || (temp_char=='c')) && (consensus(j)!=1))){
        m_i.push_back(n);
        m_j.push_back(j);
        m_x.push_back(2);
        n_snps += 1;
      } else if((((temp_char=='G') || (temp_char=='g')) && (consensus(j)!=2))){
        m_i.push_back(n);
        m_j.push_back(j);
        m_x.push_back(3);
        n_snps += 1;
      } else if((((temp_char=='T') || (temp_char=='t')) && (consensus(j)!=3))){
        m_i.push_back(n);
        m_j.push_back(j);
        m_x.push_back(4);
        n_snps += 1;
      } else if((nt.find(temp_char) == std::string::npos && (consensus(j)!=4))){
        m_i.push_back(n);
        m_j.push_back(j);
        m_x.push_back(5);
        n_snps += 1;
      }
    }
    n += 1;
  }

  
  // cout << "loading snp vecs.. " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
  // begin_time = clock();

  uvec m_x_vec = conv_to< uvec >::from(m_x);
  uvec m_i_vec = conv_to< uvec >::from(m_i);
  uvec m_j_vec = conv_to< uvec >::from(m_j);

  
  // cout << "convert to uvec.. " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
  // begin_time = clock();


  uvec idsA = find(m_x_vec == 1); // Find indices
  umat locations(2, idsA.size());
  locations.row(0) = m_i_vec.elem(idsA).t();
  locations.row(1) = m_j_vec.elem(idsA).t();
  sp_umat sparse_matrix_A(locations, ones<uvec>(idsA.size()), n_seqs, seq_length);

  uvec idsC = find(m_x_vec == 2); // Find indices
  locations.set_size(2, idsC.size());
  locations.row(0) = m_i_vec.elem(idsC).t();
  locations.row(1) = m_j_vec.elem(idsC).t();
  sp_umat sparse_matrix_C(locations, ones<uvec>(idsC.size()), n_seqs, seq_length);

  uvec idsG = find(m_x_vec == 3); // Find indices
  locations.set_size(2, idsG.size());
  locations.row(0) = m_i_vec.elem(idsG).t();
  locations.row(1) = m_j_vec.elem(idsG).t();
  sp_umat sparse_matrix_G(locations, ones<uvec>(idsG.size()), n_seqs, seq_length);

  uvec idsT = find(m_x_vec == 4); // Find indices
  locations.set_size(2, idsT.size());
  locations.row(0) = m_i_vec.elem(idsT).t();
  locations.row(1) = m_j_vec.elem(idsT).t();
  sp_umat sparse_matrix_T(locations, ones<uvec>(idsT.size()), n_seqs, seq_length);

  uvec idsN = find(m_x_vec == 5); // Find indices
  locations.set_size(2, idsN.size());
  locations.row(0) = m_i_vec.elem(idsN).t();
  locations.row(1) = m_j_vec.elem(idsN).t();
  sp_umat sparse_matrix_N(locations, ones<uvec>(idsN.size()), n_seqs, seq_length);

  // cout << "convert to sparse.. " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
  // begin_time = clock();

  umat comp_snps = umat(sparse_matrix_A * sparse_matrix_A.t());
  comp_snps = comp_snps + umat(sparse_matrix_C * sparse_matrix_C.t());
  comp_snps = comp_snps + umat(sparse_matrix_G * sparse_matrix_G.t());
  comp_snps = comp_snps + umat(sparse_matrix_T * sparse_matrix_T.t());

  // cout << "calc similarity matrix.. " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
  // begin_time = clock();
  
  if(dist){
    umat diff_n = umat(sparse_matrix_N * sparse_matrix_N.t());
    comp_snps = comp_snps + diff_n;

    sp_umat total_snps = sum(sparse_matrix_A + sparse_matrix_C + sparse_matrix_G + sparse_matrix_T + sparse_matrix_N, 1);

    umat differing_snps = umat(n_seqs, n_seqs);
    for (int i=0 ; i<n_seqs; i++){
      for(int j=0; j<n_seqs; j++){
        differing_snps(i,j) = total_snps(i) + total_snps(j);
      }
    }
    sp_umat binary_snps = spones(sparse_matrix_A + sparse_matrix_C + sparse_matrix_G + sparse_matrix_T + sparse_matrix_N);

    if(count_n){
      comp_snps = differing_snps - umat(binary_snps * binary_snps.t()) - comp_snps;
    } else {
      umat total_n = umat(sum(sparse_matrix_N, 1));

      uvec cons_idsN = find(consensus == 4); // Find indices

      umat matrix_n_cols = zeros<umat>(n_seqs, cons_idsN.size());
      for (int i=0; i<cons_idsN.size(); i++){
        sp_umat col(binary_snps.col(cons_idsN(i)));
        for (arma::sp_umat::iterator it = col.begin(); it != col.end(); ++it) {
         matrix_n_cols(it.row(), i) = 1;
        }
      }

      // cout << "up to n matrix conversion.. " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
      // begin_time = clock();

      umat  cons_snps_N = matrix_n_cols * matrix_n_cols.t();
      uvec tot_cons_snps_N = uvec(sum(matrix_n_cols, 1));

      for (int i=0 ; i<n_seqs; i++){
        for(int j=0; j<n_seqs; j++){
          diff_n(i,j) = total_n(i) + total_n(j) - 2*diff_n(i,j) + tot_cons_snps_N(i) + tot_cons_snps_N(j) - 2*cons_snps_N(i,j);
        }
      }
      // cout << "calc diff.. " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
      // begin_time = clock();
      
      comp_snps = differing_snps - umat(binary_snps * binary_snps.t()) - comp_snps - diff_n;
    }
  }
  
  // cout << "calc remaining dist.. " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
  // begin_time = clock();

  // Output the distance matrix to stdout
  for (int j=0; j < n_seqs; j++) {
    printf("%s", seq_names[j].c_str());
    for (int i=0; i < n_seqs; i++) {
      printf("%c%d", sep, int(comp_snps(i,j)));
    }
    printf("\n");
  }

  return 0;

  }







