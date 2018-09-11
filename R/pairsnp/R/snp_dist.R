#' snp_dist
#'
#' Function to calculate pairwise snp distance matrix
#'
#' @import Matrix
#'
#' @param sparse.data a sparse SNP data object returned from import_fasta_sparse
#' @param compare.n whether to count comparisons with gaps in distance calculation
#'
#' @return A pairwise snp distance matrix
#'
#' @examples
#' fasta.file.name <- system.file("extdata", "seqs.fa", package = "pairsnp")
#' sparse.data <- import_fasta_sparse(fasta.file.name)
#' snp_dist(sparse.data)
#'
#' @export
snp_dist <- function(sparse.data, compare.n=FALSE){

  # Check inputs
  if(!is.list(sparse.data)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!(class(sparse.data$snp.matrix)=="dgCMatrix")) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!is.numeric(sparse.data$consensus)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!is.logical(compare.n)) stop("compare.n must be TRUE/FALSE")

  n.isolates <- ncol(sparse.data$snp.matrix)

  shared.snps <- as.matrix(tcrossprod(t(sparse.data$snp.matrix==1)))
  shared.snps <- shared.snps + as.matrix(tcrossprod(t(sparse.data$snp.matrix==2)))
  shared.snps <- shared.snps + as.matrix(tcrossprod(t(sparse.data$snp.matrix==3)))
  shared.snps <- shared.snps + as.matrix(tcrossprod(t(sparse.data$snp.matrix==4)))
  shared.snps <- shared.snps + as.matrix(tcrossprod(t(sparse.data$snp.matrix==5)))

  total.snps <- colSums(sparse.data$snp.matrix>0)

  if(compare.n){
    differing.snps <- (matrix(rep(total.snps, n.isolates), nrow = n.isolates, byrow = TRUE) +
                         matrix(rep(total.snps, n.isolates), nrow = n.isolates, byrow = FALSE) -
                         as.matrix(tcrossprod(t(sparse.data$snp.matrix>0))) -
                         shared.snps)
  } else {
    total.n <- colSums(sparse.data$snp.matrix==5)

    n.cols.submatrix <- 1*(sparse.data$snp.matrix>0)[sparse.data$consensus==5,]
    cons.snps.N <- as.matrix(tcrossprod(t(n.cols.submatrix)))
    tot.cons.snps.N <- colSums(n.cols.submatrix)
    tot.cons.snps.N <- (matrix(rep(tot.cons.snps.N, n.isolates), nrow = n.isolates, byrow = TRUE) +
                           matrix(rep(tot.cons.snps.N, n.isolates), nrow = n.isolates, byrow = FALSE))

    diff.n <- as.matrix(tcrossprod(t(sparse.data$snp.matrix==5)))
    diff.n <- (matrix(rep(total.n, n.isolates), nrow = n.isolates, byrow = TRUE) +
                 matrix(rep(total.n, n.isolates), nrow = n.isolates, byrow = FALSE) -
                 2*diff.n +
                 tot.cons.snps.N -
                 2*cons.snps.N)
    differing.snps <- (matrix(rep(total.snps, n.isolates), nrow = n.isolates, byrow = TRUE) +
                         matrix(rep(total.snps, n.isolates), nrow = n.isolates, byrow = FALSE) -
                         as.matrix(tcrossprod(t(sparse.data$snp.matrix>0))) -
                         shared.snps -
                         diff.n)
  }

  return(differing.snps)
}
