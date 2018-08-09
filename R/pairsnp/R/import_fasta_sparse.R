#' import_fasta_sparse
#'
#' Imports a fasta file to a sparse matrix representing SNPs from the consensus
#'
#' @import Matrix
#'
#' @param fasta.file.name path to the fasta file
#'
#' @return A sparse matrix reprsentation of the SNPs along with a consensus
#'
#' @examples
#' fasta.file.name <- system.file("extdata", "seqs.fa", package = "fastbaps")
#' sparse.data <- import_fasta_sparse(fasta.file.name)
#'
#' @export
import_fasta_sparse <- function(fasta.file.name){

  # Check inputs
  if(!file.exists(fasta.file.name)) stop(paste("Can't locate file", fasta.file.name))

  snp.data <- pairsnp:::import_fasta_to_vector_each_nt(fasta.file.name)
  snp.data$seq.names <-  gsub("^>","",snp.data$seq.names)

  snp.matrix <- Matrix::sparseMatrix(i=snp.data$i,
                                     j=snp.data$j,
                                     x=snp.data$x,
                                     dims = c(snp.data$num.seqs, snp.data$seq.length),
                                     dimnames = list(snp.data$seq.names, 1:snp.data$seq.length))

  #remove conserved columns as they aren't needed
  conserved <- colSums(snp.matrix>0) == 0
  snp.matrix <- snp.matrix[,!conserved]
  snp.data$consensus <- snp.data$consensus[!conserved]

  return(list(snp.matrix=t(snp.matrix), consensus=snp.data$consensus))
}
