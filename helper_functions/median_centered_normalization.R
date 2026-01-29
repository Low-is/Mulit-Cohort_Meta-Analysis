gene_medians <- function(expr_mat, genes_of_interest = NULL) {
  
  if (!is.null(genes_of_interest)) {
    genes <- intersect(genes_of_interest, rownames(expr_mat))
    apply(expr_mat[genes, , drop = FALSE], 1, median, na.rm = TRUE)
  } else {
    apply(expr_mat, 1, median, na.rm = TRUE)
  }
}
