#' Simulate a matrix with random values
#' 
#' Simulate a \code{nrow} x \code{ncol} matrix with random values. 
#' 
#' @param nrow Number of rows in the matrix.
#' @param ncol Number of columns in the matrix.
#' @param row_prefix Prefix for row names.
#' @param col_prefix Prefix for column names.
#' @param func Function to generate random values (default is \code{stats::rnorm}).
#' @param seed Optional seed for reproducibility. 
#' @export
#' @examples
#' X <- simulate_matrix(nrow=100, ncol=10)
simulate_matrix <- function(nrow=100,
                            ncol=10,
                            row_prefix="gene_",
                            col_prefix="sample_",
                            func=stats::rnorm,
                            seed=NULL) {
    if (!is.null(seed)) {
        set.seed(seed)
    }
    X <- matrix(func(nrow*ncol), nrow=nrow, ncol = ncol)
    rownames(X) <- paste0(row_prefix, seq_len(nrow))
    colnames(X) <- paste0(col_prefix, seq_len(ncol))
    return(X)
}