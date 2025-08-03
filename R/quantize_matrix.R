#' Quantize matrix
#'
#' @description Quantize matrix
#' @param  X Input matrix.
#' @param as_sparse Convert the matrix to sparse format.
#' @inheritParams quantize_vector
#' @inheritParams to_sparse
#'
#' @export
#' @examples
#' set.seed(42)
#' X <-simulate_matrix(nrow=100, ncol=10)
#' Xq <- quantize_matrix(X)
quantize_matrix <- function(X,
                            n=40,
                            as_sparse=TRUE,
                            verbose=TRUE){
  ## Quantize matrix
  Xq <- apply(X, 2,
              FUN = quantize_vector,
              n = n) |>
    as.matrix() |>
    `row.names<-`(rownames(X))
  ## Convert to sparse matrix
  if(isTRUE(as_sparse)){
    Xq <- to_sparse(
      obj = Xq,
      verbose = verbose
    )
  }
  ## Return
  return(Xq)
}
