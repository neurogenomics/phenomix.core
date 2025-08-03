#' Iterate GSEA
#' 
#' Iterate Gene Set Enrichment Analysis (GSEA) across all traits and cell types.
#' @param use_quantiles Which quantiles in to use when filtering genes in each column of  
#' \code{xmat} and \code{ymat}. 
#' @inheritParams iterate_lm
#' @inheritParams set_cores
#' 
#' @return \code{data.table} of enrichment results. 
#' 
#' @export
#' @importFrom parallel mclapply 
#' @importFrom broom tidy
#' @importFrom data.table rbindlist
#' @importFrom dplyr mutate
#' @importFrom stats p.adjust
#' @importFrom methods slot is
#' @importFrom utils tail
#' @examples
#' set.seed(42)
#' 
#' ### x matrix (e.g. cell type expression specificity signatures)
#' xmat <- simulate_matrix(nrow = 200, ncol = 10, col_prefix="celltype_")
#'
#' ### y matrix (e.g. gene-phenotype association signatures)
#' ymat <- simulate_matrix(nrow = 100, ncol = 50, col_prefix="phenotype_") 
#'
#' res_gsea <- iterate_gsea(xmat = xmat,
#'                          ymat = ymat,
#'                          workers = 1)
iterate_gsea <- function(xmat,
                         ymat,
                         correction_method = "BH",
                         qvalue_thresh = .05,
                         quantize = list(x=10,
                                         y=10),
                         # Use the top 2 quantiles by default
                         use_quantiles = list(x=utils::tail(seq(quantize$x),2),
                                              y=utils::tail(seq(quantize$y),2)
                                              ),
                         progressbar = TRUE,
                         workers = 1,
                         verbose = TRUE) {
    requireNamespace("GeneOverlap")
    x <- y <- p <- q <- NULL;
    
    data.table::setDTthreads(threads = 1)
    BPPARAM <- set_cores(workers = workers,
                         progressbar = progressbar,
                         verbose = verbose)
    
    ## Filter data
    X_list <- filter_matrices(X_list = list(xmat=xmat,
                                            ymat=ymat),
                              verbose = verbose)
    xmat <- X_list$xmat
    ymat <- X_list$ymat
    
    ### Run lm for all combinations of xmat and ymat
    messager(
        "Running", formatC(ncol(xmat) * ncol(ymat), big.mark = ","),
        "tests:", formatC(ncol(xmat), big.mark = ","), "xmat columns x",
        formatC(ncol(ymat), big.mark = ","), "ymat columns.",v=verbose
    )
    ## Quantize data
    if(is.numeric(quantize$x)){
        xmat <- quantize_matrix(X=xmat,
                                n=quantize$x,
                                verbose = verbose)
    }
    if(is.numeric(quantize$y)){
        ymat <- quantize_matrix(X=ymat,
                                n=quantize$y,
                                verbose = verbose)
    }
    
   
    # Get feature intersect again after filtering and quantizing
    features <- feature_intersect(list(xmat=xmat,
                                       ymat=ymat),
                                    verbose = verbose)  

    
    gsea_res <- BiocParallel::bplapply(
        BPPARAM = BPPARAM,
        X = stats::setNames(seq(ncol(ymat)),
                            colnames(ymat)), 
        FUN = function(i) {
            
        tt <- colnames(ymat)[i]
        if(!progressbar){
            messager("-",tt,": (",i,"/",ncol(ymat),")",
                     parallel=TRUE)
        }
        
        lapply(colnames(xmat), function(ct) { 
            dat <- data.frame(
                x = xmat[features, ct],
                y = ymat[features, tt],
                row.names = features
            )
            res <- GeneOverlap::newGeneOverlap(
                listA = rownames(subset(dat, x %in% use_quantiles$x)),
                listB = rownames(subset(dat, y %in% use_quantiles$y))
            ) |>
                GeneOverlap::testGeneOverlap()
            res_df <- data.frame(
                xvar = ct,
                trait_genes = length(res@listA),
                celltype_genes = length(res@listB),
                intersection = length(res@intersection),
                union = length(res@union),
                genome.size = length(res@genome.size),
                odds.ratio = res@odds.ratio,
                Jaccard = res@Jaccard,
                p = res@pval
            )
            return(res_df)
        }) |> data.table::rbindlist()
    }) |> 
        data.table::rbindlist(idcol = "yvar")

    ### Multiple-testing correction
    gsea_res <- gsea_res |>
        dplyr::mutate(q = stats::p.adjust(p = p,
                                          method = correction_method))  
    ### Filter only sig results
    sig_res <- subset(gsea_res, q < qvalue_thresh)
    messager("\n", formatC(nrow(sig_res), big.mark = ","),
            "significant results @ ", 
            correction_method, "<", qvalue_thresh)
    ### Return FULL results (not just sig)
    return(gsea_res)
}
