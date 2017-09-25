#' @useDynLib vaxedemic
#' @importFrom Rcpp evalCpp
.onUnload <- function(libpath){
    library.dynam.unload("vaxedemic",libpath)
}
