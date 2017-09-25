#' Initial test function
#'
#' Takes one parameter and does something completely useless
#' @param x the parameter to take
#' @return a useless numeric value
#' @export
initial_test <- function(x){
    message(paste0("You passed: ", x))
    return(runif(1,0,1))
}
