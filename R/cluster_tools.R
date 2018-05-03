#' Cluster setup
#'
#' This function sets up an object linking to the DIDE cluster in a very crude way. Each user will need to implement their own version of this that returns a valid didehpc object to submit jobs to.
#' @param user a user identifier so that each user can specify their cluster options.
#' @return a didehpc::queue_didehpc object
#' @export
setup_cluster <- function(user){
    user_options <- get_user_options(user)
    setwd(user_options$wd)
    do.call(options, user_options$cluster_options)
    src <- provisionr::package_sources(local = user_options$package_dir, expire = 1e-10)
    sources <- NULL
    
    ## Setup contexts
    context::context_log_start()
    root <- "contexts"
    packages <- list(attached=c("vaxedemic","plyr","reshape2","data.table","ggplot2","Matrix","foreach"))
    ctx <- context::context_save(packages=packages,path=root, sources=sources,package_sources=src)
    
    ## Submit setup to cluster
    obj1 <- didehpc::queue_didehpc(ctx)
    return(obj1)
}

#' make argument list to run simulations on or off the cluster
#' 
#' @param runs data frame where each row is a set of parameters for which to run the simulation.
#' # if there are no variable parameters, set to NULL.
#' @param submit_fn character string.  name of function to run on cluster
#' @param obj if running on cluster, the output of setup_cluster.  if not running
#' on cluster, set to NULL.
#' @return an argument list to run simulations on or off the cluster
make_arg_list <- function(runs = NULL, submit_fn, obj = NULL) {
  args_submit_fn <- formalArgs(submit_fn)
  if(!is.null(runs)) {
    stopifnot(all(colnames(runs) %in% args_submit_fn))
    args_submit_fn <- args_submit_fn[!(args_submit_fn %in% colnames(runs))]
  }
  envir <- parent.frame()
  args_list <- list_vars_from_environment(args_submit_fn, envir = envir)
  if(!is.null(runs)) {
    if(is.null(obj)) {
      args_list <- lapply(seq_len(nrow(runs)), function(x) c(runs[x,], args_list))
    } else {
      args_list <- c(list(obj, runs, submit_fn, do_call = TRUE, timeout = 0), args_list)
    }
  }
  args_list
}

#' specify user options for cluster
#' 
#' @param user a user identifier so that each user can specify their cluster options.
#' @return a list with the working directory on the network drive; the directory
#' in which the package code sits; and options for didehpc
#' @export
get_user_options <- function(user) {
  # wd is the working directory to run the cluster job from. 
  # This should be the user's network home drive eg. "~/net/home/vaxedemic"
  if(user == "JH") {
    list(wd = "~/net/home/vaxedemic/",
         package_dir = "~/Documents/vaxedemic/",
         cluster_options = list(didehpc.credentials = "~/.smbcredentials",
                                        didehpc.cluster = "fi--didemrchnb"))
  } else if(user == "ayan") {
    list(wd = "~/net/home/vaxedemic2/",
         package_dir = "~/Documents/vaxedemic/",
         cluster_options = list(didehpc.username = "ayan",
                                didehpc.cluster = "fi--didemrchnb"))
  } else {
    stop("unknown user")
  }
}