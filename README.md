# vaxedemic
> Optimal global allocation of vaccinations during a pandemic

[![Build Status](https://travis-ci.org/user/pkg.svg?branch=master)](https://travis-ci.org/user/pkg)

Testing link to vignette, if [needed](https://jameshay218.github.io/vaxedemic/inst/doc/usage.html).

## Running an example simulation

1. Clone this repository to a local directory (henceforth referred to as `local_dir`).
2. If running locally, make a subdirectory in the repository directory called `outputs`.  The simulation outputs will be saved here.  
3. If running on the cluster, make a directory in the network drive (henceforth referred to as `remote_dir`).  Make two subdirectories in `remote_dir`, called `outputs` and `data` respectively.  Copy everything in the `data` directory of this repository into `remote_dir/data`.
4. If running on the cluster, and you have not previously submitted an R job to the cluster, setup the interface between R and the cluster by following [these instructions](https://mrc-ide.github.io/didehpc/vignettes/quickstart.html)
5. Choose a username (doesn't have to be your DIDE username).  Edit [R/cluster_tools.R](https://jameshay218.github.io/vaxedemic/R/cluster_tools.R) so that your `didehpc` options (as in [the instructions](https://mrc-ide.github.io/didehpc/vignettes/quickstart.html)) are configured in `get_user_options`.  Make sure that `wd` and `package_dir` are set to `remote_dir` and `local_dir` respectively.
6. Edit [scripts/example.R](https://jameshay218.github.io/vaxedemic/scripts/example.R), replacing `user <- "ayan"` with your username, and `package_dir <- "~/Documents/vaxedemic/"` with your local directory.
7. Set `cluster <- TRUE` near the top of [scripts/example.R](https://jameshay218.github.io/vaxedemic/scripts/example.R) if running on the cluster, otherwise replace with `cluster <- FALSE`.
8. Set `run_fixed <- TRUE` near the top of [scripts/example.R](https://jameshay218.github.io/vaxedemic/scripts/example.R) if running for a fixed set of parameters, otherwise replace with `run_fixed <- FALSE` to run for many combinations of parameters.
9. Source [scripts/example.R](https://jameshay218.github.io/vaxedemic/scripts/example.R).


## Installation
```r
devtools::install_github("jameshay218/vaxedemic")
```

## To do to package structure, if needed
1. Integrate with Travis
2. Write test cases
3. Write full vignettes


## License

GPL-3 © [James Hay &lt;james.hay13@imperial.ac.uk&gt;](https://github.com/jameshay218).
