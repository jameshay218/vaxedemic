library(microbenchmark)
sourceCpp("src/tmpCpp.cpp")
## Testing matrix multiplication speeds
n <- 1600
x <- matrix(runif(n*n, 0, 1000),ncol=n, nrow=n)
y <- matrix(runif(n*n - 100,0,1000),ncol=n-100,nrow=n)

re1 <- microbenchmark(matrix_mult(x,y),times=10)
