
##================================================================================================

function1 <- function(beta, Kdelta, K1, I, Mm1, n, S, gamma){
  ##Simulation:
  tend <- ndays*div
  for (i in 2:tend){
    #lambda <- matrix_mult(beta*Kdelta, matrix_mult(K1,I)*Mm1)
    lambda <- beta*Kdelta%*%((K1%*%I)*Mm1)
    Pinf <- 1-exp(-lambda)
    infect <- rbinom(8*n,S,Pinf)
    S <- S-infect
    I=I+infect
    Prec <- 1-exp(-gamma)
    recover <- rbinom(8*n,I,Prec)
    I <- I-recover
    Smat[,i] <- S
    Imat[,i] <- I
  }
}

res1 <- microbenchmark(function1(beta,Kdelta,K1,I,Mm1,n,S,gamma),times=1)
res2 <- microbenchmark(stochastic_simulation(S, I, S, beta, gamma, Kdelta, K1, Mm1, div, ndays, 8, n),times=1)


##================================================================================================
##Crude plot:
Iplot <- colSums(Imat)
plot(seq(1,ndays*div,by=1),Imat[3,])
