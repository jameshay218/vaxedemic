library(microbenchmark)
#sourceCpp("src/simulation_functions.cpp")
## Testing matrix multiplication speeds
#n <- 1600
#x <- matrix(runif(n*n, 0, 1000),ncol=n, nrow=n)
#y <- matrix(runif(n*n,0,1000),ncol=n,nrow=n)

#re1 <- microbenchmark(matrix_mult(x,y),times=10)

setwd("~/Documents/vaxedemic/")
devtools::document()
devtools::load_all()

##Inputs:
R0 <- 1
TR <- 2.6
gamma <- 1/TR
pop <- 100000*matrix(1,200,1)#Popuation
X <- pop
n <- length(X)
propRisk <- c(.2,.2,.2,.2)#Proportion of each age group at risk
risk <- matrix(1,4,1);#Risk multiplication factor - migth want to re-structure********
K <- matrix(1,n,n)#Travel coupling - assumed independent of age (but can be changed)
ndays <- 200#Number of days to simulate

#Properties:
ageProp <- matrix(c(5,14,45,16)/80,4,1)#Assumed uniform between 0 and 79 - option to change/make country dependent

##================================================================================================
##Heterogeneities:
#Age:
Cnum <- matrix(c(6.92,.25,.77,.45,.19,3.51,.57,.2,.42,.38,1.4,.17,.36,.44,1.03,1.83), nrow=4, ncol=4)
Cnum=t(Cnum)#Contact number
Cdur <- matrix(c(3.88,.28,1.04,.49,.53,2.51,.75,.5,1.31,.8,1.14,.47,1,.85,.88,1.73), nrow=4, ncol=4)
Cdur=t(Cdur)#Contact duration
C1 <- Cnum*Cdur#Age/risk matrix (4x4)

#Account for risk in dimensionality:
k2 <- matrix(1,2,2);
C2 <- kronecker(C1,k2);

#Risk (binary for each age group):
nonrisk <- matrix(1,4,1)-propRisk
Rx <- matrix(c(propRisk,nonrisk),4,2)
Rx <- t(Rx)
Rx <- matrix(Rx,8,1)*kronecker(ageProp,matrix(1,2,1))#a1r,a1x,a2r,a2x,...
X <- kronecker(Rx,X)#Population in each country/age/risk class
X <- trunc(X)
#X1 <- kronecker(t(ageProp),pop)#Population in each country/age class (in case it's useful) - collapse X back down

#Risk factor:
nonrisk <- matrix(1,4,1)
Rf <- matrix(c(risk,nonrisk),4,2)
Rf <- t(Rf)
Rf <- matrix(Rf,8,1)
Rrep <- matrix(1,1,8)
Rf <- kronecker(Rrep,Rf)
C <- C2*Rf#THE age/risk matrix

#Travel:
Krow <- rowSums(K)
K <- kronecker(matrix(1,1,n),matrix(Krow,n,1))
Kdelta <- kronecker(diag(8),K)
K1 <- kronecker(matrix(1,8,8),t(K))
KC <- kronecker(C,t(K))

#Vaccination:
V <- matrix(0,n,n)#Vaxine distribution

##================================================================================================
##Pre-simulation:
#Beta calcualtion (global):
#Don't need space if total contact isotropic - David to proove********
div=6#Time interval = day/div
gamma <- gamma/div
XX <- gamma*C
ev <- eigen(XX)
ev <- ev$values
Rstar <- max(abs(ev))
beta <- R0/Rstar

#Force of Infection - dual mobility
M <- K1%*%X#Normalisation factor (multiply)
Mm1 <- 1/M
Mm1[M==0] <- 0
M <- kronecker(matrix(1,1,8*n),M)

#Seeding:
seedC <- 1#Seed country
seedN <- 10#Seed number
sigma <- matrix(0,n*8,1)#Seed
sigma[(seedC-1)*8+3,1] <- seedN#Seed in an adult
I <- sigma
S <- X-I

Smat <- matrix(0,8*n,ndays*div)
Imat <- Smat
Smat[,1] <- S
Imat[,1] <- I
#wow <- stochastic_simulation(S, I, S, beta, gamma, Kdelta, K1, Mm1, div, ndays, 8, n)
