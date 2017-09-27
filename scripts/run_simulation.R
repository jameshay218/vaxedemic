source("~/Documents/vaxedemic/R/simulation.R")

## LIFE HISTORY PARAMETER INPUTS
life_history_params <- list(R0=1.8,TR=2.6)

## SIMULATION OPTIONS
simulation_flags <- list(ageMixing=TRUE,
                         riskGroups=TRUE,
                         normaliseTravel=TRUE,
                         spatialCoupling=TRUE)
tmax <- 100
tdiv <- 24

## SETUP FAKE COUNTRY DATA
popn_size <- 100000
n_countries <- 10

## Setup age propns
n_ages <- 4
age_propns <- rep(1/n_ages, n_ages)
age_propns <- c(5,14,45,16)/80
## Setup risk groups
n_riskgroups <- 2

risk_propns <- rep(1/n_riskgroups,n_riskgroups) ## Assume risk groups are uniformly distributed
risk_propns <- matrix(rep(risk_propns,each=n_ages),ncol=n_riskgroups) ## Assume that proportion of ages in each risk group are the same for all ages
risk_factors <- rep(1, n_riskgroups) ## Assume that each risk group has same modifier

## Enumerate out risk factors for each age group
age_specific_riskgroup_factors <- matrix(rep(risk_factors,each=n_ages),
                                         ncol=n_riskgroups)

## Seeding setting
seedCountries <- c(1)
seedSizes <- c(10)
seedAges <- 3

## Contact rates
contactRates <- c(6.92,.25,.77,.45,.19,3.51,.57,.2,.42,.38,1.4,.17,.36,.44,1.03,1.83)
contactDur <- c(3.88,.28,1.04,.49,.53,2.51,.75,.5,1.31,.8,1.14,.47,1,.85,.88,1.73)


## Travel coupling
K <- matrix(1,n_countries,n_countries)+999*diag(n_countries) #Travel coupling - assumed independent of age (but can be changed)

tmp <- setup_populations(popn_size,n_countries,age_propns, n_ages,
                         risk_propns, risk_factors,
                         n_riskgroups)
X <- tmp$X
labels <- tmp$labels
    
## Generate a contact matrix with dimensions (n_ages*n_riskgroups) * (n_ages*n_riskgroups). ie. get age specific,
## then enumerate out by risk group. If we had country specific contact rates, we would need
## a matrix with the same dimensions as X
C1 <- generate_contact_matrix(contactRates, contactDur,n_ages, simulation_flags[["ageMixing"]])
C2 <- kronecker(C1, matrix(1,n_riskgroups,n_riskgroups))

## Generate risk factor modifier. ie. modifier for each age/risk group pair, same dimensions as C2
risk <- c(t(age_specific_riskgroup_factors))
risk_matrix <- t(kronecker(risk,matrix(1,1,n_riskgroups*n_ages)))

C3 <- C2*risk_matrix

## Normalise 

sim_params <- list(n_countries=n_countries,
                   n_ages=n_ages,
                   n_riskgroups=n_riskgroups,
                   seedCs=seedCountries,
                   seedNs=seedSizes,
                   seedAges=seedAges)
travelMatrix <- K
contactMatrix <- C3

res <- run_simulation(simulation_flags, life_history_params, sim_params,
                      X, C3, K, tmax, tdiv)




plot_labels <- expand.grid("Time"=seq(0,tmax,by=1/tdiv),"Location"=1:n_countries,"Age"=1:n_ages,"RiskGroup"=1:n_riskgroups)

I <- cbind(labels[,c("Location","Age","RiskGroup")], res$I)
I <- melt(I, id.vars=c("Location","Age","RiskGroup"))
I$Age <- as.factor(I$Age)
I$RiskGroup <- as.factor(I$RiskGroup)
I$variable <- as.numeric(I$variable)
times <- seq(0,tmax,by=1/tdiv)
I$variable <- times[I$variable]
I_aggregated <- aggregate(I[,"value"], I[,c("variable","Location","Age")], FUN=sum)

N <- aggregate(data=labels, X~Location + Age,FUN=sum)
I_aggregated <- merge(I_aggregated,N,id.vars=c("Location","Age"))

p1 <- ggplot(I_aggregated,aes(x=variable,y=x/X,col=Age)) +
    geom_line() +
    facet_wrap(~Location) +
    theme_bw()

p2 <- ggplot(I, aes(x=variable,y=value,col=RiskGroup)) + geom_line() + facet_grid(Age~Location) + theme_bw()

cowplot::plot_grid(p1,p2,ncol=2,align="hv")
