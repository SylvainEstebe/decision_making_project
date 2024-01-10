install.packages("pacman")
pacman::p_load(extraDistr, R2jags, parallel, ggpubr, ggplot2, tidyverse, hesim)

set.seed(1983) 

### NB! Don't forget to set your working directory
setwd("~/Code/decision_project/decision_making_project/analysis/parameter_recovery")

# defining a function for calculating the maximum of the posterior density (not exactly the same as the mode)
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}


#------ create task environment -------------------
#"..."

#-------test CC_corr_indiv function and jags script ---------

#---set params
ntrials <- 10
ngroups <- 280
groupSize <- 4
ntokens <- 20
vals <- seq(0,ntokens,1)

# -- Defining 
alpha <- matrix(runif(ngroups*groupSize,0,20), nrow=groupSize, byrow=TRUE) 
rho <- matrix(runif(ngroups*groupSize,0,1), nrow=groupSize, byrow=TRUE) 
omega <- matrix(runif(ngroups*groupSize,0,1), nrow=groupSize, byrow=TRUE)

# Sourcing in simulation function
source("CC_fun_individual.R")

# Simulate
CC_sims <- CC_fun_individual(alpha,
                         rho,
                         omega,
                         ngroups,
                         groupSize,
                         ntrials)

# Pull things out
c <- CC_sims$c
Ga <- CC_sims$Ga
Gb <- CC_sims$Gb
p <- CC_sims$p


par(mfrow=c(2,2))
plot(CC_sims$Gb)
plot(CC_sims$Ga)
plot(CC_sims$c)
plot(CC_sims$p)

# Set up jags and run jags model
# -- DEFINING DATA

data_list <- list(
  ngroups = ngroups,
  ntrials = ntrials, 
  groupSize = groupSize,
  Ga = Ga,
  c = c)

# -- DEFINING PARAMS
params <- c("omega", "rho", "alpha") # Add any other model parameters you are interested in

start_time = Sys.time()
samples <- jags(data = data_list,
                inits=NULL,
                parameters.to.save = params,
                model.file ="CC_individual.txt", 
                n.chains = 4,
                n.iter=15000, n.burnin=5000, n.thin=1) #, n.cluster=3)

end_time = Sys.time()
end_time - start_time

#print(samples)

alpha_recov <- array(NA, c(groupSize, ngroups))
rho_recov <- array(NA, c(groupSize, ngroups))
omega_recov <- array(NA, c(groupSize, ngroups))

X <- samples$BUGSoutput$sims.list

for (g in 1:ngroups) {
  
  for (s in 1:groupSize) {
    
    alpha_recov[s,g] <- MPD(X$alpha[,s,g])
    rho_recov[s,g] <- MPD(X$rho[,s,g])
    omega_recov[s,g] <- MPD(X$omega[,s,g])
    
  }
}

#vectorize
true_alpha <- as.vector(alpha)
true_rho <- as.vector(rho)
true_omega <- as.vector(omega)

recov_alpha <- as.vector(alpha_recov)
recov_rho <- as.vector(rho_recov)
recov_omega <- as.vector(omega_recov)

source("recov_plot.R")
# straight up recovery
pl1 <- recov_plot(true_alpha, recov_alpha, c("true alpha", "recov alpha"), 'smoothed linear fit')
pl2 <- recov_plot(true_rho, recov_rho, c("true rho", "recov rho"), 'smoothed linear fit')
pl3 <- recov_plot(true_omega, recov_omega, c("true omega", "recov omega"), 'smoothed linear fit')
plog = ggarrange(pl1, pl2, pl3)
plog
# recov as effect of other true parameters...
pl1 <- recov_plot(true_alpha, recov_rho, c("true alpha", "recov rho"), 'smoothed linear fit')
pl2 <- recov_plot(true_alpha, recov_omega, c("true alpha", "recov omega"), 'smoothed linear fit')
pl3 <- recov_plot(true_rho, recov_alpha, c("true rho", "recov alpha"), 'smoothed linear fit')
pl4 <- recov_plot(true_rho, recov_omega, c("true rho", "recov omega"), 'smoothed linear fit')
pl5 <- recov_plot(true_omega, recov_alpha, c("true omega", "recov alpha"), 'smoothed linear fit')
pl6 <- recov_plot(true_omega, recov_rho, c("true omega", "recov rho"), 'smoothed linear fit')
plot2 = ggarrange(pl1, pl2, pl3, pl4, pl5, pl6)
plot2
# recov as effect of other recov parameters...
pl1 <- recov_plot(recov_alpha, recov_rho, c("recov alpha", "recov rho"), 'smoothed linear fit')
pl2 <- recov_plot(recov_alpha, recov_omega, c("recov alpha", "recov omega"), 'smoothed linear fit')
pl3 <- recov_plot(recov_rho, recov_alpha, c("recov rho", "recov alpha"), 'smoothed linear fit')
pl4 <- recov_plot(recov_rho, recov_omega, c("recov rho", "recov omega"), 'smoothed linear fit')
pl5 <- recov_plot(recov_omega, recov_alpha, c("recov omega", "recov alpha"), 'smoothed linear fit')
pl6 <- recov_plot(recov_omega, recov_rho, c("recov omega", "recov rho"), 'smoothed linear fit')
plot3 = ggarrange(pl1, pl2, pl3, pl4, pl5, pl6)
plot3
