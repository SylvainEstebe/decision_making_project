seed_id = 1983
set.seed(seed_id)

install.packages("pacman")
pacman::p_load(R2jags, parallel, polspline, ggplot2, glue)

# defining a function for calculating the maximum of the posterior density (not exactly the same as the mode)
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}


##### Preprocessing #####

groupSize <- 4
ntrials <- 10
pi <- 1.6 # used to be 1.4, but the original paper (and Josh' preprint both say 1.6)
ntokens <- 20
vals <- seq(0,ntokens,1)

rawDat <- read.csv("~/Code/decision_project/decision_making_project/data/public_good/HerrmannThoeniGaechterDATA.csv", skip = 3) # Public goods game

#- create covariates in raw data matrix
# nation index
rawDat$nation <- c()
rawDat$nation[rawDat$city=="Melbourne"]=1
rawDat$nation[rawDat$city=="Minsk"]=2
rawDat$nation[rawDat$city=="Chengdu"]=3
rawDat$nation[rawDat$city=="Copenhagen"]=4
rawDat$nation[rawDat$city=="Bonn"]=5
rawDat$nation[rawDat$city=="Athens"]=6
rawDat$nation[rawDat$city=="Seoul"]=7
rawDat$nation[rawDat$city=="Samara"]=8
rawDat$nation[rawDat$city=="Zurich"]=9
rawDat$nation[rawDat$city=="St. Gallen"]=9
rawDat$nation[rawDat$city=="Istanbul"]=10
rawDat$nation[rawDat$city=="Nottingham"]=11
rawDat$nation[rawDat$city=="Dnipropetrovs'k"]=12
rawDat$nation[rawDat$city=="Boston"]=13


# Variable for worry index Data from Scored out of 1, a higher value indicates either greater worry or a greater experience of harm.
# https://wrp.lrfoundation.org.uk/2021-risk-indexes/
rawDat$worry <- c()
rawDat$worry[rawDat$city=="Melbourne"]=0.30
rawDat$worry[rawDat$city=="Minsk"]=0.30
rawDat$worry[rawDat$city=="Chengdu"]=0.31
rawDat$worry[rawDat$city=="Copenhagen"]=0.23
rawDat$worry[rawDat$city=="Bonn"]=0.37
rawDat$worry[rawDat$city=="Athens"]=0.47
rawDat$worry[rawDat$city=="Seoul"]=0.47
rawDat$worry[rawDat$city=="Samara"]=0.40
rawDat$worry[rawDat$city=="Zurich"]=0.32
rawDat$worry[rawDat$city=="St. Gallen"]=0.32
rawDat$worry[rawDat$city=="Istanbul"]=0.49
rawDat$worry[rawDat$city=="Nottingham"]=0.30
rawDat$worry[rawDat$city=="Dnipropetrovs'k"]=0.41
rawDat$worry[rawDat$city=="Boston"]=0.35

# extract every third line - data file has lines representing others responses and we don't need that
redDat <- rawDat[seq(1,length(rawDat$sessionid),3),]
group_names <- unique(redDat$groupid)
ngroups <- length(group_names)

# THIS WILL REMOVE SUBJECTS WITH MISSING DATA IN NO PUNISHMENT CONDITION
ngroups <- 269

subject_names <- unique(redDat$subjectid)
nsubjects <- length(subject_names)

# data for no punishment condition #
c_no_punish <- array(0,c(groupSize,ntrials,ngroups)) # choices
Gga_no_punish <- array(0,c(ntrials,ngroups)) # group-averaged contribution (only 1 entry per group - hence the G + ga)
Ggas_no_punish <- array(0,c(groupSize,ntrials,ngroups)) # same as Gga, but specified for each subject (cuz that's how the JAGS-code wants it) - hence the s
Gc_no_punish <- array(0,c(groupSize,ntrials,ngroups)) # summed ("cumulated" hence the c) group contribution not including oneself (therefore specified for each subject) - we don't use this - this refers to the sum-command in the loop-structure
Ga_no_punish <- array(0,c(groupSize,ntrials,ngroups)) # group-averaged contribution without oneself - we don't use this cuz the participants don't see this (hence, they have to do quite a bit of mental arithmetics to represent this), and thus we're modeling their conditional preference relative to the averaged group contribution (including their own)

missing <- array(0,ngroups)

for (g in 1:ngroups) {
  c_no_punish[,,g] <- rbind(redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][1:10],
                            redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][11:20],
                            redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][21:30],
                            redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][31:40])
  
  Gga_no_punish[,g] <- colMeans(c_no_punish[,,g])
  Gga_no_punish
  missing[g] <- is.na(c_no_punish[1,1,g])
  
  for (s in 1:groupSize) {
    Gc_no_punish[s,,g] <- colSums(c_no_punish[-s,,g])
    Ga_no_punish[s,,g] <- colMeans(c_no_punish[-s,,g])
    Ggas_no_punish[s,,g] <- colMeans(c_no_punish[,,g])
  }
}

# data for punishment condition #
c_punish <- array(0,c(groupSize,ntrials,ngroups)) # choices
Gga_punish <- array(0,c(ntrials,ngroups)) # group-averaged contribution (only 1 entry per group - hence the G + ga)
Ggas_punish <- array(0,c(groupSize,ntrials,ngroups)) # same as Gga, but specified for each subject (cuz that's how the JAGS-code wants it) - hence the s
Gc_punish <- array(0,c(groupSize,ntrials,ngroups)) # summed ("cumulated" hence the c) group contribution not including oneself (therefore specified for each subject) - we don't use this - this refers to the sum-command in the loop-structure
Ga_punish <- array(0,c(groupSize,ntrials,ngroups)) # group-averaged contribution without oneself - we don't use this cuz the participants don't see this (hence, they have to do quite a bit of mental arithmetics to represent this), and thus we're modeling their conditional preference relative to the averaged group contribution (including their own)

for (g in 1:ngroups) {
  c_punish[,,g] <- rbind(redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="P-experiment"][1:10],
                         redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="P-experiment"][11:20],
                         redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="P-experiment"][21:30],
                         redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="P-experiment"][31:40])
  
  Gga_punish[,g] <- colMeans(c_punish[,,g])
  
  for (s in 1:groupSize) {
    Gc_punish[s,,g] <- colSums(c_punish[-s,,g])
    Ga_punish[s,,g] <- colMeans(c_punish[-s,,g])
    Ggas_punish[s,,g] <- colMeans(c_no_punish[,,g])
  }
}

# compile data from each condition into 4D matrix
c <- array(0,c(groupSize,ntrials,ngroups,2))
c[,,,1] <- c_no_punish
c[,,,2] <- c_punish
Gga <- array(0,c(ntrials,ngroups,2))
Gga[,,1] <- Gga_no_punish
Gga[,,2] <- Gga_punish
Ggas <- array(0,c(groupSize,ntrials,ngroups,2))
Ggas[,,,1] <- Ggas_no_punish
Ggas[,,,2] <- Ggas_punish

Gc <- array(0,c(groupSize,ntrials,ngroups,2))
Gc[,,,1] <- Gc_no_punish
Gc[,,,2] <- Gc_punish


Ga <- array(0,c(groupSize,ntrials,ngroups,2))
Ga[,,,1] <- Ga_no_punish
Ga[,,,2] <- Ga_punish

c_choice_index <- c

worry <- array(0, ngroups)
Nation <- array(0, ngroups)

for (g in 1:ngroups) {
  worry[g] <- mean(redDat$worry[redDat$groupid==group_names[g]&redDat$p=="N-experiment"])
  Nation[g] <- mean(redDat$nation[redDat$groupid==group_names[g]&redDat$p=="N-experiment"])
}

c_win <- c_no_punish[,,!is.na(worry)]
c_keep <- rep()-c_win

Gga_punish <- Gga_punish[,!is.na(worry)]
Gga_no_punish <- Gga_no_punish[,!is.na(worry)]

c <- c[,,!is.na(worry),]
Gga <- Gga[,!is.na(worry),]
Ggas <- Ggas[,,!is.na(worry),]
Gc <- Gc[,,!is.na(worry),]
Ga <- Ga[,,!is.na(worry),]
worry <- worry[!is.na(worry)]
Nation <- Nation[!is.na(Nation)]

#redefine number of groups after removing those without civic scores
ngroups <- length(worry)

# aggregate worry to just 1 number per Nation-index (using the mean here should be unproblematic since all groups within a given nation should have been given the same worry-coefficient)
worry <- aggregate(worry~Nation, FUN=mean)[,2]
worry
nnations <- length(worry)

# calculate the winnings (i.e. apply the multiplication-factor to the sum of each groups contributions)
winnings <-  array(0, ngroups)
for (g in 1:ngroups) {
  winnings[g] <- sum(colSums(c_win[,,g])*pi)
}

################################################################################
########################### Conditional cooperation model ######################
################################################################################

# JZS priors for partial correlation. Method described here
# https://link.springer.com/article/10.3758/s13423-012-0295-x
# Code available here
# https://github.com/MicheleNuijten/BayesMed/blob/master/R/jzs_corSD.R
# Paper where code is used here (mediation paper)
# https://link.springer.com/article/10.3758/s13428-014-0470-2

#################################################################
#------------------ Winnings analysis ---------------------------
#################################################################

#-------------------  Regress worry on winnings ---------------

# standardise variables

X <- worry
X <- (X-mean(X))/sd(X)

invSigma <- solve(t(X)%*%X) # required for JZS priors

Y <- (winnings-mean(winnings))/sd(winnings)

data <- list("ngroups", "Y", "nnations","X","Nation","invSigma") 
params <- c("beta0","betaX") 

# - run jags code
win.samples <- jags.parallel(data, inits=NULL, params,
                             model.file ="~/Code/decision_project/decision_making_project/analysis/win_corr.txt",
                             n.chains=3, n.iter=15000, n.burnin=5000, n.thin=1, n.cluster=3)

#traceplot(win.samples)

# diagnostics for convergence
#gelman.diag(as.mcmc(CC.samples))

#################################################################
#------------------ CC model analysis ---------------------------
#################################################################

#-------------------  Regress worry on belief weights and slope of prefs in CC model ---------------

# standardise covariate
X <- worry
X <- (X-mean(X))/sd(X)
invSigma <- solve(t(X)%*%X) # required for JZS priors

data <- list("groupSize", "ngroups", "ntrials", "nnations","c","Ga","X","Nation","invSigma") 
params <- c("beta0_alpha","betaX_alpha","beta0_rho","betaX_rho","beta0_omega","betaX_omega") 

# - run jags code
start_time = Sys.time()
CC.samples <- jags.parallel(data, inits=NULL, params,
                            model.file ="~/Code/decision_project/decision_making_project/analysis/CC_corr.txt",
                            n.chains=3, n.iter=15000, n.burnin=5000, n.thin=1, n.cluster=3)
end_time = Sys.time()
end_time - start_time

save(CC.samples, file= "~/Code/decision_project/decision_making_project/analysis/saved_JAGS/CC_worry.samples.Rdata")

#load("~/Code/decision_project/decision_making_project/analysis/saved_JAGS/CC_worry.samples.Rdata")

# diagnostics for convergence
gelman.diag(as.mcmc(CC.samples))
# Potential scale reduction factors:
#Point est. Upper C.I.
#beta0_alpha       1.10       1.25
#beta0_omega       1.03       1.04
#beta0_rho         1.05       1.14
#betaX_alpha       1.20       1.56
#betaX_omega       1.03       1.03
#betaX_rho         1.09       1.25
#deviance          1.02       1.08

#Multivariate psrf

#1.16
#################################################################
#------------------ Plotting ---------------------------
#################################################################

# ------ Create empirical and parameter arrays for plots -------
# empirical group winnings data - means and standard deviations
empirical.win <- array(0,c(3,length(worry)))
for (i in 1:length(worry)) {
  empirical.win[1,i] <- mean(winnings[Nation==i]) - sd(winnings[Nation==i]) 
  empirical.win[2,i] <- mean(winnings[Nation==i]) 
  empirical.win[3,i] <- mean(winnings[Nation==i]) + sd(winnings[Nation==i])
}


### PLOT
#----- empirical correlation - worry and group winnings --------------------
plot(c(0,0.6), c(0,1200), type = "n", main = "F: Data - Winnings", 
     xlab = "Worry National Index", ylab = "Group winnings",axes=FALSE)
for (i in 1:length(worry)) {
  lines(c(worry[i],worry[i]),c(empirical.win[1,i],empirical.win[3,i]))
  points(worry[i],empirical.win[2,i])
}
axis(1)
axis(2)
grid()


#----- posterior of effect - worry and group winnings --------------------
plot(density(win.samples$BUGSoutput$sims.list$betaX),frame=FALSE,lwd=2,ylim=c(0,7),
     cex=2,xlab = "Standardised Effect of Worry National Index",ylab = " Posterior distribution",main="G: Winnings")
COL <- adjustcolor(c("red"))
lines(c(win.samples$BUGSoutput$summary[2,3],win.samples$BUGSoutput$summary[2,7]),c(-.1,-.1),col=COL,lwd=2)
points(win.samples$BUGSoutput$summary[2,1],c(-.1),pch=19,col=COL)
abline(h=0,v=0,col="gray60")
par(mar=c(5,3,1,2)) # reset margin because no title



######

dens <- density(win.samples$BUGSoutput$sims.list$betaX)
x_range <- c(-0.10, 0.10)

plot(dens, frame=FALSE, ylim=c(0, 12),
     cex=2, xlab="Standardised Effect of Worry National Index", ylab="Posterior Distribution",
     main=expression(paste("G: Correlation")))

polygon(c(dens$x, rev(dens$x)), c(dens$y, rep(0, length(dens$y))), col=adjustcolor("skyblue", alpha.f=0.3), border=NA)

COL <- adjustcolor("red")
lines(c(win.samples$BUGSoutput$summary[2,3], win.samples$BUGSoutput$summary[2,7]), c(-.1,-.1), col=COL, lwd=2)
points(win.samples$BUGSoutput$summary[2,1], c(-.1), pch=19, col=COL)

legend("topright", legend=c("Posterior Distribution"," ", "95% Bayesian Credible Intervals"),
       col=c("black", "skyblue", "red"), lty=c(1, 1, 1), lwd=c(2, NA, 2))


# Réinitialiser la marge
par(mar=c(5,3,1,2))

######

#---------- Initial Belief -------------------------------------------------------

plot(density(CC.samples$BUGSoutput$sims.list$betaX_alpha),frame=FALSE,lwd=2,ylim=c(0,22),
     cex=2,xlab = "Standardised Effect of Worry National Index",ylab = " Posterior Distribution ",main=expression(paste("H: Initial Belief  ", alpha, "")))
COL <- adjustcolor(c("red"))
lines(c(CC.samples$BUGSoutput$summary[4,3],CC.samples$BUGSoutput$summary[4,7]),c(-.1,-.1),col=COL,lwd=2)
points(CC.samples$BUGSoutput$summary[4,1],c(-.1),pch=19,col=COL)
abline(h=0,v=0,col="gray60")

par(mar=c(5,3,1,2)) # reset margin because no title

######

dens <- density(CC.samples$BUGSoutput$sims.list$betaX_alpha)
x_range <- c(-0.10, 0.10)

plot(dens, frame=FALSE, ylim=c(0, 12),
     cex=2, xlab="Standardised Effect of Worry National Index", ylab="Posterior Distribution",
     main=expression(paste("H: Initial Belief  ", alpha, "")))

polygon(c(dens$x, rev(dens$x)), c(dens$y, rep(0, length(dens$y))), col=adjustcolor("skyblue", alpha.f=0.3), border=NA)

COL <- adjustcolor("red")
lines(c(CC.samples$BUGSoutput$summary[4,3], CC.samples$BUGSoutput$summary[4,7]), c(-.01,-.01), col=COL, lwd=2)
points(CC.samples$BUGSoutput$summary[4,1], c(-.01), pch=19, col=COL)

legend("topright", legend=c("Posterior Distribution"," ", "95% Bayesian Credible Intervals"),
       col=c("black", "skyblue", "red"), lty=c(1, 1, 1), lwd=c(2, NA, 2))


# Réinitialiser la marge
par(mar=c(5,3,1,2))

######

#---------- Belief Learning -------------------------------------------------------

plot(density(CC.samples$BUGSoutput$sims.list$betaX_omega),frame=FALSE,lwd=2,ylim=c(0,10),
     cex=2,xlab = "Standardised Effect of Worry National Index",ylab = " Posterior Distribution ",main="I: Belief Learning Weight - model parameter (omega)")
COL <- adjustcolor(c("red"))
lines(c(CC.samples$BUGSoutput$summary[5,3],CC.samples$BUGSoutput$summary[5,7]),c(-.1,-.1),col=COL,lwd=2)
points(CC.samples$BUGSoutput$summary[5,1],c(-.1),pch=19,col=COL)
abline(h=0,v=0,col="gray60")

par(mar=c(5,3,1,2))


####
dens <- density(CC.samples$BUGSoutput$sims.list$betaX_omega)

plot(dens, frame=FALSE, ylim=c(0, 12),
     cex=2, xlab="Standardised Effect of Worry National Index", ylab="Posterior Distribution",
     main=expression(paste("I: Belief Learning Weight  ", omega, "")))

polygon(c(dens$x, rev(dens$x)), c(dens$y, rep(0, length(dens$y))), col=adjustcolor("skyblue", alpha.f=0.3), border=NA)

COL <- adjustcolor("red")
lines(c(CC.samples$BUGSoutput$summary[5,3], CC.samples$BUGSoutput$summary[5,7]), c(-.1,-.1), col=COL, lwd=2)
points(CC.samples$BUGSoutput$summary[5,1], c(-.1), pch=19, col=COL)

legend("topright", legend=c("Posterior Distribution","", "95% Bayesian Credible Intervals"),
       col=c("black", "skyblue", "red"), lty=c(1, 1, 1), lwd=c(2, NA, 2))


# Réinitialiser la marge
par(mar=c(5,3,1,2))

####


#---------- Preference Slope -------------------------------------------------------

plot(density(CC.samples$BUGSoutput$sims.list$betaX_rho),frame=FALSE,lwd=2,ylim=c(0,10),
     cex=2,xlab = "Standardised Effect of Worry National Index",ylab = "Posterior Distribution ",main="J: Conditional Preferences - model parameter (rho)")
COL <- adjustcolor(c("black"))
lines(c(CC.samples$BUGSoutput$summary[6,3],CC.samples$BUGSoutput$summary[6,7]),c(-.1,-.1),col=COL,lwd=2)
points(CC.samples$BUGSoutput$summary[6,1],c(-.1),pch=19,col=COL)
abline(h=0,v=0,col="gray60")

par(mar=c(5,3,1,2))


##

####
dens <- density(CC.samples$BUGSoutput$sims.list$betaX_rho)

plot(dens, frame=FALSE, ylim=c(0, 12),
     cex=2, xlab="Standardised Effect of Worry National Index", ylab="Posterior Distribution",
     main=expression(paste("J: Conditional Preferences  ", rho, "")))

polygon(c(dens$x, rev(dens$x)), c(dens$y, rep(0, length(dens$y))), col=adjustcolor("skyblue", alpha.f=0.3), border=NA)

COL <- adjustcolor("red")
lines(c(CC.samples$BUGSoutput$summary[6,3], CC.samples$BUGSoutput$summary[6,7]), c(-.1,-.1), col=COL, lwd=2)
points(CC.samples$BUGSoutput$summary[6,1], c(-.1), pch=19, col=COL)

legend("topright", legend=c("Posterior Distribution","", "95% Bayesian Credible Intervals"),
       col=c("black", "skyblue", "red"), lty=c(1, 1, 1), lwd=c(2, NA, 2))


# Réinitialiser la marge
par(mar=c(5,3,1,2))

####

###