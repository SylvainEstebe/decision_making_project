CC_fun_individual <- function(alpha,
                          rho,
                          omega,
                          ngroups,
                          groupSize,
                          ntrials) {
  
  Gb <- array(NA, c(groupSize,ntrials,ngroups))
  Ga <- array(NA, c(groupSize,ntrials,ngroups))
  p <- array(NA, c(groupSize,ntrials,ngroups))
  c <- array(NA, c(groupSize,ntrials,ngroups))
  
  for (g in 1:ngroups) {
    
    for (s in 1:groupSize) {
      
      #beliefs about others on first trial - gamma-poisson distribution
      Gb[s,1,g] <- rpois(1,alpha[s,g])
      
      # modelled preference and first contribution - see below
      p[s,1,g] <- (rho[s,g]*Gb[s,1,g])
      c[s,1,g] <- rpois(1,p[s,1,g])
      
    }
    
    Ga[s,1,g] <- mean(c[-s,1,g]) # not really necessary as we specify this in the trial loop
    
    #--------------- Implementation of CC model --------------------------------
    
    for (t in 2:ntrials) {
      
      for (s in 1:groupSize) {
        #- Actual group contribution (i.e. except oneself)
        Ga[s,t-1,g] <- mean(c[-s,t-1,g])
        
        #- Belief about group contribution
        Gb[s,t,g] <- ((1-omega[s,g])*(Gb[s,t-1,g]))+(omega[s,g]*(Ga[s,t-1,g]))
        
        #- Contribution preference, given belief and matching preference rho  
        p[s,t,g] <- rho[s,g]*Gb[s,t,g]
        
        #- Contribution as discrete sample from preferences
        c[s,t,g] <- rpois(1,p[s,t,g])
      }
      
    }
    
  }
  
  # calculating average actual group-wise contribution for the very last trial
  for (s in 1:groupSize) {
    
    Ga[s,t,] <- colMeans(c[-s,t,])
    
  }
  
  result <- list(c=c,
                 p=p,
                 Gb=Gb,
                 Ga=Ga)
  
  return(result)
  
}

