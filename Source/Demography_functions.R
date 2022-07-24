############################################################################
### Demography functions ###

### Survival
apply_mortality <- function(canopy.layer, pft, canopy.mort.params, 
                           understory.mort.params){
  prob <-  ifelse(canopy.layer == 1, 
                  canopy.mort.params[pft],
                  understory.mort.params[pft])
  mortality <- rbinom(length(pft), 1, prob)
  return(mortality)
}

### Growth
apply_growth <- function(canopy.layer, pft, canopy.growth.params, 
                         understory.growth.params){
  growth <- ifelse(canopy.layer == 1, 
                   canopy.growth.params[pft], 
                   understory.growth.params[pft])
  return(growth)
}

### Recruitment
# In this simplest of forms recruitment is independent of adult tree abundance - 
# i.e. we assume a seed source outside of the plot. 
apply_recruitment <- function(recruit.params, max.x, max.y, max.x.sub, t){
 
  # number of recruits for each pft
  Ns <- rep(NA, no.pfts)
  for(pft in 1:no.pfts){
    Ns[pft] <- recruit.params[pft] * max.x * max.y
  }
  N <- sum(Ns)
  
  x <- runif(N, min = 0, max = max.x)
  y <- runif(N, min = 0, max = max.y)
  pft <- rep(seq(no.pfts), unlist(Ns))
  
  # find the subplot for each tree
  x.sub <- cut(x, breaks = seq(0,max.x,31.25), 
               labels = FALSE, include.lowest = TRUE)
  y.sub <- cut(y, breaks = seq(0,max.y,31.25), 
               labels = FALSE, include.lowest = TRUE)
  subplot <- find_subplot(x.sub, y.sub, max.x.sub)
  
  padding <- matrix(NA, nrow = N, ncol = (t*2)+1)
  recruit.mat <- cbind(pft, x, y, subplot, padding)
  recruit.mat[ ,ncol(recruit.mat)] <- min.dbh
  
  return(recruit.mat)
  
}
