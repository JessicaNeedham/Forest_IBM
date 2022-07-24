####################################################################
### JN July 2022 This is the body of the forest IBM. It simulates 
### forest dynamics from bare ground or stand initialisation using 
### user supplied growth, survival and recruitment algorithms. 
### Canopy dynamics use the PPA scheme. 
#################################################################################
rm(list = ls())

### LIBRARIES ### 
library(RColorBrewer)

### FUNCTIONS ### 
source('Source/Allometry_functions.R')
source('Source/PPA_scheme.R')
source('Source/Demography_functions.R')
#################################################################################
### USER SETTINGS ###

# Should we spin up from bare ground? If not user must supply data for
# stand initialisation. This data must be in the forestnet format.
bare.ground <- TRUE

# Is the user supplying a recruitment function? If not then we use the default. 
set.recruitment <- FALSE

# How many PFTs to use?
no.pfts <- 2

# Plot size (in ha) and dimensions (in m)
plot.size <- 10
max.x <- 500
max.y <- 250
# subplot size (here 31.25 m to match Bohlman and Pacala xxxx)
subplot.length <- 31.25
subplot.size <- subplot.length^2

# Minimum dbh to simulate
min.dbh <- 1

# For each PFT supply mean canopy and understory growth rates (dbh increment in 
# mm per year). 
# (Eventually there will be more flexibility to accommodate different 
# growth algorithms).
canopy.growth.params <- c(6, 3)
understory.growth.params <- c(2, 1)

if(length(canopy.growth.params) != no.pfts){
  stop('A mean canopy growth rate must be supplied for each pft.')
}
if(length(understory.growth.params) != no.pfts){
  stop('A mean understory growth rate must be supplied for each pft.')
}

# For each PFT supply mean canopy and understory mortality rates (% of trees dead per year)
# This is currently stochastic - each tree has a draw from a binomial with prob given by 
# the supplied rates. 
canopy.mort.params <- c(0.04, 0.01)
understory.mort.params <- c(0.07, 0.03)

if(length(canopy.mort.params) != no.pfts){
  stop('A mean canopy mortality rate must be supplied for each pft.')
}
if(length(understory.mort.params) != no.pfts){
  stop('A mean understory mortality rate must be supplied for each pft.')
}

# if set recruitment is true then the user needs to add a recruitment 
# algorithm and function of there own 
if(set.recruitment == FALSE){
  # For each PFT supply mean recruitment rates 
  # (number of new recruits per m2)
  recruit.params <- c(1, 1)
} else {
  # User supplies recruitment algorithm here.... 
}

### End of user settings ### 
######################################################################################
t <- 1 # year 0
t.max <- 100 # how many years to simulate

# If we are starting from bare ground then build the initial matrix
Ns <- rep(NA, no.pfts)
for(pft in 1:no.pfts){
  Ns[pft] <- recruit.params[pft] * max.x * max.y
}
N <- sum(Ns)

x <- runif(N, min = 0, max = max.x)
y <- runif(N, min = 0, max = max.y)
pft <- rep(seq(no.pfts), unlist(Ns))
dbh.1 <- rep(min.dbh, N)

# create the matrix that holds dbh through time
#mat <- cbind(pft, x, y, dbh.1)

# find the subplot for each tree
max.x.sub <- max.x/subplot.length
x.sub <- cut(x, breaks = seq(0,max.x,31.25), 
             labels = FALSE, include.lowest = TRUE)
y.sub <- cut(y, breaks = seq(0,max.y,31.25), 
             labels = FALSE, include.lowest = TRUE)
subplot <- find_subplot(x.sub, y.sub, max.x.sub)
mat <- cbind(pft, x, y, subplot, dbh.1)

# Assign canopy layers
dbh.list <- tapply(mat[ ,'dbh.1'], mat[ ,'subplot'], function(x) x)
canopy.layer.1 <- do.call(c, lapply(dbh.list, assign_canopy_layers_ppa, subplot.size = subplot.size))
mat <- cbind(mat, canopy.layer.1)
rownames(mat) <- NULL

########################################################################
# Now loop through each time step. Calculate canopy layer for each tree, 
# apply survival and then growth and recruitment. 

# To stop the matrix becoming too large from all the recruits, when a tree dies 
# we remove it from the 'alive' matrix and put it in the 'dead' matrix. New 
# recruits are added to the 'alive' matrix.
matrix.died <- vector('list', t.max)

for(t in 1:t.max){
  
  # apply mortality (NB 0s mean alive, 1 means dead)  
  mortality <- apply_mortality(mat[ ,paste0('canopy.layer.', t)], mat[ ,'pft'], 
                               canopy.mort.params, understory.mort.params)
  matrix.died[[t]] <- mat[which(mortality == 1), ]
  mat <- mat[which(mortality == 0), ] 
  
  # apply growth
  ddbh <- apply_growth(mat[ ,paste0('canopy.layer.', t)], mat[ ,'pft'], 
                       canopy.growth.params, understory.growth.params)  
  dbh.new <- mat[ ,paste0('dbh.', t)] + ddbh 
  mat <- cbind(mat, dbh.new)
  colnames(mat)[ncol(mat)] <- paste0('dbh.', t+1)
  
  # apply recruitment
  new.recruits <- apply_recruitment(recruit.params, max.x, max.y, max.x.sub, t)
  mat <- rbind(mat, new.recruits)
  
  # assign canopy layer
  dbh.list <- tapply(mat[ ,paste0('dbh.', t+1)], mat[ ,'subplot'], function(x) x)
  canopy.layer.new <- do.call(c, lapply(dbh.list, assign_canopy_layers_ppa, 
                                        subplot.size = subplot.size))
  mat <- cbind(mat, canopy.layer.new)
  colnames(mat)[ncol(mat)] <- paste0('canopy.layer.', t+1)
  
  print(t)
  t <- t+1
  
}



