####################################################################
### JN July 2022 This is the body of the forest IBM. It simulates 
### forest dynamics from bare ground or stand initialisation using 
### user supplied growth, survival and recruitment algorithms. 
### Canopy dynamics use the PPA scheme. 
#################################################################################
rm(list = ls())

### LIBRARIES ### 
library(RColorBrewer)
library(plyr)

### FUNCTIONS ### 
source('Source/Allometry_functions.R')
source('Source/PPA_scheme.R')
source('Source/Demography_functions.R')
#################################################################################
### USER SETTINGS ###
cols <- brewer.pal(3, 'Dark2')

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
understory.mort.params <- c(0.1, 0.06)

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

# Wood density for each pft
wd <- c(0.3, 0.6)
if(length(wd) != no.pfts){
  stop('An estimate of wood density must be supplied for each pft.')
}


### End of user settings ### 
######################################################################################
t <- 1 # year 0
t.max <- 25 # how many years to simulate

# If we are starting from bare ground then build the initial matrix
Ns <- rep(NA, no.pfts)
for(pft in 1:no.pfts){
  Ns[pft] <- recruit.params[pft] * max.x * max.y
}
N <- sum(Ns)

x <- runif(N, min = 0, max = max.x)
y <- runif(N, min = 0, max = max.y)
pft <- rep(seq(no.pfts), unlist(Ns))
dbh <- rep(min.dbh, N)

# create the matrix that holds dbh through time
#mat <- cbind(pft, x, y, dbh.1)

# find the subplot for each tree
max.x.sub <- max.x/subplot.length
x.sub <- cut(x, breaks = seq(0,max.x,31.25), 
             labels = FALSE, include.lowest = TRUE)
y.sub <- cut(y, breaks = seq(0,max.y,31.25), 
             labels = FALSE, include.lowest = TRUE)
subplot <- find_subplot(x.sub, y.sub, max.x.sub)
mat <- cbind(pft, x, y, subplot, dbh)

# Assign canopy layers
dbh.list <- tapply(mat[ ,'dbh'], mat[ ,'subplot'], function(x) x)
canopy.layer <- do.call(c, lapply(dbh.list, assign_canopy_layers_ppa, 
                                    subplot.size = subplot.size))
mat <- cbind(mat, canopy.layer)
rownames(mat) <- NULL


########################################################################
# Now loop through each time step. Calculate canopy layer for each tree, 
# apply survival and then growth and recruitment. 

# To stop the matrix becoming too large from all the recruits, when a tree dies 
# we remove it from the 'alive' matrix and put it in the 'dead' matrix. New 
# recruits are added to the 'alive' matrix.
matrix.died <- vector('list', t.max)

# To prevent us growing the matrix we will overwrite dbh and canopy layer as we go but
# save dbh and canopy layer from each time point 
matrix.alive <- vector('list', t.max)

for(t in 1:t.max){
  
  matrix.alive[[t]] <- mat
  
  # apply mortality (NB 0s mean alive, 1 means dead)  
  mortality <- apply_mortality(mat[ ,'canopy.layer'],
                               mat[ ,'pft'], 
                               canopy.mort.params, understory.mort.params)
  matrix.died[[t]] <- mat[which(mortality == 1), ]
  mat <- mat[which(mortality == 0), ] 
  
  # apply growth
  ddbh <- apply_growth(mat[ ,'canopy.layer'], mat[ ,'pft'], 
                       canopy.growth.params, understory.growth.params)  
  mat[ ,'dbh'] <- mat[ ,'dbh'] + ddbh 
  
  # apply recruitment
  new.recruits <- apply_recruitment(recruit.params, max.x, max.y, max.x.sub, t)
  mat <- rbind(mat, new.recruits)
  
  # assign canopy layer
  dbh.list <- tapply(mat[ ,'dbh'], mat[ ,'subplot'], function(x) x)
  mat[ ,'canopy.layer'] <- do.call(c, lapply(dbh.list, assign_canopy_layers_ppa, 
                                        subplot.size = subplot.size))
  
  print(t)
  t <- t+1
  
}

################################################################################
### Combine all the trees that died into a single matrix ###
died <- do.call(rbind.fill.matrix, matrix.died)


par(mfrow = c(2,1), mar = c(3,2,1,1), oma = c(2,2,2,0))
# Number of trees through time total
N <- unlist(lapply(matrix.alive, nrow))
plot(seq(t.max), N, type = 'l', xlab = '', ylab = '')
mtext('Total Population Size', side = 3, line = 0.1, cex = 1.3)
# Number of trees through time x pft
N_x_pft <- lapply(matrix.alive, function(x) tapply(x[ ,'dbh'], x[ ,'pft'], length))
plot(seq(t.max), lapply(N_x_pft, '[[', 2), type = 'l', col = cols[2])  # ST
points(seq(t.max), lapply(N_x_pft, '[[', 1), type = 'l', col = cols[1])  # LD
mtext('Population Size by PFT', side = 3, line = 0.1, cex = 1.3)
mtext('Year', side = 1, line = 0, outer = TRUE, cex = 1.3)
mtext('N. plants', side = 2, line = 0, outer = TRUE, cex = 1.3)
legend('topleft', legend = c('LD', 'ST'), col = c(cols[2], cols[1]), lwd = 2)

# AGB through time x pft
agb.matrix <- matrix(NA, nrow = t.max, ncol = 2)

for(i in 1:length(matrix.alive)){
  tmp.mat <- matrix.alive[[i]]
  h <- d2h_nomax_martcano(tmp.mat[ ,'dbh']/10)
  wds <- wd[tmp.mat[ ,'pft']]
  agb <- d2agb(dbh = tmp.mat[ ,'dbh']/10, h = h, wd = wds)
  agb.matrix[i, ] <- tapply(agb, tmp.mat[ ,'pft'], sum)/plot.size/1000 # make it Mg per ha
  print(i)
}

par(mfrow = c(2,1), mar = c(3,2,1,1), oma = c(2,2,2,0))
# AGB through time total
AGB <- rowSums(agb.matrix)
plot(seq(t.max), AGB, type = 'l', xlab = '', ylab = '')
mtext('Total Aboveground biomass (Mg C ha-1)', side = 3, line = 0.1, cex = 1.3)
# AGB through time x pft
plot(seq(t.max), agb.matrix[ ,1], type = 'l', col = cols[1])  # LD
points(seq(t.max), agb.matrix[ ,2], type = 'l', col = cols[2])  # ST
mtext('Aboveground biomass by PFT (Mg C ha-1)', side = 3, line = 0.1, cex = 1.3)
mtext('Year', side = 1, line = 0, outer = TRUE, cex = 1.3)
mtext('N. plants', side = 2, line = 0, outer = TRUE, cex = 1.3)
legend('topleft', legend = c('LD', 'ST'), col = c(cols[2], cols[1]), lwd = 2)


