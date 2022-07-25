###################################################################
### Perferct Plasticity Approximation (Purves et al. 2009) ###

find_subplot <- function(x.sub, y.sub, max.x.sub){
  subplot <- (max.x.sub * (y.sub - 1)) + x.sub
  return(subplot)
}

assign_canopy_layers_ppa <- function(dbh, subplot.size){
  # Calculate height
  h <- d2h_nomax_martcano(dbh)
  # Calculate crown area
  CA <- d2CA(dbh)
  # Sort by height
  h.order <- order(h, decreasing = TRUE)
  # Cumulative crown area
  cumulative.CA <- cumsum(h.order)
  canopy.layer <- ifelse(cumulative.CA < subplot.size, 1, 2)
  canopy.layer[1] <- 1  # ensure tallest tree is in the canopy
  return(canopy.layer)
}

