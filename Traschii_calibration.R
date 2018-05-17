######################################################################
##                                                                  ##
## Calibration procedure with Price algorithm                       ##
##                                                                  ##
## $Id$                                                             ##
######################################################################

require(matrixStats)
#require(ecolMod)

krill_cal_rasch <- function(pars) {
  #--------------------------------------------------------------------#
  # Load the required functions

  #setwd("/my/directory/")
  source("Traschii.R")
  source("Pricefit.R")
  
  #--------------------------------------------------------------------#
  # Define the cost function
  
  sce_fun <- function(pars) {
    
    pars_sens = 0
    lg_ini    = 0
    
    result <- krill_fun_rasch(pars, pars_sens, lg_ini)
    
    cw    <- result$cw_ind
    cmez  <- result$cp_mez
    xplot <- result$xplot[1:4380]

    end   <- dim(cw)[1]
    
    sce    <- sum( ( (colMeans(cw[(end-3279):end,], na.rm=T)   - 17  ) / 17  )^2, na.rm=T ) +
              sum( ( (colMeans(cmez[(end-3279):end,], na.rm=T) - 0.28) / 0.28)^2, na.rm=T )
    
    if (sum(sce)==0) sce = 1e9
    return(sce)
  } 
  
  #--------------------------------------------------------------------#
  # Calibration - Price algorithm
  
  # Define the parameters needed by preicefit() (see details in Pricefit.R)
  #pref_food <- FALSE
  
  # T. raschii
  if(sum(unlist(pars))==0) pars <- list( k0=c(0.1,0.05), h0=50, ei=0.2, er=0.5)
  minpars  <- list( c(1e-2,1e-2), 1e1, 0.2, 0.5 )
  maxpars  <- list( c(   1,   1), 1e3, 0.4, 0.5 )
  
  npop     <- max(2*length(pars),20)
  numiter  <- 1000
  centroid <- 3
  varleft  <- 1e-4
  
  calib_rasch <- pricefit(pars, minpars, maxpars, sce_fun, npop, numiter, centroid, varleft)
  write.csv2(calib_rasch$par, file="params_calib_Traschii.csv")
  
  return(calib_rasch)
  
}
