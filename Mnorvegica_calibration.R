######################################################################
##                                                                  ##
## Calibration procedure with Price algorithm                       ##
##                                                                  ##
## $Id$                                                             ##
######################################################################

require(matrixStats)
#require(ecolMod)

krill_cal_norv <- function(pars) {
  
  #--------------------------------------------------------------------#
  # Argument : pars
  #
  # Named vector that contains the initial parameters to calibrate :
  #

  #--------------------------------------------------------------------#
  # Load the required functions

  #setwd("/my/directory/")
  source("Mnorvegica.R")
  source("Pricefit.R")
  
  #--------------------------------------------------------------------#
  # Define the cost function
  
    sce_fun <- function(pars) {
    
    pars_sens = 0
    lg_ini    = 0
    
    result <- krill_fun_norv(pars, pars_sens, lg_ini)
    
    cw    <- result$cw_ind
    cmez  <- result$cp_mez
    xplot <- result$xplot

    end   <- dim(cw)[1]
    
    sce    <- sum( ( (colMeans(cw[(end-3279):end,], na.rm=T)   - 130 ) / 130 )^2, na.rm=T ) +
              sum( ( (colMeans(cmez[(end-3279):end,], na.rm=T) - 0.69) / 0.69)^2, na.rm=T )
    
    if (sum(sce)==0) sce = 1e9
    return(sce)
  }
  #--------------------------------------------------------------------#
  # Calibration - Price algorithm
  
  # Define the parameters needed by preicefit() (see details in Pricefit.R)
  #pref_food <- FALSE 
  
  # M. norvegica 
  if(sum(unlist(pars))==0) pars <- list( k0=c(0.5,0.1), h0=50, ei=0.2, er=0.5 )
  minpars  <- list( c(1e-2,1e-2), 1e1, 0.2, 0.5 )
  maxpars  <- list( c(   1,   1), 1e3, 0.4, 0.5 )
  
  npop     <- max(2*length(pars),20)
  numiter  <- 1000
  centroid <- 3
  varleft  <- 1e-4
  
  calib_norv <- pricefit(pars, minpars, maxpars, sce_fun, npop, numiter, centroid, varleft)
  write.csv2(calib_norv$par, file="params_calib_Mnorvegica.csv")
    
  return(calib_norv)

}
