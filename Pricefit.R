############################################### 
## Soetaert and Herman (2008)                ## 
## A practical guide to ecological modelling ## 
## Chapter 4.                                ## 
## Parameterisation                          ## 
## Chapter 4.4.3.                            ## 
## Pseudo-random search, a random-based      ## 
## minimisation routine                      ## 
###############################################

# The Price algorithm...
pricefit <- function ( par,                            # initial par estimates
                       minpar = rep(-1e8,length(par)), # minimal parameter values
                       maxpar = rep( 1e8,length(par)), # maximal parameter values
                       func,                           # function to minimise
                       npop = max(5*length(par),50),   # nr elements in population
                       numiter  = 100,                 # number of iterations
                       centroid = 3,                   # number of points in centroid
                       varleft  = 1e-8,                # relative variation upon stopping  
                       ...)                      
{
  
  # Initialization
  
cost    <- function (pars) func(pars,...)
minpar  <- unlist(minpar)
maxpar  <- unlist(maxpar)
npar    <- length(minpar)
tiny    <- 1e-8
varleft <- max(tiny,varleft)

populationpar <- matrix(nrow = npop, ncol = npar, byrow = TRUE,
                        data = minpar + runif(npar*npop)*rep((maxpar-minpar),npop))

#colnames(populationpar) <- names(par)
populationpar[1,] <- unlist(par)
populationcost    <- apply(populationpar, FUN=cost, MARGIN=1)
iworst            <- which.max(populationcost)
worstcost         <- populationcost[iworst]  


# Hybridisation phase

iter<-0

while ( iter<numiter & (max(populationcost)-min(populationcost))>(min(populationcost)*varleft) ) 
  {
  iter<-iter+1
  
  selectpar <- sample(1:npop,size=centroid)     # for cross-fertilisation
  mirrorpar <- sample(1:npop,size=1)            # for mirroring
  
  newpar    <- colMeans(populationpar[selectpar,]) # centroid
  newpar    <- 2*newpar-populationpar[mirrorpar,]  # mirroring
  
  newpar    <- pmin( pmax(newpar,minpar) ,maxpar)
  newcost   <- cost(newpar)
  
  if (newcost < worstcost)
  {
    populationcost[iworst] <- newcost
    populationpar[iworst,] <- newpar
    
    iworst    <- which.max(populationcost) # new worst member 
    worstcost <- populationcost[iworst]    
  }
} # end j loop

ibest    <- which.min(populationcost)
bestpar  <- populationpar[ibest,]
bestcost <- populationcost[ibest]

return (list(par=bestpar, cost=bestcost, poppar=populationpar, popcost=populationcost))
}
  
  
  
  
  
  
  
  
  
  