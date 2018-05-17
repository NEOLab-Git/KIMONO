############################################################################################
##                                                                                        ##
##              Individual Growth Model of Meganyctiohanes norvegica                      ##
##                  based on the energetic budget of the organism                         ##
##                                                                                        ##
## $Id$                                                                                   ##
############################################################################################

krill_fun_rasch <- function (pars, pars_sens=0, lg_ini=0) {
  
  #----------------------------#
  #   Parameters description   #
  #----------------------------#
  
  # pars        : list of parameters
  # nb_day      : Time duration of the simulation in day  
  # time_moult  : Time of the moult in day (Relation to the temperatrure from Cuzin and Buchholz 1999)
  # food        : Concentration of phytoplankton or zooplankton in mgC.l^-1 (netCDF data)
  # temp        : Temperature in °C (netCDF data)
  # aw          : The regression constant of length(mm)/mass(mgC) relationship
  # bw          : The regression coefficient of length(mm)/mass(mgC) relationship 
  # A           : Assimilation efficiency coefficent in % (From Fowler 1971, McClatchie 1985, Tarling 2000)
  # w_molt      : Percentage of mass loss due to moulting exuvie in % of mass (From Sameoto 1976)
  # ei          : The activation energy in eV (From Vettes 1997, Salomon et Saborowski 2006) 
  # T0          : Reference temperature (the frozing point of water) in K (From Gilloly et al. 2001)
  # k           : The Boltzmann's constant in eV.K^-1 (From Gilloly et al. 2001)
  # r0          : A scaling constant at T0 in mgC^3/4.h^-1
  # rb          : The allometric coefficient of the respiration/mass relationship
  # k0          : A scaling constant at T0 in l.h^–1.mgC^-3/4  
  # h0          : A scaling constant at T0 in h.mgCfood^-1.mgC^-3/4
  # nb_ind      : Number of individual simulated 

  
  #----------------------#
  #   Model parameters   #
  #----------------------#
  
  # Open data files 
  
  data  <- read.csv("Forcing_Riki_2009-2010_8h.csv")                  # data forcing 

  # Read temperature and food concentration data
  tp_surf <- rep(data$Tmix, 5)
  tp_deep <- rep(data$Trasch, 5)
  
  phyto  <- rep(data$PhyC, 5)
  zoo    <- rep(data$ZC, 5)
  
  food   <- rbind(phyto, zoo) * 1e-3 # mg C L^-1
  
  #tp     <- data$temp 
  #diat   <- data$diat
  #flag   <- data$flag
  #micro  <- data$micro
  #meso   <- data$meso
  
  tsteps <- length(tp_surf)
  
  #C2N    <- c(6.625, 6.625, 7, 7) # carbon to nitrogen ratio for each category
  #food   <- rbind( diat*C2N[1] + flag*C2N[2] + micro*C2N[3] + meso*C2N[4], meso*0) * 14 * 1e-3 # mg C L^-1
  
  
  dt     <- 8        # time step of the model in h
  days   <- tsteps*dt/24   
  start  <- 1         #1*24/dt   # start at the end of June
  end    <- min( tsteps, start + days*24/dt )
  time   <- seq( start, end) # total number of time steps
  
  # shift time series to start 1st of August
  #tp     <- c( tp[(213*24/dt):length(tp)], tp[1:(213*24/dt-1)] )
  #food   <- cbind( food[,(213*24/dt):dim(food)[2]], food[,1:(213*24/dt-1)] )
  
  # Number of individuals simulated
  nb_ind <- 10
  
  # Parameters of allometric relationship between body mass and length 
  # M(mgC) = aw L(mm)^bw from actual experimental results (Cabrol's)
  aw     <- 0.0025 
  bw     <- 2.98        
  
  # Arrhenius parameters 
  T0     <- 273         # constant for both species
  k      <- 8.62e-5     # constant for both species
  
  # Parameters of reproduction costs from Laurie-Emma experiments
  a_rep  <- exp(-1.4552) 
  b_rep  <- 2.197         
  cw_egg <- 1.2e-3        
  
  # Developement time parameters
  a_molt  <- 15.3   # T. raschii data from T. inermis
  b_molt  <- -0.91  # T. raschii data from T. inermis 
  
  # Proportion of mass lost during moulting (From Sameoto 1976)
  w_molt  <- 0.005 
   
  # Threshhold of food concentration to caracterised a bloom
  food_bloom <- 0.2
  
  
  ### PARAMETERS FOR SENSIBILITY ANALYSE
  
  if (sum(unlist(pars_sens))==0){
    
  # Respiration parameters from actual experimental results (Angelique's)
  r0      <- exp(-7.59383) 
  rb      <- 0.958
    
  # Ingestion parameters 
  A       <- 0.7
  
  # Gonad parameters 
  Q  <- 0.4
    
  ### CALIBRATED PARAMETERS
  
  if(sum(unlist(pars))==0){
    
    # Ingestion parameters
    k0      <- c(0.4897, 0.1204) # T. raschii
    h0      <- 269.465           # T. raschii
    
    # Activation energy the same for all temperature-dependant processes 
    # from actual experimental results (Angelique's)
    ei      <- 0.2 
    er      <- 0.5
    
  } else if (is.list(pars)) {

    # Ingestion parameters
    k0 <- pars$k0
    h0 <- pars$h0
    
    # Activation energy 
    ei      <- pars$ei
    er      <- pars$er

  } else {

    # Ingestion parameters
    k0 <- pars[1:2]
    h0 <- pars[3]
    # Activation energy 
    ei <- pars[4]
    er <- pars[5]
  }
  
} else {
  A  <- pars_sens[1]
  k0 <- pars_sens[2:3]
  h0 <- pars_sens[4]
  r0 <- pars_sens[5]
  rb <- pars_sens[6]
  Q  <- pars_sens[7]
  ei <- pars_sens[8]
  er <- 0.5
}
  
  ### END CALIBRATED PARAMETERS
 
  #----------------------#
  #   Model functions    #
  #----------------------#
  
  
  # Photoperiod function (for feeding)
  sun <- function(lat=48.5, lon=-63, tzone=-5, day=180) {
    
    # angle conversion
    r2d    <- 180/pi
    d2r    <- pi/180
    
    # 'simple' equation from NOAA
    g      <- 2*pi/365 * day
    
    eqtime <- 229.18 * (  0.000075 + 0.001868*cos(g) - 0.032077*sin(g) - 0.014615*cos(2*g) - 0.040849*sin(2*g) )
    
    decl   <-  0.006918 - 0.399912*cos(g) + 0.070257*sin(g) - 0.006758*cos(2*g) + 0.000907*sin(2*g) - 0.002697*cos(3*g) + 0.001480*sin(3*g)
    
    ha     <- acos( cos(90.833*d2r) / (cos(lat*d2r)*cos(decl)) - tan(lat*d2r)*tan(decl) )
    
    SunR   <- ( trunc( 720+4.*(lon-ha*r2d)-eqtime )-tzone*60 )%%1440 / 60
    
    SunS   <- ( trunc( 720+4.*(lon+ha*r2d)-eqtime )-tzone*60 )%%1440 / 60
    
    night  <- 1 - (SunS - SunR) / 24
    
    return(list(sunrise=SunR, sunset=SunS, night=night))
  }
  
  
  # Arrhenius function of temperature, unitless
  arrhenius <- function(t, ei=0.2, T0=273, k=8.62e-5) { 
    
    at <- exp( ei * t / (k*(t+T0)*T0) )
    
    return(at)
  }
  
  
  # Development time (intermoult period). Default parameters from Sameoto 1976 Journal of the Fisheries Board of Canada (CJFAS) 33:2568-2576 (GSL)
  # Cuzin & Buchholz 1999 doi:10.1007/s002270050466 (Mediterranean) provide a_molt=19.1 & b_molt=-0.843
  devel <- function(t, a_molt, b_molt) {
    
    d <- a_molt + b_molt * t
    
    return(d)
  }
  
  
  # Ingestion : Holling II * allometry(mass) * Arrhenius(T)
  # expressed in mgC ind^-1 day^-1
  ingest <- function(w, f, t, dt=1, k0=k0, h0=h0, ei=0.2) {
    
    w[w==0] <- 1e-6
    
    a     <- outer( k0, w^(3/4) ) * arrhenius(t,ei) # parameter used to compute effective encounter rate. Increases with temperature
    
    b     <- h0 * w^(-3/4) * arrhenius(t,-ei)       # Handling time in h. Decreases with temperature.
    b     <- matrix(b, nrow=length(f), ncol=length(b), byrow=TRUE)

    alpha <- a * f / ( 1 + a * b * f )              # varying effective encounter rate in L h^-1(for switching function)
    
    ing   <- matrix(0, nrow=length(f), ncol=length(w))
    
    for (i in 1:length(f)) {
      #ing[i,]  <- a[i,]     * f[i]^2 / ( 1 + colSums( a * b * f^2 ) )   * dt # Tian 2006 Eq. 60: constant preference
      ing[i,]  <- alpha[i,] * f[i]   / ( 1 + colSums( alpha * b * f ) ) * dt # Tian 2006 Eq. 67: preference follows abundance
    }
    
    return( list(ing=ing, a=a, b=b, alpha=alpha) )
  }
  
  
  # Respiration: basal metabolism + activity (swimming) metabolism 
  # but NOT specific dynamic action (SDA) as it is included in assimilation coeffcient
  # expressed in mgC.ind^-1.day^-1
  # Data from Angelique (TODO)
  respiration <- function (w, temp, dt=1, r0, rb, er=0.5) {
    
    resp  <- r0 * w^rb * arrhenius(temp,er) * dt
    return(resp)
  }
  
  
  # Reproduction cost
  
  egg_prod <- function(lg_ind, a_rep, b_rep) {
    
    nb_egg <- a_rep*(lg_ind^(b_rep))
    return(nb_egg)
    
  }
  
  #----------------------#
  # RUNNING the model:   #
  #----------------------#
  
  # Initialization
 
  # Nighttime daily proportion
  night <- sun( day=(time*dt/24) )$night

  # actual temperature experienced by migrating krill
  tp <- night * tp_surf + ( 1 - night ) * tp_deep
    
  # Individuals' length (mm)
  lg_ind <- matrix(0, nrow=end, ncol=nb_ind) # Individuals' length in mm
  if (sum(unlist(pars_sens))==0){
    lg_ind[start,] <- round( runif(nb_ind, min(11, na.rm=TRUE), max(27, na.rm=TRUE)), digits=2) 
  } else {
    lg_ind[start,]     <- lg_ini
  }
  
  #lg_ind[start,] <- round( rnorm(nb_ind, mean=19, sd=2.6), digits=2) 
  
  # Individual carbon mass (mg C) according to allometric relationship
  cw_ind         <- matrix(0, nrow=end, ncol=nb_ind) # Individuals' dry weight in mgC
  cw_ind[start,] <- aw * lg_ind[start,]^bw
  
  # Proportion of individual carbon coming from mesozooplankton (copepods)
  cp_mez  <- 0.26 + 0*cw_ind # T. raschii
  
  ing_ind <- matrix(0, nrow=end, ncol=nb_ind)
  ing_f   <- array(0,  dim=c(dim(food)[1], nb_ind, end) ); a <- ing_f; b <- ing_f; alpha <- ing_f
  
  res_ind <- matrix(0, nrow=end, ncol=nb_ind)
  
  # Development fraction (from 0 to 1 = moulting)
  fr_molt    <- rep(0,nb_ind)
  i_ready    <- seq(1,nb_ind)
  i_molt     <- c()
  
  # Reproduction fraction 
  g_cum     <- rep(0,nb_ind)
  i_repro   <- rep(T,nb_ind)
  
  egg_cum   <- matrix(0, nrow=end, ncol=nb_ind)
  egg_batch <- matrix(0, nrow=end, ncol=nb_ind)
  
  gonad     <- 0.05 * cw_ind
  max_gonad <- 0.1  * cw_ind
  pgonad    <- rep(0.05,nb_ind)
  mature    <- matrix(FALSE, nrow=end, ncol=nb_ind)
  
  # Individuals' length when moulting
  lg_molt  <- lg_ind[start,]
  nb_moult <- rep(0, nb_ind)
  
  # Loop on time
  for (t in time[2:length(time)]) {

    # growth
    ingestion   <- ingest(cw_ind[t-1,], food[,t-1], tp[t-1], dt, k0, h0, ei)
    ing         <- ingestion$ing
    ing_ind[t,] <- colSums(ing)

    ing_f[,,t]  <- ingestion$ing # save diagnostics data
    a[,,t]      <- ingestion$a
    b[,,t]      <- ingestion$b
    alpha[,,t]  <- ingestion$alpha
    
    res_ind[t,] <- respiration(cw_ind[t-1,], tp[t-1], dt, r0, rb, er)
    
    grow        <- (A * night[t] * ing_ind[t,]) - res_ind[t,]
    i           <- cw_ind[t-1,] + grow < 0
    grow[i]     <- -cw_ind[t-1,i]
    
    g_cum       <- g_cum + grow
    
    pgonad      <- gonad[t-1,] / cw_ind[t-1,]
    cw_ind[t,]  <- cw_ind[t-1,] + grow
    
    max_gonad[t,] <- 0.1 * cw_ind[t,]
    
    # gonad growth
    gonad[t,]  <- gonad[t-1,]
    
    mature[t,] <- mature[t-1,]
    
    i          <- grow <= 0
    gonad[t,i] <- pgonad[i] * cw_ind[t,i]
    
    i          <- !mature[t,] & grow>0 & gonad[t,]<max_gonad[t,]
    gonad[t,i] <- gonad[t-1,i] + Q*grow[i]
    
    
    # proportion of C from mesozooplankton (copepods)
    i           <- grow>0
    cp_mez[t,]  <-   cp_mez[t-1,]
    cp_mez[t,i] <- ( cp_mez[t-1,i] * cw_ind[t-1,i] + grow[i] * ing[2,i]/ing_ind[t,i] ) / cw_ind[t,i]
    
    # moulting is irreversibly engaged after 40% of the IMP has passed
    # the new individuals' length is decided at that time according to the allometric relationship
    
    lg_ind[t,]  <- lg_ind[t-1,]
    
    fr_molt     <- fr_molt + dt / (devel(tp[t-1], a_molt, b_molt)*24)
    
    # check whether 40% of the IMP has passed
    # if so compute new lengths and move individual indices around 
    i <- fr_molt[i_ready]-0.4 >=0 
    j <- i_ready[i]
    
    lg_molt[j]  <- (cw_ind[t,j]/aw)^(1/bw)
    
    i_ready     <- i_ready[ which(!i) ]
    
    i_molt      <- c(i_molt,j)
    
    # check whether 100% of the IMP has passed
    # if so compute new development fractions, new masses, record new lengths and move individual indices around 
    i <- fr_molt[i_molt]-1 >=0
    j <- i_molt[i]
    
    # Moulting
    if ( length(j) != 0 ) {
      
      for (k in j) {
        
        fr_molt[k]  <- fr_molt[k]-1
        nb_moult[k] <- nb_moult[k]+1
        
        egg <- 0 * nb_moult[k]
        reproduction <- egg
        
        # Maturation decision
        if ( !mature[t,k] & gonad[t,k] >= max_gonad[t,k] & grow[k] >= 0) {
          mature[t,k] <- TRUE
          nb_moult[k] <- 0
        }
        
        # Egg production decision
        if ( mature[t,k] & nb_moult[k]%%1==0 & food[1,t]>food_bloom ) {
          
          egg <- egg_prod(lg_ind[t,k], a_rep, b_rep) 
          
          egg_batch[t,k] <- egg
          reproduction   <- cw_egg * egg

          cw_ind[t,k]    <- max(0,  cw_ind[t,k] - reproduction) # cost of reproduction on body mass
          
          gonad[t,k]     <- max(0, gonad[t-1,k] - reproduction) # cost of reproduction on gonad
          
          if( reproduction > gonad[t-1,k] ) { # if gonad depleted -> immature
            mature[t,k] <- FALSE
          }

        }
        
        # Stop maturation decision
        if ( mature[t,k] & ( food[1,t]<food_bloom )) { # | ( nb_moult[k]%%2==0 & grow[k]<0 ) ) ) {
          mature[t,k] <- FALSE
        }
        
        # Accounting of exuvia loss
        if( cw_ind[t,k]!=0 ) {
          pgonad <- gonad[t,k] / cw_ind[t,k]
        } else {
          pgonad <- 0
        }

        cw_ind[t,k] <- cw_ind[t,k] * (1-w_molt) # cost of exuvia loss
        
        gonad[t,k]  <- pgonad * cw_ind[t,k] # adjust gonad mass to keep proportion constant

        g_cum[k]    <- 0 
        
        lg_ind[t,k] <- lg_molt[k]
        
        i_molt      <- i_molt[which(!i)]
        
        i_ready     <- c(i_ready,k)
        
        reproduction <- 0
      }
    }

  }
  
  result <- list(food=food, tp_surf=tp_surf, tp_deep=tp_deep, lg_ind=lg_ind, cw_ind=cw_ind, gonad=gonad, max_gonad=max_gonad, mature=mature, cp_mez=cp_mez, ing_ind=ing_ind, ing_f=ing_f, a=a, b=b, alpha=alpha, res_ind=res_ind, night=night, egg_batch=egg_batch)

  #write.csv(as.data.frame(result), file="results_Traschii.csv")

  return(result)
}

