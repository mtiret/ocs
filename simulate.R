###########################################################
## libraries ##############################################
###########################################################
require( MASS      , quietly = T, warn.conflicts = F ) # load before dplyr.
require( dplyr     , quietly = T, warn.conflicts = F )
require( data.table, quietly = T, warn.conflicts = F ) # fread, fwrite.

###########################################################
## script parameters ######################################
###########################################################
## number of individuals.
npar = 1009
noff = npar

seed = sample( 1:2^15, 1 )
set.seed( seed )

## reading arguments.
args = commandArgs(trailingOnly = TRUE)

## number of generations.
ngeneration = eval( parse( text = args[1] ) )

## equation parameters.
alpha = eval( parse( text = args[2] ) )
beta  = eval( parse( text = args[3] ) )
gamma = eval( parse( text = args[4] ) )

## massal selection option.
# mass_selection = masstobereplaced
# ndpar = 100

## avoiding scientific notations.
options(scipen=999)

# reading allele effects.
qtn_option = args[5]
qtn.effect = fread( qtn_option )$V3

# heritability to randomise GEBV.
h2 = 0.5134 #eval( parse( text = args[6] ) )

# evaluation boolean.
eval = (h2 != 1)

mate_allocation = eval( parse( text = args[6] ) )

shuffle = eval( parse( text = args[7] ) )


###########################################################
## impl functions #########################################
###########################################################
# constructing individual wise genotypes from biallelic data.
# -1, 0, 1
get_ind_genotype = function( file_name )
{
  .genotype_data = fread( file_name ) %>% as.matrix()
  
  ind_genotype =
    .genotype_data[, which(1:ncol(.genotype_data) %% 2 == 0) ] +
    .genotype_data[, which(1:ncol(.genotype_data) %% 2 == 1) ] - 1
  
  return( ind_genotype )
}

# computing gebv.
compute_gbv = function( .genotype_data, .h2 )
{
  # computing the true breeding values (after removing label and positions).
  .TBV = ( t(.genotype_data) + 1 ) %*% qtn.effect
  
  # randomising the GEBV (phenotypic score)
  .PS = .TBV + rnorm( length(.TBV), 0, sqrt((1-.h2)/.h2 * var(.TBV)) )
  
  return( data.frame( TBV = .TBV, PS = .PS ) )
}

# computing grm. (from -1 0 1).
compute_grm = function( .genotype_data, .normalise = T )
{
  grm_data = .genotype_data %>% as.matrix()
  
  if( !.normalise || (beta == 0 && gamma == 0) )
  {
    # normalising by the number of snps (for numeric stability).
    grm_data = (t(grm_data) %*% grm_data) / nrow(grm_data)
    
  } else {
    
    # only couples of 0/1 or 0/-1 have negative values. 
    Q = nrow(grm_data) * 0.5 - 2 * t(abs(grm_data)-0.5) %*% (abs(grm_data)-0.5)
    
    # couples of 0x0.
    W = (grm_data - 1) * (grm_data + 1)
    
    ## computing offsprings' Jaccard matrix. ##################
    grm_joff = t(grm_data) %*% grm_data + gamma * t(W) %*% W + beta * Q
    
    ## projection on the set of PSD matrices. #################
    # @dev : using eigenvalue decomposition. 
    .eigen = eigen( grm_joff, symmetric = T )
    .eigen$values[ .eigen$values < 0 ] = 1e-6
    
    # normalising by the number of snps (for numeric stability).
    # orthogonal matrix because %grm_joff is symetric. 
    grm_data = .eigen$vectors %*% diag(.eigen$values) %*% t(.eigen$vectors) / nrow(.genotype_data)
  }
  
  return( grm_data )
}

# Newton method.
ocs_newton = function( .G, .Y, .alpha, 
                       cgcmin = 0, cgcmax = 1, gmin = 0, gmax = 1 )
{
  # parameters
  a   = 0.01
  b   = 0.8
  eps = 1e-7
  lb  = 2
  mu  = 10
  
  N = ncol(.G)
  
  # objective function.
  f = function( .ct )
  {
    if( !all(.ct > 0) || !all(.ct < 0.5)  )
      return(NA)
    else
      return( 
        .alpha * ( t(.ct) %*% .G %*% .ct - cgcmin)/(cgcmax - cgcmin) +
          - ( 1 - .alpha ) * (t(.Y) %*% .ct - gmin)/(gmax - gmin)+ 
          - 1/lb * sum(log(.ct)) + 
          - 1/lb * sum(log(0.5-.ct)) 
      )
  }
  
  # taking care of the equality constraint.
  ct = rep( 1 / N, N )
  
  while( N/lb > eps )
  {
    # centering step.
    lambda2 = 1
    while( lambda2/2 > eps )
    {
      # derivatives.
      H = .alpha * .G * 2 / (cgcmax - cgcmin)
      diag(H) = diag(H) + 1/lb*( 1/c(ct)^2 + 1/c(0.5-ct)^2 )
      
      
      Df = .alpha * .G %*% ct * 2 /(cgcmax - cgcmin) + 
        - (1-.alpha) * .Y /(gmax - gmin) + 
        - 1/lb/ct + 1/lb/(0.5-ct)
      
      # compute newton step (KKT system)
      # Hm1 = solve(H)
      # nu = - sum( Hm1 %*% Df )/( sum(Hm1) )
      # Dx = Hm1 %*% ( - nu - Df )
      # @dev : backward/forward substitution.
      Hchol   = H %>% chol() # limiting operation.
      sum_hm1 = backsolve( Hchol,
                           forwardsolve(t(Hchol), rep(1,npar)) ) %>% sum()
      sum_hm1df = backsolve( Hchol,
                             forwardsolve(t(Hchol), Df) ) %>% sum()

      nu = - sum_hm1df / sum_hm1

      Dx = backsolve( Hchol,
                      forwardsolve(t(Hchol), - nu - Df) )
      
      # compute newton decrement.
      lambda2 = (Hchol %*% Dx) %>% crossprod()
      
      # line search
      t = 1
      while( (!all(ct+t*Dx > 0) || !all(ct+t*Dx < 0.5)) ||
             ( f(ct+t*Dx) - f(ct) - a * t * sum( Df * Dx ) > 0 &&
               abs(f(ct+t*Dx) - f(ct) - a * t * sum( Df * Dx )) > eps ) )
        t = b * t
      
      # update
      ct = ct + t * Dx
    }
    
    # update.
    lb = mu * lb
  }
  
  return( ct )
}

discrete_contrib = function( .ct )
{
  # rounding the contributions. 
  sum_ct = sum(.ct)
  merit = .ct - floor(.ct)
  .ct = floor(.ct)
  while( sum(.ct) != round(sum_ct) )
  {
    # deterministic.
    index = which.max( merit )
    merit[index] = 0
    
    # probabilistic.
    # a = which( rmultinom(1, 1, ct.counter) > 0 ) 
    # ct.counter[a] = max( ct.counter[a] - 1, 0 )
    
    .ct[index] = .ct[index] + 1
  }
  
  return(.ct)
}

write_pedigree = function( .ct )
{
  # constructing pedigree.
  index    = which( .ct > 0 )
  parents  = rep( index, .ct[index] ) %>% sample
  pedigree = matrix( parents, ncol = 2 )
  
  # correcting selfing.
  selfing = which( pedigree[,1] == pedigree[,2] )
  for( i in selfing )
  {
    if( pedigree[i,1] == pedigree[i,2] )
    {
      for( j in 1:nrow(pedigree) )
      {
        if( pedigree[j,1] != pedigree[i,1] && 
            pedigree[j,2] != pedigree[i,1] )
        {
          tmp = pedigree[j,1]
          pedigree[j,1] = pedigree[i,1]
          pedigree[i,1] = tmp 
          break
        }
      }
    }
  }
  
  return( pedigree)
}

discret_mc = function( .C )
{
  # rounding the contributions. 
  sum_ct = sum(.C)
  merit = .C - floor(.C)
  .C = floor(.C)
  while( sum(.C) != round(sum_ct) )
  {
    # deterministic.
    index = which.max( merit )
    merit[index] = 0
    
    # probabilistic.
    # a = which( rmultinom(1, 1, ct.counter) > 0 ) 
    # ct.counter[a] = max( ct.counter[a] - 1, 0 )
    
    .C[index] = .C[index] + 1
  }
  
  .C
}

write_mpedigree = function( .C )
{
  pedigree = matrix( NA, nrow = sum(.C)/2, ncol = 2 )
  k = 1
  for( i in 2:nrow(.C) ) 
    for( j in 1:(i-1) )
    {
      if( .C[i,j] > 0 )
      {
        pedigree[k:(k+.C[i,j]-1),] = 
          matrix( rep( c(i,j), .C[i,j] ), ncol = 2, byrow = T )
        k = k + .C[i,j]
      }
    }
  
  pedigree
}

mate_newton = function( .ct, Gij, A )
{
  # objective function.
  f = function( .cij )
  {
    if( !all(.cij >= 0) || !all(.cij <= 1) )
      return(NA)
    else
      return(c(
        crossprod(.cij, Gij) - (sum(log(.cij)) + sum(log(1-.cij)))/lb
      ))
  }
  
  N = npar 
  
  # Newton method. --------------------------------
  # parameters
  a   = 0.01
  b   = 0.8
  eps = 1e-7
  lb  = 2
  mu  = 10
  
  while( N/lb > eps )
  {
    # centering step.
    lambda2 = 1
    while( lambda2/2 > eps )
    {
      # derivatives.
      # H = diag( c(1/.ct^2 + 1/(.5 - .ct)^2) )/lb
      Hdiag = c(1/.ct^2 + 1/(1 - .ct)^2)/lb
      
      Df = t(t(Gij - 1/lb/.ct + 1/lb/(1-.ct)))
      
      # # @dev : should be backward/forward substitution instead.
      # Hchol.up.inv = H %>% chol() %>% solve()
      # Hm1 = Hchol.up.inv %*% t(Hchol.up.inv)
      
      Hdiagm1 = 1/Hdiag
      
      # compute newton step (KKT system)
      sq = Hdiagm1 %>% sqrt
      
      ## @dev: limiting step. ----
      AHAt = tcrossprod( sweep(A, 2, sq, `*`) ) # A %*% diag(sq)
      
      # @dev: ginv: handle singular matrix. 
      nu = - ginv( AHAt ) %*% A %*% (Hdiagm1*Df)
      # nu = - solve( A %*% Hm1 %*% t(A) ) %*% (A %*% Hm1 %*% Df)
      ## -------------------------
      
      Dx = Hdiagm1 * ( - crossprod(A,nu) - Df )
      # Dx = Hm1 %*% ( - t(A) %*% nu - Df )
      
      # compute newton decrement.
      lambda2 = crossprod( Dx, Hdiag * Dx )
      # lambda2 = t(Dx) %*% H %*% Dx
      
      # line search
      t = 1
      while( (!all(.ct+t*Dx > 0) || !all(.ct+t*Dx < 1)) ||
             ( f(.ct+t*Dx) >= f(.ct) + a * t * t(Df) %*% Dx &&
               abs(f(.ct+t*Dx) - f(.ct) - a * t * t(Df) %*% Dx) >= eps ) )
        t = b * t
      
      # update
      .ct = .ct + t * Dx
    }
    
    # update.
    lb = mu * lb
  }
  
  return( .ct )
}

random_pedigree = function( .ct, .G, .nburn )
{
  best.ped = write_pedigree(.ct)
  
  for( i in 1:.nburn )
  {
    new = write_pedigree(.ct)
    if( sum(.G[new]) < sum(.G[best.ped]) )
      best.ped = new
  }
  
  best.ped
}

sa_pedigree = function(.ct, .G, kmax, kept.proportion = .9, restart = 10 )
{
  best.ped = write_pedigree(.ct)
  
  for( i in 1:restart )
  {
    ped = write_pedigree(.ct)
    for( k in 0:(kmax-1) )
    {
      # temperature.
      temperature = (k+1)/kmax
      
      # neighbours.
      index = order(.G[ped]) < floor(kept.proportion * length(.ct))
      rm.ped = ped[!index, ]
      rm.dct = rm.ped %>% factor( levels = 1:length(.ct) ) %>% table
      
      new = rbind(
        ped[index,],
        write_pedigree(rm.dct)
      )
      
      # replacing.
      enew = sum(.G[new])
      eold = sum(.G[ped])
      if( enew < eold || exp(-(enew-eold)/temperature) >= runif(1) )
        ped = new
    }
    
    if( sum(.G[ped]) < sum(.G[best.ped]) )
      best.ped = ped
  }
  
  best.ped
}

date = function()
  paste0( "[",format(Sys.time(), "%Y/%m/%d-%H:%M:%S"), "]")


#####################################################################
## Simulation #######################################################
#####################################################################
# initializing. #################################
if( shuffle > 0 )
{
  system( "cp in.gen.ran in.gen.cp" )
} else
  system( "cp in.gen in.gen.cp" )
out.p = tibble( cgc.p = NA, cy.p = NA )
out.r = tibble()

Xprev = NULL
BVprev = NULL

# mate allocation variables. 
# pos = combn(npar,2)
# 
# # @dev: constant population size. 
# A = matrix( 0, nrow = npar, ncol = npar*(npar-1)/2 )
# for( i in 1:ncol(pos) ) A[pos[,i],i] = 1

cat( paste0( "-- seed: ", seed, " --\n" ) )

for( i in 0:ngeneration )
{
  # inputting. ##################################
  cat( paste0( "-- generation n.", i, " -----------------------------\n") )
  ind_genotype = get_ind_genotype( "in.gen.cp" )
  BV = compute_gbv( ind_genotype, h2 )
  G  = compute_grm( ind_genotype )
  mG = compute_grm( ind_genotype, .normalise = F )
  
  # constructing the contributions. #############
  if( i < ngeneration )
  {
    # estimating the EBV. ###################################
    ## Estimating BV (for GS; otherwise it is phenotypic selection).
    ## @opt: using matrix inversion lemma. 
    ## Searle, 2006: correcting bias: lambda = 1/h2 - 1, instead of 1/h2. 
    ## RR: 
    ##    EBV = t(X) * solve( X * t(X) + lambda * diag(nrow(X)) ) * X * Y
    ##        = lambda * tXX * Y - lambda^2 * tXX * solve(diag(ncol(X)) + lambda * tXX) * tXX * Y
    ##
    ## RR-BLUP (Piepho 2012):
    ##    Z = rbind(1,X)
    ##    D = rbind( 0, cbind( 0, diag(nrow(X))) )
    ##    EBV = t(Z) * solve( ZtZ + 1/h2 * D ) * Z * Y
    ##        
    ## cor(RR, RR-blup) ~ 1 if no fixed effect -> using RR, much faster.  
    if( eval )
    {
      cat( paste0(date(), ": estimating BV...\n") )
      X = cbind( Xprev, ind_genotype %>% as.matrix())
      BVnew = rbind( BVprev, BV )
      
      lambda = (1/h2-1)^-1
      tXX    = crossprod(X)
      tXXY   = tXX %*% BVnew$PS 
      
      # RR.
      EBV = lambda * tXXY - 
        lambda^2 * tXX %*% solve( diag(ncol(X)) + lambda * tXX, tXXY )
    
      Xprev = X
      BVprev = BVnew
      EBV = EBV[(length(EBV)-npar+1):length(EBV)]
      
      cat( paste0(date(), ": estimating done.\n" ) )
      
    } else
      EBV = BV$TBV
    
    # constructing the contribution vector. #################
    # only selection, max inbreeding.
    ct0 = matrix( 0, nrow = npar, ncol = 1 )
    ct0[order(EBV, decreasing = T)[1:2]] = 0.5

    # no selection, min inbreeding.
    ct1 = ocs_newton(G, EBV, 1)

    # computing boundaries.
    cgcmin = t(ct1) %*% G %*% ct1 %>% c
    cgcmax = t(ct0) %*% G %*% ct0 %>% c
    gmin = t(ct1) %*% EBV %>% c
    gmax = t(ct0) %*% EBV %>% c
    
    # mass selection.
    # ct = array( 0, length(EBV) ) 
    # ct[ EBV %>% order( decreasing = T ) %>% head(ndpar) ] = noff / ndpar 
    
    # OCS. 
    cat( paste0(date(), ": OCS...\n" ))
    ct = 2 * noff * 
      ocs_newton( G, EBV, alpha, cgcmin, cgcmax, gmin, gmax )
      # ocs_newton( G, EBV, alpha )
    ct = discrete_contrib(ct)
    cat( paste0(date(), ": OCS done.\n") )
    
    # mate allocation. --------------------------
    if( mate_allocation > 0 )
    {
      cat( paste0(date(), ": mate allocation...\n" ) )
      
      index = which( ct > 0 )
      
      if( length(index) > 2 )
      {
        pos = combn( length(index), 2 )
        A = matrix( 0, nrow = length(index), ncol = ncol(pos) )
        for( j in 1:ncol(pos) ) A[pos[,j],j] = 1
        
        # Yij = apply( pos, 2, function(i) sum(EBV[index[i]]))
        Gij = G[index,index]
        Gij = Gij[lower.tri(Gij)]
        
        ped.old = write_pedigree(ct)
        C = matrix( 0, nrow = length(index), ncol = length(index) )
        for( j in 1:nrow(ped.old) )
          C[which( index == ped.old[j,1]), which( index == ped.old[j,2])] =
          C[which( index == ped.old[j,2]), which( index == ped.old[j,1])] = 
          C[which( index == ped.old[j,1]), 
            which( index == ped.old[j,2])] + 1
        
        Cij = C[lower.tri(C)] + 1e-7
        Cij = Cij / sum(Cij)
        
        Copt = mate_newton( Cij, Gij, A )
        C.opt = matrix(0, nrow = length(index), ncol = length(index) )
        C.opt[lower.tri(C.opt)] = Copt
        C.opt = (C.opt + t(C.opt))
        C.opt = C.opt/sum(C.opt) * 2 * noff
        
        ped.new = write_mpedigree( discret_mc(C.opt) )
        ped = cbind( index[ped.new[,1]], index[ped.new[,2]] )
        
        mate.gain = sum( G[ped.old] ) - sum( G[ped] ) 
        mate.gain.t = sum( mG[ped.old] ) - sum( mG[ped] ) 
      } else {
        ped = write_pedigree(ct)
        
        mate.gain   = 0
        mate.gain.t = 0
      }
      
      
      fwrite( data.frame( ped ), "in.ped",
              col.names = F, sep = " " )
      
      cat( paste0(date(), ": mate allocation done.\n" ) )
      
    } else {
      
      ped = write_pedigree( ct )
      
      fwrite( data.frame( ped ), "in.ped",
              col.names = F, sep = " " )
      
      mate.gain   = 0
      mate.gain.t = 0
    }
    
    
    # storing. ----------------------------------
    nparent = sum(ct > 0)
  }
  
  # outputting. #################################
  # prediction.
  if( i < ngeneration )
  {
    out.p = 
      bind_rows(
        out.p, 
        tibble( cgc.p = t(ct) %*% G %*% ct / sum(ct)^2, 
                cy.p  = t(EBV) %*% ct / sum(ct))
      )
  }
  
  # characterisation. 
  out.freq = (ind_genotype %>% rowMeans() + 1)/2 
  out.r =
    bind_rows(
      out.r,
      tibble(
        cgc.r = mean(G),
        cy.r  = mean(BV$TBV),
        
        mate.gain = mate.gain, 
        mate.gain.t = mate.gain.t, 
        
        # other diversity indices.
        cgc.t = mean(mG), # non normalised coancestry. 
        nparent,    # number of parent
        he = 2 * mean(out.freq * (1-out.freq)), # heterozygosity
        ho = 1 - he, # homozygosity
        
        vtot = var(BV$TBV),
        vgen = sum(out.freq * (1 - out.freq) * 2 * qtn.effect^2),
        
        lal.b = 
          sum( (qtn.effect > 0) * (out.freq == 0) ),
        lal.b.effect = 
          sum( (qtn.effect) * (qtn.effect > 0) * (out.freq == 0)),
        lal.d = 
          sum((qtn.effect < 0) * (out.freq == 1)),
        lal.d.effect = 
          sum( (qtn.effect) * (qtn.effect < 0) * (out.freq == 1) ),
        gal.b = 
          sum((qtn.effect > 0) * (out.freq == 1)), 
        gal.b.effect = 
          sum( (qtn.effect) * (qtn.effect > 0) * (out.freq == 1) ),
        gal.d = 
          sum( (qtn.effect < 0) * (out.freq == 0) ),
        gal.d.effect = 
          sum( (qtn.effect) * (qtn.effect < 0) * (out.freq == 0) ),
      )
    )
  
  # simulating. #################################
  if( i < ngeneration )
  {
    system( "./meiosis -m in.map.small -g in.gen.cp -p in.ped > out.gen" )
    system( "cp out.gen in.gen.cp" )
  }
}

cat( paste0( date(), ": recording done.\n" ) )

## final outputs. #########################################
bind_cols(
  out.p, out.r
) %>% 
  mutate( alpha = alpha, 
          gen   = seq_along(cy.r) - 1, 
          beta  = beta,
          gamma = gamma,
          mate  = mate_allocation,
          shuffle = shuffle, 
          
          # inbreeding rate (Falconer 1996).
          ir = c(NA, diff(ho)/(1-ho[-length(ho)])) 
  ) %>% 
  select( alpha, beta, gamma, mate, shuffle, gen, everything() ) %>% 
  fwrite( "out.cgc" )




