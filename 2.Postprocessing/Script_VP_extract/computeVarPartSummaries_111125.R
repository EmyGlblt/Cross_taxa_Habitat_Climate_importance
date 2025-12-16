#.
#.              FULL POSTERIOR DISTRIBUTION - VARIANCE PARTITIONING
#.              ---------------------------------------------------
#.
#... The below code extend the variance partitioning R function available in github 
#... computeVariancePartitioning from the HMSC package (Tikhonov et al., 2022). Downloaded in June 2022.
#... This code focus only on the variance partition thus other computation have been 
#... disregarded and commented out of the script.
#... 
#... Code extended by Emy Guilbault, September 2022. Last changes September 2025.
#... Extension: calculate the pairwise Normalized, diagonal, marginal, partial and covariance partition from Schulz et al. (2025)
#... 
#..................................................................................



# implement main VP summarise function
#.........................................................................................................................

computeVarPartSummaries = function(hM, group=NULL, groupnames=NULL, start=1, na.ignore=FALSE, exclude = NULL, 
                                   conditional = NULL, cond=NULL, wb=NULL, marginal = NULL, Global = NULL)
{
  
  ## if we want to calculate the conditional VP .:
  # the conditions: conditional = NULL, cond=cond, wb=NULL 
  #  between and within groups:conditional = NULL, cond=cond, wb=TRUE 
  if(conditional == TRUE){
    VP = computeVariancePartitioningbwcond(hM, group=group, groupnames=groupnames, start=1, na.ignore=FALSE, 
                                           exclude = exclude, cond=cond, wb=wb, marginal= marginal) 
  }else{
    
    # Otherwise: conditional = NULL, cond=NULL, wb=NULL
    
    ## get the dimension information
    ny = hM$ny
    nc = hM$nc
    nt = hM$nt
    ns = hM$ns
    np = hM$np
    nr = hM$nr
    
    ## because dimension is important when removing a covariate such as effort to calculate the VP we need to keep track of the names and posiition of the linear predictor below:
    # important for the end and names allocation
    if(is.null(group)){gp.name = FALSE}else{gp.name=TRUE}
    
    
    ## This consider if the user has pre entered grouping for linear predictors and if not it creates grouping for each individual linear predictor
    gp = group # keep track on if we had group submitted or not
    gpname = groupnames
    if(is.null(group)){
      ## default: use terms
      if(nc > 1){
        if (is.null(hM$XFormula))
          stop("no XFormula: you must give 'group' and 'groupnames'")
        
          group = attr(hM$X, "assign")
          if (group[1] == 0) # assign (Intercept) to group
            group[1] <- 1
          
          groupnames = attr(terms(hM$XFormula, data = hM$XData), "term.labels")
          gpname = groupnames
          
      } else {
        group = c(1)
        groupnames = hM$covNames[1]
        gpname = groupnames
      }
      
    }
    
    ngroups = max(group)
    
    ## this is not of use in the new version*********
    R2T.Y = 0
    R2T.Beta = rep(0,nc)
    ## this is not of use in yhe new version*********
    
    
    ## info over MCMC iterations
    postList = poolMcmcChains(hM$postList, start=start)
    
    nsamp = length(postList)#  # of mcmc samples
    
    ## Initialize arrays used later
    
    ## we initiate the normalized, marginal, diagonal, partition finalized VP matrix for each species, mcmc samples (Schulz et al., 2025) and covariance partition
    if(is.null(exclude)){
      Vnorm = array(NA, dim = c(1, ngroups+nr, ns, nsamp))
      Vdiag = array(NA, dim = c(1, ngroups+nr, ns, nsamp))
      Cnorm = array(NA, dim = c(ngroups+nr, ngroups+nr, ns, nsamp))   
      
      if(Global==TRUE){
        Vmarg = array(NA, dim = c(1, ngroups+nr, ns, nsamp))
        Vpart = array(NA, dim = c(1, ngroups+nr, ns, nsamp))
        
      }else{
        Vmarg = array(NA, dim = c(ngroups+nr, ngroups+nr, ns, nsamp))
        Vpart = array(NA, dim = c(ngroups+nr, ngroups+nr, ns, nsamp))
        
      }
      
    }else{
      Vnorm = array(NA, dim = c(1, ngroups+nr-length(exclude), ns, nsamp))   # if we exclude certain covariates like effort for the calculation
      Vdiag = array(NA, dim = c(1, ngroups+nr-length(exclude), ns, nsamp))
      Cnorm = array(NA, dim = c(ngroups+nr-length(exclude), ngroups+nr-length(exclude), ns, nsamp))  
      
      if(Global==TRUE){
        Vmarg = array(NA, dim = c(1, ngroups+nr-length(exclude), ns, nsamp))
        Vpart = array(NA, dim = c(1, ngroups+nr-length(exclude), ns, nsamp))
        
      }else{
        Vmarg = array(NA, dim = c(ngroups+nr-length(exclude), ngroups+nr-length(exclude), ns, nsamp))
        Vpart = array(NA, dim = c(ngroups+nr-length(exclude), ngroups+nr-length(exclude), ns, nsamp))
        
      }
    }
    
    ## prepare covariance - for fixed part only (This is technically already done from the new HMSC code line 120)
    ## not done before in Jarno's version
    #cMA = cov(hM$X)
    
    ## ........................................................run over MCMC chains 
    for (i in 1:nsamp){
      print(i)
      #i
      
      ### Preparation of the K matrix linear terms and covariance *******************************************************************************
      
      # 1 #++++++++++++++++++++++++++++++++++++++   FIXED part 
      
      ## estimates from the posterior sample
      Beta.i = postList[[i]]$Beta
      Lambdas.i = postList[[i]]$Lambda
      Etas.i = postList[[i]]$Eta
      
      ## calculation for all species
      
      # temporary names to use
      Beta.x = Beta.i
      
      LF.array = array(NA, c(dim(hM$X)[1], dim(hM$X)[2], ns))
      for (kk in 1:dim(hM$X)[2]){
        LF.temp = as.vector(hM$X[,kk])%*%t(as.vector(Beta.x[kk,]))
        LF.array[,kk,] = LF.temp
      }
      # 
      # for (jj in 1:dim(LF.array)[3]){
      #   cov.lf = cov(LF.array[,,jj])  ## if one random effect needs to add as.matrix
      #   colnames(cov.lf) = rownames(cov.lf) = colnames(hM$X)
      # }
      
      
      # 2 #++++++++++++++++++++++++++++++++++++++   RANDOM part 
      
      RF.array = array(NA, dim=c(dim(hM$Pi)[1], dim(hM$Pi)[2], ns))
      for (kk in 1:length(Lambdas.i)){ ## across random effect
        Pi.temp = matrix(0, nrow=dim(hM$Pi)[1], ncol=length(unique(hM$Pi[,kk])))
        for (jj in 1:dim(hM$Pi)[1]){  ## across sampling occasions
          Pi.temp[jj,hM$Pi[jj,kk]] = 1
        }
        RF.temp = Pi.temp%*%Etas.i[[kk]]%*%Lambdas.i[[kk]]
        RF.array[,kk,] = RF.temp
      }
      
      # for (jj in 1:dim(RF.array)[3]){
      #   cov.rf = cov(as.matrix(RF.array[,,jj]))  ## if one random effect needs to add as.matrix
      #   colnames(cov.rf) = rownames(cov.rf) = colnames(hM$Pi)
      # }
      
      
      
      # 3 #++++++++++++++++++++++++++++++++++++++   FIXED - RANDOM covar part (only exist when considering Update = TRUE)
      
      
      # Combined covariance partition matrix of all fixed effects and random effects
      dim(LF.array)
      dim(RF.array)
      
      SIGMA = array(NA, dim=c(dim(hM$X)[2]+dim(hM$Pi)[2], dim(hM$X)[2]+dim(hM$Pi)[2], ns))
      for (sp in 1:ns) {
        SIGMA[,,sp] = cov(cbind(LF.array[,,sp], RF.array[,,sp]))
      }
      dimnames(SIGMA) = list(element=c(colnames(hM$X), colnames(hM$Pi)), 
                             element=c(colnames(hM$X), colnames(hM$Pi)),
                             Species=hM$sp)
      
      # 4 # -------- Do we need grouping?
      # In any case we group the intercept so:
      ## adapted from Torsti's example
      
      linear_terms = c(colnames(hM$XScaled), colnames(hM$Pi))
      group.tmp = c(group, seq(from = max(group)+1, length.out= dim(hM$dfPi)[2]))
      groups <- split(linear_terms, group.tmp)
      names(groups) = c(groupnames, colnames(hM$Pi))
      
      B <- sapply(groups, function(x) { as.integer(linear_terms %in% x) })
      dimnames(B) <- list(linear_terms = linear_terms, groups = colnames(B))
      
      n_group = max(group.tmp)
      K_B <- apply(SIGMA, "Species", function(x, B) { t(B) %*% x %*% B }, B)
      dim(K_B) <- c(n_group, n_group, ns)
      dimnames(K_B) <- list(rows = c(groupnames, colnames(hM$Pi)), cols = c(groupnames, colnames(hM$Pi)), species = hM$sp)
      
      SIGMA.final = K_B
      
      ### Matrix K linear terms  ready *******************************************************************************
      ### Similar to Torsti's paper *******************************************************************************
      
      
      # Total and partition calculation ******************************************************************************************
      
      if(is.null(exclude)){
        for (sp in 1:ns) {
          #Var_sp_tot = Var_sp_tot  # or
          SIGMA_sp = SIGMA.final[,,sp]
          Var_sp_tot = sum(SIGMA.final[,,sp]) #### take the whole triangle
          
          if(marginal == TRUE){
            ## normalized variance
            Vnorm_sp = diag(SIGMA.final[,,sp])/Var_sp_tot  # 
            Vnorm[,,sp,i] = Vnorm_sp
            
            ## marginal variance
            if(Global == TRUE){
              Vmarg_sp = rowSums(SIGMA.final[,,sp])/Var_sp_tot  # same as Torsti
              Vmarg[,,sp,i] =  Vmarg_sp
              
            }else{
              Vmarg_sp = matrix(NA, n_group, n_group)
              groups = colnames(SIGMA_sp)
              Vmarg_sp_diag = rowSums(SIGMA.final[,,sp])/Var_sp_tot  # same as Torsti, will be the diag elements
              
              S.tmp = SIGMA.final[,,sp]
              
              Vmarg_sp_nondiag = lapply(sapply(1:nrow(S.tmp), function(i) {
                row <- S.tmp[i, ]
                diag_element <- row[i]
                non_diag_elements <- row[-i]
                diag_element + non_diag_elements
              }, simplify=F, USE.NAMES = F), function(x){x/Var_sp_tot})   # same but for individual groups
              
              diag(Vmarg_sp) = Vmarg_sp_diag
              
              delta <- col(Vmarg_sp) - row(Vmarg_sp) ## select the non diagonal elements
              Vmarg_sp[delta > 0 | delta < 0] <- do.call(cbind, Vmarg_sp_nondiag) ## fill it in the columns
              
              colnames(Vmarg_sp) = rownames(Vmarg_sp) = colnames()
              Vmarg[,,sp,i] = t(Vmarg_sp)  ## transpose to get by row
            }
            
            
            ## covariance
            SIGMA_sp = SIGMA.final[,,sp]
            SIGMA_sp[! upper.tri(SIGMA_sp)] = NA
            
            Cnorm_sp = SIGMA_sp/Var_sp_tot
            Cnorm[,,sp,i] = Cnorm_sp  # diag + cov
            
            ## partial variance
            ## this may change if we have done the grouping before / to think later about it
            groups = colnames(SIGMA_sp)
            B_groups <- sapply(groups, function(x) { as.integer(colnames(SIGMA_sp) %in% x) })
            dimnames(B_groups) <- list(linear_terms = colnames(SIGMA_sp), groups = colnames(SIGMA_sp))
            
            if(Global == TRUE){
              Vpart_sp = partial_variances(SIGMA.final[,,sp], B = B_groups)/Var_sp_tot
              Vpart[,,sp,i] = Vpart_sp
              
            }else{
              Vpart_sp = matrix(NA, n_group, n_group)
              groups = colnames(SIGMA_sp)
              Vpart_sp_diag = partial_variances(SIGMA.final[,,sp], B = B_groups)/Var_sp_tot  # same as Torsti, will be the diag elements
              
              S.tmp = as.data.frame(SIGMA.final[,,sp])
              
              # Get all permutations of 2 elements (order matters)
              all_ordered_pairs <- expand.grid(groups, groups, stringsAsFactors = FALSE)
              all_ordered_pairs <- all_ordered_pairs[all_ordered_pairs$Var1 != all_ordered_pairs$Var2, ]
              
              Vpart_sp_nondiag = do.call(rbind, sapply(1:nrow(all_ordered_pairs), function(i){
                S.tmp.i = S.tmp[as.character(all_ordered_pairs[i,]), as.character(all_ordered_pairs[i,])]
                B_groups.i = B_groups[as.character(all_ordered_pairs[i,]), as.character(all_ordered_pairs[i,])]
                Vpart_sp = partial_variances(S.tmp.i, B = B_groups.i)/Var_sp_tot
              }, simplify = F))[,2]  # the two columns are the same but in a different order
              
              
              diag(Vpart_sp) = Vpart_sp_diag
              
              delta <- col(Vpart_sp) - row(Vpart_sp) ## select the non diagonal elements
              Vpart_sp[delta > 0 | delta < 0] <- Vpart_sp_nondiag ## fill it in the columns
              
              colnames(Vpart_sp) = rownames(Vpart_sp) = groups
              Vpart[,,sp,i] = t(Vpart_sp)  ## transpose to get by row
              
            }
          }
          
          ## diagonal variance 
          Vdiag_sp = diag(SIGMA.final[,,sp]) / sum(diag(SIGMA.final[,,sp]))
          Vdiag[,,sp,i] = Vdiag_sp
          
        }
      }else{
        
        #exclude
        name.exclude = colnames(hM$X)[exclude]
        pos = which(groupnames == name.exclude) # new position within our groups
        
        for (sp in 1:ns) {
          #Var_sp_tot = Var_sp_tot  # or
          SIGMA_sp = SIGMA.final[-pos,-pos,sp]
          Var_sp_tot = sum(SIGMA.final[-pos,-pos,sp])
          
          if(marginal == TRUE){
            ## normalized variance
            Vnorm_sp = diag(SIGMA.final[-pos,-pos,sp])/Var_sp_tot  # same as Torsti
            Vnorm[,,sp,i] = Vnorm_sp
            
            
            ## marginal variance
            if(Global == TRUE){
              Vmarg_sp = rowSums(SIGMA.final[-pos,-pos,sp])/Var_sp_tot  # same as Torsti
              Vmarg[,,sp,i] =  Vmarg_sp
              
            }else{
              Vmarg_sp = matrix(NA, n_group-length(pos), n_group-length(pos))
              groups = colnames(SIGMA_sp)
              Vmarg_sp_diag = rowSums(SIGMA.final[-pos,-pos,sp])/Var_sp_tot  # same as Torsti, will be the diag elements
              
              S.tmp = SIGMA.final[-pos,-pos,sp]
                
              Vmarg_sp_nondiag = lapply(sapply(1:nrow(S.tmp), function(i) {
                row <- S.tmp[i, ]
                diag_element <- row[i]
                non_diag_elements <- row[-i]
                diag_element + non_diag_elements
              }, simplify=F, USE.NAMES = F), function(x){x/Var_sp_tot})   # same but for individual groups
              
              diag(Vmarg_sp) = Vmarg_sp_diag
              
              delta <- col(Vmarg_sp) - row(Vmarg_sp) ## select the non diagonal elements
              Vmarg_sp[delta > 0 | delta < 0] <- do.call(cbind, Vmarg_sp_nondiag) ## fill it in the columns
              
              colnames(Vmarg_sp) = rownames(Vmarg_sp) = groups
              Vmarg[,,sp,i] = t(Vmarg_sp)  ## transpose to get by row
            }
            
            
            ## covariance
            SIGMA_sp = SIGMA.final[-pos,-pos,sp]
            SIGMA_sp[! upper.tri(SIGMA_sp)] = NA
            
            Cnorm_sp = SIGMA_sp/Var_sp_tot
            Cnorm[,,sp,i] = Cnorm_sp  # diag + cov
            
            ## partial variance
            ## this may change if we have done the grouping before / to think later about it
            groups = colnames(SIGMA_sp)
            B_groups <- sapply(groups, function(x) { as.integer(colnames(SIGMA_sp) %in% x) })
            dimnames(B_groups) <- list(linear_terms = colnames(SIGMA_sp), groups = colnames(SIGMA_sp))
            
            if(Global == TRUE){
              Vpart_sp = partial_variances(SIGMA.final[-pos,-pos,sp], B = B_groups)/Var_sp_tot
              Vpart[,,sp,i] = Vpart_sp
              
            }else{
              Vpart_sp = matrix(NA, n_group-length(pos), n_group-length(pos))
              groups = colnames(SIGMA_sp)
              Vpart_sp_diag = partial_variances(SIGMA.final[-pos,-pos,sp], B = B_groups)/Var_sp_tot  # same as Torsti, will be the diag elements
              
              S.tmp = as.data.frame(SIGMA.final[-pos,-pos,sp])
              
              # Get all permutations of 2 elements (order matters)
              all_ordered_pairs <- expand.grid(groups, groups, stringsAsFactors = FALSE)
              all_ordered_pairs <- all_ordered_pairs[all_ordered_pairs$Var1 != all_ordered_pairs$Var2, ]

              Vpart_sp_nondiag = do.call(rbind, sapply(1:nrow(all_ordered_pairs), function(i){
                S.tmp.i = S.tmp[as.character(all_ordered_pairs[i,]), as.character(all_ordered_pairs[i,])]
                B_groups.i = B_groups[as.character(all_ordered_pairs[i,]), as.character(all_ordered_pairs[i,])]
                Vpart_sp = partial_variances(S.tmp.i, B = B_groups.i)/Var_sp_tot
              }, simplify = F))[,2]  # the two columns are the same but in a different order
              
              
              diag(Vpart_sp) = Vpart_sp_diag
              
              delta <- col(Vpart_sp) - row(Vpart_sp) ## select the non diagonal elements
              Vpart_sp[delta > 0 | delta < 0] <- Vpart_sp_nondiag ## fill it in the columns
              
              colnames(Vpart_sp) = rownames(Vpart_sp) = groups
              Vpart[,,sp,i] = t(Vpart_sp)  ## transpose to get by row
              
            }
            
            
          }
          
          ## diagonal variance 
          Vdiag_sp = diag(SIGMA.final[-pos,-pos,sp]) / sum(diag(SIGMA.final[-pos,-pos,sp]))
          Vdiag[,,sp,i] = Vdiag_sp
          
        }
      }
      
      
    } # MCMC loop
    
    #........................................................................................
    # Preparing the output
    
    #......................................................................................
    ##  end results from the function
    
    VP = list()
    #VP$R2T = list(Beta=R2T.Beta,Y=R2T.Y)
    VP$group = group
    VP$groupnames = groupnames
    
    if(is.null(exclude)){
      gpnames_X = groupnames
    }else{
      gpnames_X = groupnames[-pos]
      
    }
    
    if(marginal == TRUE){
      if(Global == TRUE){
        dimnames(Vmarg) = list(Value='values', 
                               element=c(gpnames_X, colnames(hM$dfPi)),
                               Species=hM$sp)
        dimnames(Vpart) = list(Value='values', 
                               element=c(gpnames_X, colnames(hM$dfPi)),
                               Species=hM$sp)
      }else{
        dimnames(Vmarg) = list(Value=c(gpnames_X, colnames(hM$dfPi)), 
                               element=c(gpnames_X, colnames(hM$dfPi)),
                               Species=hM$sp)
        dimnames(Vpart) = list(Value=c(gpnames_X, colnames(hM$dfPi)), 
                               element=c(gpnames_X, colnames(hM$dfPi)),
                               Species=hM$sp)
      }
      
      dimnames(Vnorm) = list(Value='values', 
                             element=c(gpnames_X, colnames(hM$dfPi)),
                             Species=hM$sp)
      
      dimnames(Vdiag) = list(Value='values', 
                             element=c(gpnames_X, colnames(hM$dfPi)),
                             Species=hM$sp)
      
      dimnames(Cnorm) = list(element=c(gpnames_X, colnames(hM$dfPi)),
                             element=c(gpnames_X, colnames(hM$dfPi)),
                             Species=hM$sp)
      
      
      VP$Vnorm = Vnorm
      VP$Vmarg = Vmarg
      VP$Cnorm = Cnorm
      VP$Vpart = Vpart
      VP$Vdiag = Vdiag
      
    }else{
      VP$Vnorm = NULL
      VP$Vmarg = NULL
      VP$Cnorm = NULL
      VP$Vpart = NULL
      
      dimnames(Vdiag) = list(Value='values', 
                             element=c(gpnames_X, colnames(hM$dfPi)),
                             Species=hM$sp)
      
      
      VP$Vdiag = Vdiag
      
    }
    
    #. ###### ####  This part ends the normal and / or marginal way to partition the VP
    
  }
  
  return(VP = VP) 
}


#.........................................................................................................................
## function focusing on conditional and within / between partitions summaries

computeVariancePartitioningbwcond = function(hM, group=NULL, groupnames=NULL, start=1, na.ignore=FALSE,
                                             exclude = NULL, cond=NULL, wb=NULL, all=TRUE, conditional=FALSE, marginal=FALSE)
{
  ny = hM$ny
  nc = hM$nc
  nt = hM$nt
  ns = hM$ns
  np = hM$np
  nr = hM$nr
  
  ## This consider if the user has pre entered grouping for linear predictors and if not it creates grouping for each individual linear predictor
  gp = group # keep track on if we had group submitted or not
  gpname = groupnames
  if(is.null(group)){
    ## default: use terms
    if(nc > 1){
      if (is.null(hM$XFormula))
        stop("no XFormula: you must give 'group' and 'groupnames'")
      
      group = attr(hM$X, "assign")
      if (group[1] == 0) # assign (Intercept) to group
        group[1] <- 1
      
      groupnames = attr(terms(hM$XFormula, data = hM$XData), "term.labels")
      gpname = groupnames
      
    } else {
      group = c(1)
      groupnames = hM$covNames[1]
      gpname = groupnames
    }
    
  }
  
  ngroups = max(group)
  
  ## info over MCMC iterations
  postList=poolMcmcChains(hM$postList, start=start)
  
  ## make sure conditions are submitted for the grouping
  if(is.null(cond)){
    stop("You need to define grouping conditions to use the conditional VP summaries")
  }
  
  # +++++++ if conditions => separate
  Gcond = unique(cond)
  n_cond = length(Gcond)
  
  ## within / between prep
  if(!is.null(wb)){
    Hcond <- sapply(Gcond, function(x, g) { as.integer(g == x) }, cond)
    colnames(Hcond) <- Gcond
  }
  
  nsamp = length(postList)#hM$samples
  
  
  ## Initialize arrays ......................................................
  
  if(is.null(wb)){  # conditional
    
    if(is.null(exclude)){
      Vnorm = array(NA, dim = c(1, ngroups+nr, ns, nsamp, n_cond))
      Vmarg = array(NA, dim = c(1, ngroups+nr, ns, nsamp, n_cond))
      Vdiag = array(NA, dim = c(1, ngroups+nr, ns, nsamp, n_cond))
      Vpart = array(NA, dim = c(1, ngroups+nr, ns, nsamp, n_cond))
      Cnorm = array(NA, dim = c(ngroups+nr, ngroups+nr, ns, nsamp, n_cond))
    }else{
      Vnorm = array(NA, dim = c(1, ngroups+nr-length(exclude), ns, nsamp, n_cond))   # if we exclude certain covariates like effort for the calculation
      Vmarg = array(NA, dim = c(1, ngroups+nr-length(exclude), ns, nsamp, n_cond))
      Vdiag = array(NA, dim = c(1, ngroups+nr-length(exclude), ns, nsamp, n_cond))
      Vpart = array(NA, dim = c(1, ngroups+nr-length(exclude), ns, nsamp, n_cond))
      Cnorm = array(NA, dim = c(ngroups+nr-length(exclude), ngroups+nr-length(exclude), ns, nsamp, n_cond))
    }
  }else{ # within / between
    
    if(is.null(exclude)){
      Vnorm = array(NA, dim = c(1, ngroups+nr, ns, nsamp, 2))
      Vmarg = array(NA, dim = c(1, ngroups+nr, ns, nsamp, 2))
      Vdiag = array(NA, dim = c(1, ngroups+nr, ns, nsamp, 2))
      Vpart = array(NA, dim = c(1, ngroups+nr, ns, nsamp, 2))
      Cnorm = array(NA, dim = c(ngroups+nr, ngroups+nr, ns, nsamp, 2))
    }else{
      Vnorm = array(NA, dim = c(1, ngroups+nr-length(exclude), ns, nsamp, 2))   # if we exclude certain covariates like effort for the calculation
      Vmarg = array(NA, dim = c(1, ngroups+nr-length(exclude), ns, nsamp, 2))
      Vdiag = array(NA, dim = c(1, ngroups+nr-length(exclude), ns, nsamp, 2))
      Vpart = array(NA, dim = c(1, ngroups+nr-length(exclude), ns, nsamp, 2))
      Cnorm = array(NA, dim = c(ngroups+nr-length(exclude), ngroups+nr-length(exclude), ns, nsamp, 2))
    }
  }
  
  
  ##run over MCMC chains NB
  for (i in 1:nsamp){   # to test for now
    print(i)
    
    # 1 #++++++++++++++++++++++++++++++++++++++   FIXED part
    
    if(is.null(wb)){   # 1 # calculating posterior sample variance or conditional variance
      
      ### under conditions
      ## estimates posterior sample
      Beta.i = postList[[i]]$Beta
      Lambdas.i = postList[[i]]$Lambda
      Etas.i = postList[[i]]$Eta
      
      
      ## calculation for all species
      
      # temporary names to use
      Beta.x = Beta.i
      
      LF.array = array(NA, c(dim(hM$X)[1], dim(hM$X)[2], ns))
      for (kk in 1:dim(hM$X)[2]){
        LF.temp = as.vector(hM$X[,kk])%*%t(as.vector(Beta.x[kk,]))
        LF.array[,kk,] = LF.temp
      }
      
      
      # 2 #++++++++++++++++++++++++++++++++++++++   RANDOM part
      
      RF.array = array(NA, dim=c(dim(hM$Pi)[1], dim(hM$Pi)[2], ns))
      for (kk in 1:length(Lambdas.i)){ ## across random effect
        Pi.temp = matrix(0, nrow=dim(hM$Pi)[1], ncol=length(unique(hM$Pi[,kk])))
        for (jj in 1:dim(hM$Pi)[1]){  ## across sampling occasions
          Pi.temp[jj,hM$Pi[jj,kk]] = 1
        }
        RF.temp = Pi.temp%*%Etas.i[[kk]]%*%Lambdas.i[[kk]]
        RF.array[,kk,] = RF.temp
      }
      
      
      # 3 #++++++++++++++++++++++++++++++++++++++   FIXED - RANDOM covar part  ##does not exist  if conditional == FALSE in this branch
      
      # Combined covariance partition matrix of all fixed effects and random effects
      dim(LF.array)
      dim(RF.array)
      
      SIGMA = array(NA, dim=c(dim(hM$X)[2]+dim(hM$Pi)[2], dim(hM$X)[2]+dim(hM$Pi)[2], n_cond, ns))
      for (jj in 1:ns) {
        
        LRF.array = cbind(LF.array[,,jj], RF.array[,,jj])
        SIGMA[,,,jj] = sapply(Gcond, function(x, A, g) { cov(A[g == x,])}, LRF.array, cond, simplify = "array")
        dimnames(SIGMA[,,,jj]) <- list(colnames(SIGMA[,,,jj]), colnames(SIGMA[,,,jj]), condition = Gcond)
      }
      dimnames(SIGMA) = list(element=c(colnames(hM$X), colnames(hM$Pi)), 
                             element=c(colnames(hM$X), colnames(hM$Pi)),
                             cond = Gcond, Species=hM$sp)
      
      SIGMA = aperm(SIGMA, c(1,2,4,3)) ## condition at the end
      
      # 4 # -------- Do we need grouping?
      # In any case we group the intercept so:
      ## adapted from Torsti's example
      
      linear_terms = c(colnames(hM$XScaled), colnames(hM$Pi))
      group.tmp = c(group, seq(from = max(group)+1, length.out= dim(hM$dfPi)[2]))
      groups <- split(linear_terms, group.tmp)
      names(groups) = c(groupnames, colnames(hM$Pi))
      
      B <- sapply(groups, function(x) { as.integer(linear_terms %in% x) })
      dimnames(B) <- list(linear_terms = linear_terms, groups = colnames(B))
      
      
      n_group = max(group.tmp)
      K_B <- apply(SIGMA, c("Species", 'cond'), function(x, B) { t(B) %*% x %*% B }, B)
      dim(K_B) <- c(n_group, n_group, ns, n_cond)
      dimnames(K_B) <- list(rows = c(groupnames, colnames(hM$Pi)), cols = c(groupnames, colnames(hM$Pi)), 
                            species = hM$sp, cond = Gcond)
      
      SIGMA.final = K_B
      
      #.........................................................................................................
      # Last part of the processing
      
      ### Matrix K linear terms  ready *******************************************************************************
      ### Similar to Torsti's paper *******************************************************************************
      
      
      # Total and partition calculation
      for (ncd in 1:n_cond){
        if(is.null(exclude)){
          for (sp in 1:ns) {
            #Var_sp_tot = Var_sp_tot  # or
            SIGMA.final_sp = SIGMA.final[,,sp,ncd]
            Var_sp_tot = sum(SIGMA.final[,,sp,ncd])
            
            if(marginal == TRUE){
              ## normalized variance
              Vnorm_sp = diag(SIGMA.final[,,sp,ncd])/Var_sp_tot  # same as Torsti
              Vnorm[,,sp,i,ncd] = Vnorm_sp
              
              ## marginal variance
              Vmarg_sp = rowSums(SIGMA.final[,,sp,ncd])/Var_sp_tot  # same as Torsti
              Vmarg[,,sp,i,ncd] =  Vmarg_sp
              
              ## covariance
              SIGMA.final_sp = SIGMA.final[,,sp,ncd]
              SIGMA.final_sp[! upper.tri(SIGMA.final_sp)] = NA
              
              Cnorm_sp = SIGMA.final_sp/Var_sp_tot
              Cnorm[,,sp,i,ncd] = Cnorm_sp
              
              ## partial variance
              ## this may change if we have done the grouping before / to think later about it
              
              remove = which(colSums(SIGMA.final[,,sp,ncd] != 0) == 0)
              if(!is.na(remove)){
                groups = colnames(SIGMA.final_sp[-c(remove), -c(remove)])
                B_groups <- sapply(groups, function(x) { as.integer(colnames(SIGMA.final_sp)[-(remove)] %in% x) })
                dimnames(B_groups) <- list(linear_terms = colnames(SIGMA.final_sp)[-(remove)], groups = colnames(SIGMA.final_sp)[-(remove)])
                
                Vpart_sp = partial_variances(SIGMA.final[-c(remove),-c(remove),sp,ncd], B = B_groups)/Var_sp_tot
                Vpart[-c(remove),-c(remove),sp,i,ncd] = Vpart_sp
                
              }
            }
            
            ##
            Vdiag_sp = diag(SIGMA.final[,,sp,ncd]) / sum(diag(SIGMA.final[,,sp,ncd]))
            Vdiag[,,sp,i,ncd] = Vdiag_sp
          }
          
        }else{
          
          #exclude
          name.exclude = colnames(hM$X)[exclude]
          pos = which(groupnames ==name.exclude) # new position within our groups
          
          for (sp in 1:ns) {
            #Var_sp_tot = Var_sp_tot  # or
            SIGMA.final_sp = SIGMA.final[-pos,-pos,sp,ncd]
            Var_sp_tot = sum(SIGMA.final[-pos,-pos,sp,ncd])
            
            if(marginal==TRUE){
              ## normalized variance
              Vnorm_sp = diag(SIGMA.final[-pos,-pos,sp,ncd])/Var_sp_tot  # same as Torsti
              Vnorm[,,sp,i,ncd] = Vnorm_sp
              
              ## marginal variance
              Vmarg_sp = rowSums(SIGMA.final[-pos,-pos,sp,ncd])/Var_sp_tot  # same as Torsti
              Vmarg[,,sp,i,ncd] =  Vmarg_sp
              
              ## covariance
              SIGMA.final_sp = SIGMA.final[-pos,-pos,sp,ncd]
              SIGMA.final_sp[! upper.tri(SIGMA.final_sp)] = NA
              
              Cnorm_sp = SIGMA.final_sp/Var_sp_tot
              Cnorm[,,sp,i,ncd] = Cnorm_sp
              
              ## partial variance
              ## this may change if we have done the grouping before / to think later about it
              
              remove = which(colSums(SIGMA.final[,,sp,ncd] != 0) == 0)
              if(!is.na(remove)){
                groups = colnames(SIGMA.final_sp[-c(remove-length(pos)), -c(remove-length(pos))])
                B_groups <- sapply(groups, function(x) { as.integer(colnames(SIGMA.final_sp)[-(remove-length(pos))] %in% x) })
                dimnames(B_groups) <- list(linear_terms = colnames(SIGMA.final_sp)[-(remove-length(pos))], groups = colnames(SIGMA.final_sp)[-(remove-length(pos))])
                
                Vpart_sp = partial_variances(SIGMA.final[-c(pos, remove),-c(pos, remove),sp,ncd], B = B_groups)/Var_sp_tot
                Vpart[-c(pos, remove),-c(pos, remove),sp,i,ncd] = Vpart_sp
                
              }  
            }
            
            ## diagonal variance should be the same as original Hmsc
            Vdiag_sp = diag(SIGMA.final[-pos,-pos,sp,ncd]) / sum(diag(SIGMA.final[-pos,-pos,sp,ncd]))
            Vdiag[,,sp,i,ncd] = Vdiag_sp
          }
        }
        
      } # Finalized for each condition
      
      
    }else{  # 2 # calculating posterior sample btw - within variance
      
      
      Beta.i = postList[[i]]$Beta
      Lambdas.i = postList[[i]]$Lambda
      Etas.i = postList[[i]]$Eta
      
      
      ## 1    fixed part   .....................................
      # temporary names to use
      Beta.x = Beta.i
      
      LF.array = array(NA, c(dim(hM$X)[1], dim(hM$X)[2], ns))
      for (kk in 1:dim(hM$X)[2]){
        LF.temp = as.vector(hM$X[,kk])%*%t(as.vector(Beta.x[kk,]))
        LF.array[,kk,] = LF.temp
      }
      
      ##  2-  random  part ............................................
      
      RF.array = array(NA, dim=c(dim(hM$Pi)[1], dim(hM$Pi)[2], ns))
      for (kk in 1:length(Lambdas.i)){ ## across random effect
        Pi.temp = matrix(0, nrow=dim(hM$Pi)[1], ncol=length(unique(hM$Pi[,kk])))
        for (jj in 1:dim(hM$Pi)[1]){  ## across sampling occasions
          Pi.temp[jj,hM$Pi[jj, kk]] = 1
        }
        RF.temp = Pi.temp%*%Etas.i[[kk]]%*%Lambdas.i[[kk]]
        RF.array[,kk,] = RF.temp
      }
      
      # 3 #++++++++++++++++++++++++++++++++++++++   FIXED - RANDOM covar part
      
      #Combined covariance partition matrix of all fixed effects and random effects
      dim(LF.array)
      dim(RF.array)
      
      Hcond <- sapply(Gcond, function(x, g) { as.integer(g == x) }, cond)
      colnames(Hcond) <- Gcond
      
      SIGMA = array(NA, dim=c(dim(hM$X)[2]+dim(hM$Pi)[2], dim(hM$X)[2]+dim(hM$Pi)[2], 2, ns))
      for (sp in 1:ns) {
        LRF.array = cbind(LF.array[,,sp], RF.array[,,sp])
        
        A_cm <- Hcond %*% diag(1 / colSums(Hcond)) %*% t(Hcond) %*% LRF.array
        K_wit <- cov(LRF.array - A_cm)
        K_btw <- cov(A_cm)
        
        SIGMA[,,,sp] <- abind::abind(K_wit, K_btw, rev.along = 0)
        
      }
      dimnames(SIGMA) = list(element=c(colnames(hM$X), colnames(hM$Pi)), 
                             element=c(colnames(hM$X), colnames(hM$Pi)),
                             partition = c("within", "btw"),
                             Species=hM$sp)
      
      SIGMA = aperm(SIGMA, c(1,2,4,3)) ## partition at the end
      
      # 4 # -------- Do we need grouping?
      # In any case we group the intercept so:
      ## adapted from Torsti's example
      
      linear_terms = c(colnames(hM$XScaled), colnames(hM$Pi))
      group.tmp = c(group, seq(from = max(group)+1, length.out= dim(hM$dfPi)[2]))
      groups <- split(linear_terms, group.tmp)
      names(groups) = c(groupnames, colnames(hM$Pi))
      
      B <- sapply(groups, function(x) { as.integer(linear_terms %in% x) })
      dimnames(B) <- list(linear_terms = linear_terms, groups = colnames(B))
      
      
      n_group = max(group.tmp)
      K_B <- apply(SIGMA, c("Species", 'partition'), function(x, B) { t(B) %*% x %*% B }, B)
      dim(K_B) <- c(n_group, n_group, ns, 2)
      dimnames(K_B) <- list(rows = c(groupnames, colnames(hM$Pi)), cols = c(groupnames, colnames(hM$Pi)), species = hM$sp, partition = c("within", "btw"))
      
      SIGMA.final = K_B
      
      #.........................................................................................................
      # Last part of the processing
      
      ### Matrix K linear terms  ready *******************************************************************************
      ### Similar to Torsti's paper *******************************************************************************
      
      for (wbind in 1:2) {
        # Total and partition calculation
        
        if(is.null(exclude)){
          for (sp in 1:ns) {
            SIGMA.final_sp = SIGMA.final[,,sp,wbind]
            Var_sp_tot = sum(SIGMA.final[,,sp,wbind]) #### take the whole triangle
            
            if(marginal == TRUE){
              ## normalized variance
              Vnorm_sp = diag(SIGMA.final[,,sp,wbind])/Var_sp_tot  # same as Torsti
              Vnorm[,,sp,i,wbind] = Vnorm_sp
              
              
              ## marginal variance
              Vmarg_sp = rowSums(SIGMA.final[,,sp,wbind])/Var_sp_tot  # same as Torsti
              Vmarg[,,sp,i,wbind] =  Vmarg_sp
              
              
              ## covariance
              SIGMA.final_sp = SIGMA.final[,,sp,wbind]
              SIGMA.final_sp[! upper.tri(SIGMA.final_sp)] = NA
              
              Cnorm_sp = SIGMA.final_sp/Var_sp_tot
              Cnorm[,,sp,i,wbind] = Cnorm_sp
              
              ## partial variance
              ## If the grouping is based on the random effect it creates problems so the solution (remove that random effect from the matrix (values were 0 anyways))
              In.random = apply(hM$studyDesign, 2, function(x){Gcond %in% x})[1,]
              
              remove = which(colnames(SIGMA.final[,,sp,wbind]) %in% names(which(In.random==TRUE)))
              if(!is.na(remove)){
                groups = colnames(SIGMA.final_sp[-c(remove), -c(remove)])
                B_groups <- sapply(groups, function(x) { as.integer(colnames(SIGMA.final_sp)[-(remove)] %in% x) })
                dimnames(B_groups) <- list(linear_terms = colnames(SIGMA.final_sp)[-(remove)], groups = colnames(SIGMA.final_sp)[-(remove)])
                
                res <- try(partial_variances(SIGMA.final[-c(remove),-c(remove),sp,wbind], B = B_groups)/Var_sp_tot)
                if(inherits(res, "try-error"))
                {
                  Vpart_sp = rep(NA, max(group)+nr-1)
                  next
                }else{
                  Vpart_sp = partial_variances(SIGMA.final[-c(remove),-c(remove),sp,wbind], B = B_groups)/Var_sp_tot
                }
                Vpart[-remove,-remove,sp,i,wbind] = Vpart_sp
                
              }
            }
            
            ## diagonal variance should be the same as original Hmsc
            Vdiag_sp = diag(SIGMA.final[,,sp,wbind]) / sum(diag(SIGMA.final[,,sp,wbind]))
            Vdiag[,,sp,i,wbind] = Vdiag_sp
            
          }
        }else{
          
          #exclude
          name.exclude = colnames(hM$X)[exclude]
          pos = which(groupnames ==name.exclude) # new position within our groups
          
          for (sp in 1:ns) {
            #Var_sp_tot = Var_sp_tot  # or
            SIGMA.final_sp = SIGMA.final[-pos,-pos,sp,wbind]
            #Var_sp_tot = sum(diag(SIGMA.final[-pos,-pos,sp])) + sum(SIGMA.final_sp[upper.tri(SIGMA.final_sp)], na.rm=T)
            Var_sp_tot = sum(SIGMA.final[-pos,-pos,sp,wbind])
            
            if(marginal ==TRUE){
              ## normalized variance
              Vnorm_sp = diag(SIGMA.final[-pos,-pos,sp,wbind])/Var_sp_tot  # same as Torsti
              Vnorm[,,sp,i,wbind] = Vnorm_sp
              
              
              ## marginal variance
              Vmarg_sp = rowSums(SIGMA.final[-pos,-pos,sp,wbind])/Var_sp_tot  # same as Torsti
              Vmarg[,,sp,i,wbind] =  Vmarg_sp
              
              
              ## covariance
              SIGMA.final_sp = SIGMA.final[-pos,-pos,sp,wbind]
              SIGMA.final_sp[! upper.tri(SIGMA.final_sp)] = NA
              
              Cnorm_sp = SIGMA.final_sp/Var_sp_tot
              Cnorm[,,sp,i,wbind] = Cnorm_sp
              
              ## partial variance
              ## this may change if we have done the grouping before / to think later about it
              In.random = apply(hM$studyDesign, 2, function(x){Gcond %in% x})[1,]
              
              remove = which(colnames(SIGMA.final[,,sp,wbind]) %in% names(which(In.random==TRUE)))
              if(!is.na(remove)){
                groups = colnames(SIGMA.final_sp[-c(remove-length(pos)), -c(remove-length(pos))])
                B_groups <- sapply(groups, function(x) { as.integer(colnames(SIGMA.final_sp)[-(remove-length(pos))] %in% x) })
                dimnames(B_groups) <- list(linear_terms = colnames(SIGMA.final_sp)[-(remove-length(pos))], groups = colnames(SIGMA.final_sp)[-(remove-length(pos))])
                
                res <- try(partial_variances(SIGMA.final[-c(pos, remove),-c(pos, remove),sp,wbind], B = B_groups)/Var_sp_tot)
                if(inherits(res, "try-error"))
                {
                  Vpart_sp = rep(NA, max(group)+nr-1)
                  next
                }else{
                  Vpart_sp = partial_variances(SIGMA.final[-c(pos, remove),-c(pos, remove),sp,wbind], B = B_groups)/Var_sp_tot
                }
                
                Vpart[-c(pos, remove),-c(pos, remove),sp,i,wbind] = Vpart_sp
                
              }
              
            }
            
            ## diagonal variance should be the same as original Hmsc
            Vdiag_sp = diag(SIGMA.final[-pos,-pos,sp,wbind]) / sum(diag(SIGMA.final[-pos,-pos,sp,wbind]))
            Vdiag[,,sp,i,wbind] = Vdiag_sp
            
          }
        }
        
      } # wind 1 or 2 for btw wit
      
    } # for btw and within conditions
    
    
  } # MCMC loop
  
  #........................................................................................
  # Preparing the output
  
  #......................................................................................
  ##  end results from the function
  
  VP = list()
  #VP$R2T = list(Beta=R2T.Beta,Y=R2T.Y)
  VP$group = group
  VP$groupnames = groupnames
  
  if(is.null(exclude)){
    gpnames_X = groupnames
  }else{
    gpnames_X = groupnames[-pos]
    
  }
  
  rd_names = colnames(hM$dfPi)
  
  
  if(marginal==TRUE){
    if(is.null(wb)){
      dimnames(Vnorm) = list(Value='values',
                             element=c(gpnames_X, rd_names),
                             Species=hM$sp, iteration = 1:nsamp,
                             group = Gcond)
      dimnames(Vmarg) = list(Value='values',
                             element=c(gpnames_X, rd_names),
                             Species=hM$sp, iteration = 1:nsamp,
                             group = Gcond)
      dimnames(Vpart) = list(Value='values',
                             element=c(gpnames_X, rd_names),
                             Species=hM$sp, iteration = 1:nsamp,
                             group = Gcond)
      dimnames(Vdiag) = list(Value='values',
                             element=c(gpnames_X, rd_names),
                             Species=hM$sp, iteration = 1:nsamp,
                             group = Gcond)
      
      dimnames(Cnorm) = list(element=c(gpnames_X, rd_names),
                             element=c(gpnames_X, rd_names),
                             Species=hM$sp, iteration = 1:nsamp,
                             group = Gcond)
    }else{
      dimnames(Vnorm) = list(Value='values',
                             element=c(gpnames_X, rd_names),
                             Species=hM$sp, iteration = 1:nsamp,
                             group = c('wit', 'btw'))
      dimnames(Vmarg) = list(Value='values',
                             element=c(gpnames_X, rd_names),
                             Species=hM$sp, iteration = 1:nsamp,
                             group = c('wit', 'btw'))
      dimnames(Vpart) = list(Value='values',
                             element=c(gpnames_X, rd_names),
                             Species=hM$sp, iteration = 1:nsamp,
                             group = c('wit', 'btw'))
      dimnames(Vdiag) = list(Value='values',
                             element=c(gpnames_X, rd_names),
                             Species=hM$sp, iteration = 1:nsamp,
                             group = c('wit', 'btw'))
      
      dimnames(Cnorm) = list(element=c(gpnames_X, rd_names),
                             element=c(gpnames_X, rd_names),
                             Species=hM$sp, iteration = 1:nsamp,
                             group = c('wit', 'btw'))
    }
    
    
    VP$Vnorm = Vnorm
    VP$Vmarg = Vmarg
    VP$Cnorm = Cnorm
    VP$Vpart = Vpart
    VP$Vdiag = Vdiag
  }else{
    if(is.null(wb)){
      
      dimnames(Vdiag) = list(Value='values',
                             element=c(gpnames_X, rd_names),
                             Species=hM$sp, iteration = 1:nsamp,
                             group = Gcond)
      
    }else{
      
      dimnames(Vdiag) = list(Value='values',
                             element=c(gpnames_X, rd_names),
                             Species=hM$sp, iteration = 1:nsamp,
                             group = c('wit', 'btw'))
      
    }
    
    
    VP$Vnorm = NULL
    VP$Vmarg = NULL
    VP$Cnorm = NULL
    VP$Vpart = NULL
    VP$Vdiag = Vdiag
  }
  
  
  return(VP = VP)
}


#.........................................................................................................................
# partial variance function from Schulz's et al (2025) paper but with changes
partial_variance <- function(K, idx) {
  P <- sum(K[idx,idx] - K[idx,!idx] %*% solve(K[!idx,!idx], tol = 1e-40) %*% K[!idx,idx])
}

partial_variances <- function(K, B) {
  P <- apply(B, 1, function(x, K) { partial_variance(K, x == 1 )}, K)
}
