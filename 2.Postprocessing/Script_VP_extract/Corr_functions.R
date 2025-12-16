#
#
#    Extra functions to prepare correlations
#
#....................................................................

## extract linear terms correlations - similar way to do as for variance partition from HMSC
Extract_corr = function(hM, group=NULL, groupnames=NULL, start=1, na.ignore=FALSE){
  
  ngroups = max(group)
  #ny = hM$ny
  #nc = hM$nc
  #nt = hM$nt
  ns = hM$ns
  #np = hM$np
  nr = hM$nr
  
  ## info over MCMC iterations
  postList = poolMcmcChains(hM$postList)
  nsamp = length(postList)
  
  XDATA_bf = hM$X
  
  #fixed.cor = array(NA, dim = c(ngroups, ngroups, ns, nsamp))
  #rd.cor = array(NA, dim = c(ngroups, 1, ns, nsamp))
  FRcor = array(NA, dim = c(ngroups+nr, ngroups+nr, ns, nsamp))
  
  # over posterior samples
  for (i in 1:nsamp){
    
    #keeping track of the chains
    print(i)
    i
    
    ## estimates posterior sample
    Beta.i = postList[[i]]$Beta
    Lambdas.i = postList[[i]]$Lambda
    Etas.i = postList[[i]]$Eta
    
    ns = hM$ns
    nr = hM$nr
    
    ## prepare covariance - for fixed part only
    cMA = cov(hM$X)
    
    
    # 1 #++++++++++++++++++++++++++++++++++++++   FIXED part 
    
    fixedsplit1.marg = array(0, dim=c(ngroups, ngroups, ns))

    ## calculation for each species
    for (j in 1:ns){
      
      Beta.x = Beta.i
      cM.x = cMA
      
      for (k in 1:ngroups){
        sel = (group==k)
        selgp =  k
        alltsel = unique(group)#[-k]  # we actually need all groups not only the non selected ones to get also the diagonal (VP)
        
        
        fpart.marg = sapply(alltsel, function(x, A, B, g, y) {(t(B[g == x, j])%*%A[g == x, g == y]%*%B[g == y, j])}, 
                            cM.x, Beta.x, group, selgp, simplify = "array")
        
        
        fixedsplit1.marg[k,,j] = fixedsplit1.marg[k,,j] + fpart.marg
        dimnames(fixedsplit1.marg) <- list(groups = unique(group), groups = unique(group), species=colnames(Beta.i))
        
        
      }
      #fixed.cor[,,j,i] = cov2cor(fixedsplit1.marg[,,j])      
      
    }
    
    
    
    # 2 #++++++++++++++++++++++++++++++++++++++   RANDOM - fixed part 
    
    ## for now I repeat depending on how how we choose took at VP(random)
    H.i = GetRandparam2(hM = hM, Lambdas = Lambdas.i, Etas = Etas.i)
    
    H.i.lev = H.i$H
    dim(H.i.lev)
    
    dimnames(H.i.lev) = list(occasions=rownames(H.i$H), random=colnames(hM$dfPi),
                             Species=hM$sp)
    
    # variance group of covariates - random effects
    
    rd1fix.marg = array(0, dim=c(length(unique(group)), nr, ns))
    # random
    
    for (j in 1:ns) {
      
      if(nr > 1){
        colnames(H.i.lev[,,j]) = colnames(hM$dfPi)
        H.marg = cbind(hM$X, H.i.lev[,,j])
      }else{
        H.marg = cbind(hM$X, H.i.lev[,,j])
      }
      
      dim(H.marg)
      G = unique(group)
      
      #***** Testing
      r.tmp = cov(H.marg)
      
      colnames(r.tmp) = c(colnames(hM$X), colnames(hM$dfPi))
      rownames(r.tmp) = c(colnames(hM$X), colnames(hM$dfPi))
      
      r.tmp2 = r.tmp[1:dim(hM$X)[2], -c(1:dim(hM$X)[2])]
      
      
      Beta.x = Beta.i
      Beta.xj = Beta.x[, j]
      
      # dimensions of lambdas.i vary so need to adapt:
      # Step 1: Determine the maximum dimensions
      max_rows <- max(sapply(Lambdas.i, nrow))
      max_cols <- max(sapply(Lambdas.i, ncol))
      
      # Step 2: Pad the matrices
      padded_mat_list <- lapply(Lambdas.i, function(mat) {
        padded_mat <- matrix(NA, nrow = max_rows, ncol = max_cols)
        padded_mat[1:nrow(mat), 1:ncol(mat)] <- mat
        return(padded_mat)
      })
      
      # Step 3: Combine into an array
      lambda.arr <- abind::abind(padded_mat_list, along = 3)
      
      lambda.arr2 = aperm(lambda.arr, c(3,1,2))
      
      #gtest = G[1]
      #t(Beta.xj[group == gtest])%*%r.tmp2[group == gtest,]%*%lambda.arr2[,,j]
      #vect.text = t(Beta.xj[group == gtest])%*%r.tmp2[group == gtest,]
      #rowSums(sweep(lambda.arr2[,,j], MARGIN=1, vect.text, '*'), na.rm=T)
      
      
      # with coef
      rtotal.marg = matrix(NA, ncol = nr, nrow = ngroups)
      for (g in 1:length(G)) {
        vect.text = t(Beta.xj[group == g])%*%r.tmp2[group == g,]
        rtotal.marg[g,] = rowSums(sweep(lambda.arr2[,,j], MARGIN=1, vect.text, '*'), na.rm=T)
      }
      colnames(rtotal.marg) = colnames(r.tmp2)
      rownames(rtotal.marg) = groupnames
      
      rd1fix.marg[,,j] = rtotal.marg
      #rd.cor[,,j,i] = cov2cor(rd1fix.marg[,,j])  
    }
    
    
    # 2 #++++++++++++++++++++++++++++++++++++++   RANDOM part 
    
    
    ## New way to look at random factors VP
    # evaluate the random parameters
    H.i = GetRandparam2(hM=hM, Lambdas = Lambdas.i, Etas = Etas.i)
    
    H.i.lev = H.i$H
    dim(H.i.lev)
    dimnames(H.i.lev) = list(occasions=rownames(H.i$H), random=colnames(hM$dfPi),
                             Species=hM$sp)
    
    rdsplit1 = array(0, dim=c(nr, nr, ns))
    
    if(nr > 1){
      for (j in 1:ns) {
        
        # Covariance matrix - Random effects
        r.tmp = cov(H.i.lev[,,j])
        dimnames(r.tmp) = list(rand=colnames(hM$dfPi), rand=colnames(hM$dfPi))
        
        # Get lambdas into a good format for multiplication
        # dimensions of lambdas.i vary so need to adapt:
        # Step 1: Determine the maximum dimensions
        max_rows <- max(sapply(Lambdas.i, nrow))
        max_cols <- max(sapply(Lambdas.i, ncol))
        
        # Step 2: Pad the matrices (dimensions match the maximum latent factor levels number but gets NA is a latent factor has less levels)
        padded_mat_list <- lapply(Lambdas.i, function(mat) {
          padded_mat <- matrix(NA, nrow = max_rows, ncol = max_cols)
          padded_mat[1:nrow(mat), 1:ncol(mat)] <- mat
          return(padded_mat)
        })
        
        # Step 3: Combine into an array
        lambda.arr <- abind::abind(padded_mat_list, along = 3)
        lambda.arr2 = aperm(lambda.arr, c(3,1,2))  ## dim rd effects x max latent factor levels x species
        
        
        # for a species at a time
        LB.j = lambda.arr2[,,j]
        
        rdsplt.tmp = matrix(ncol=3, nrow = 3)
        # Calculate the aramater covariance multiplication for each latent factor level and then sums them
        for (r in 1:dim(r.tmp)[1]) {
          selrd =  r ## the group focus/selected
          allrd = 1: dim(r.tmp)[1]
          
          # we do the same as what we did for the fixed effect groups but for each level of latent factor and then sum them
          rdsplt.tmp[r,] = rowSums(sapply(1:ncol(LB.j), function(i, A, C) {
            LB <- as.matrix(A[, i, drop = FALSE])
            
            sapply(allrd, function(x, A, B, g, y) {(t(B[g == x])%*%A[g == x, g == y]%*%B[g == y])}, 
                   C, LB, allrd, selrd, simplify = "array")
            
          }, LB.j, r.tmp, simplify = 'array'), na.rm=T)
          
        }
        
        rdsplit1[,,j] = rdsplt.tmp
        dimnames(rdsplit1)[[1]] = colnames(hM$dfPi)
        dimnames(rdsplit1)[[2]] = colnames(hM$dfPi)  
        
      }
      dimnames(rdsplit1)[[3]] = hM$sp

    }else{ # not done yet for one rd effect
      for (j in 1:ns) {
        r.tmp = cov(as.matrix(H.i.lev[,,j]))
        #names(r.tmp) = list(rand=names(hM$ranLevels), rand=names(hM$ranLevels))
        
        #rtotal = apply(r.tmp, 'condition', function(x) {sum(x)}, simplify = "array")
        #rd1[j,] = unlist(rtotal)
        
        rdsplit1[,,j] = r.tmp
        
      }
      dimnames(rdsplit1)[[3]] = hM$sp
    }
    
    
    
    # 4 #++++++++++++++++++++++++++++++++++++++   combine for corr
    
    dim(rdsplit1)
    dim(rd1fix.marg)
    dim(fixedsplit1.marg)
    
    FRmarg = array(NA, dim = c(ngroups+nr, ngroups+nr, ns))
    for (j in 1:ns) {
      FRmarg[,,j] = rbind(cbind(fixedsplit1.marg[,,j], rd1fix.marg[,,j]),
                          cbind(t(rd1fix.marg[,,j]), rdsplit1[,,j]))
    }
    dimnames(FRmarg) = list(element=c(groupnames, colnames(rdsplit1)), 
                            element=c(groupnames, colnames(rdsplit1)),
                            Species=hM$sp)
    
    ## 5 -----------------------------------------  Correlation like in Torsti's paper
    
    for (j in 1:ns) {
      FRcor[,,j, i] = cov2cor(FRmarg[,,j])
    }
    
  }
  
  return(FRcor = FRcor)
}

# function to get random parameters
GetRandparam2 = function(hM, Lambdas, Etas){
  
  nr = hM$nr
  np = hM$np
  ns = hM$ns
  
  # random information from design 
  mat.random = cbind(1:dim(hM$dfPi)[1], hM$studyDesign[1:dim(hM$dfPi)[2]])
  rand_eff = colnames(hM$dfPi)
  mat.rand = list()
  
  for(x in rand_eff) {
    mat.random = cbind(mat.random, value=1) 
    mat.random[,x] = paste0(x, mat.random[,x])
    mat.rand[[x]] = dcast(mat.random, 
                          as.formula(paste0(paste(colnames(mat.random)[1]), '~ ', x)), fill=0)
    mat.rand[[x]] = mat.rand[[x]][, -1]
  }
  
  lapply(mat.rand, dim) # to verify
  
  
  H = array(NA, dim=c(nrow(hM$X), length(rand_eff), ns))
  Lrand_list = list()
  for (level in seq_len(nr)) {
    Lambda = Lambdas[[level]]
    nf = dim(Lambda)[[1]]
    Eta = Etas[[level]]
    
    # random effect from HMSC
    Lrand = matrix(NA, ncol = ns, nrow = np[level])
    Lr.fact = lapply(1:ns, matrix, data= NA, nrow=np[level], ncol=nf)
    for (j in 1:ns){
      for (factor in 1:nf) {
        Lr.fact[[j]][, factor] = Eta[,factor]#*Lambda[factor,j]
      }
      Lrand[, j] = rowSums(Lr.fact[[j]])
    }
    
    Lrand_list[[level]] = Lrand
    # getting H
    H[, level, ] = (as.matrix(mat.rand[[level]]) %*% Lrand_list[[level]])
    ### 3random effects list of matrices (xx occasions x yy species)
  }
  return(list(H=H, Lrand_list = Lrand_list))
}


