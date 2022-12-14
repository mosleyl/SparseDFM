
##  Estimate a Sparse Dynamic Factor Model  (see doc for details)
##
##  Main function to allow estimation of a standard DFM and a sparse DFM
##    - option for IID or AR(1) errors 
##    - option for Kalman Filter/Smoother estimation using standard multivariate equations or fast univariate filtering equations
##
##  X_t = Lambda.tilde*F.tilde_t + eta_t,  eta_t ~ N(0,Sigma.eta)
##  F.tilde_t = A.tilde*F.tilde_{t-1} + u.tilde_t,  u.tilde_t ~ N(0,Sigma.u.tilde)
##
##  If IID errors:  
##
##    - Latent state of dim n x r: F.tilde_t = F_t for r factors F_t, t = 1,...,n.
##    - Loadings of dim p x r: Lambda.tilde = Lambda where Lambda are the loadings (possibly) assumed to be sparse. 
##    - State transition matrix of dim r x r: A.tilde = A.
##    - Idio error covariance of dim p x p: Sigma.eta = Sigma_epsilon.
##    - State equation error covariance of dim r x r: Sigma.u.tilde = Sigma_u. 
##
##
##  If AR(1) errors: 
##
##    - k = r + p 
##    - Latent state of dim n x k: F.tilde_t = [F_t e_t]^T for r factors F_t and p idio errors e_t., t = 1,...,n.
##    - Loadings of dim k x r: Lambda.tilde = [Lambda I_p] where Lambda are the loadings assumed to be sparse. 
##    - State transition matrix of dim k x k: A.tilde = [A 0; 0 Phi].
##    - Idio error covariance of dim p x p: Sigma.eta = kappa*I_p, kappa set to 1e-4.
##    - State equation error covariance of dim k x k: Sigma.u.tilde = [Sigma_u 0; 0 Sigma_epsilon]

## Inputs:
#
# X: n x p matrix of (stationary) time series.
# r: integer, number of factors.
# q: integer, the first q variables should not be made sparse. Default q = 0.
# alphas: vector of lasso regularisation parameters. Default is alphas = logspace(-2,3,100).
# alg: character, option for algorithm. 'PCA', '2Stage', 'EM' or 'EM-sparse'. Default is 'EM-sparse'.
# err: character, option for idiosyncratic errors. 'AR1' or 'IID'. Default is 'AR1'.
# kalman: character, option for Kalman filter and smoother equations. 'multivariate' or 'univariate'. Default is 'univariate'. 
# standardize: character, option to standardize data. TRUE or FALSE. Default is TRUE  
# max_iter: maximum number of EM iterations. Default is 100.
# threshold: tolerance on EM iterates. Default is 1e-4.



library(Matrix)

Sparse_DFM <- function(X, r, q = 0, alphas = logspace(-2,3,100), alg = 'EM-sparse', err = 'AR1', kalman = 'univariate', standardize = TRUE, max_iter=100, threshold=1e-4) {
  
  
  ## Correct input checks
  
  if(alg != 'PCA' && alg != '2Stage' && alg != 'EM' && alg != 'EM-sparse'){
    stop("Incorrect alg input")
  }
  if(err != 'IID' && err != 'AR1'){
    stop("Incorrect err input")
  }
  if(kalman != 'multivariate' && kalman != 'univariate'){
    stop("Incorrect kalman input")
  }
  if(!is.numeric(r) || r <= 0L){
    stop("r needs to be an integer > 0")
  }
  if(!is.numeric(q) || q < 0){
    stop("q needs to be an integer >= 0")
  }
  if(!is.numeric(alphas) || anyNA(alphas)){
    stop("alphas must be a numeric vector with no missing values")
  }
  if(!is.numeric(max_iter) || max_iter <= 0){
    stop("max_iter must be an integer > 0")
  }
  if(!is.numeric(threshold) || threshold <= 0){
    stop("threshold must be > 0")
  }
  
  
  # make sure lasso parameter grid is low to high 
  alphas = sort(alphas)
  
  # dimensions 
  n = dim(X)[1]
  p = dim(X)[2]
  k = r + p 
  
  # standardize if TRUE 
  X.raw = X
  X.scale = scale(X)
  X.mean = attr(X.scale, "scaled:center")
  X.sd = attr(X.scale, "scaled:scale")
  
  if(standardize){
    X = X.scale
  }
  
  
  ## Apply algorithm: PCA, 2Stage, EM, EM-sparse

  if(alg == 'PCA'){       # PCA algorithm applied 
  
    ## PCA and VAR(1) 
  
      initialise <- initPCA(X,r,err)
    
        a0_0 = initialise$a0_0
        P0_0 = initialise$P0_0
        A.tilde = initialise$A.tilde
        Lambda.tilde = initialise$Lambda.tilde
        Sigma.u.tilde = initialise$Sigma.u.tilde
        Sigma.eta = initialise$Sigma.eta
        factors.PCA = initialise$factors.pca
        loadings.PCA = initialise$loadings.pca
        
        fore_X = factors.PCA %*% t(loadings.PCA)
        if(standardize){
          fore_X = kronecker(t(X.sd),rep(1,n))*fore_X + kronecker(t(X.mean),rep(1,n))
        }
        
        errors.PCA = initialise$X.bal - factors.PCA %*% t(loadings.PCA)
        
    ## Output for PCA - depends on if err = 'AR1' or 'IID'
        
      if(err == 'AR1'){
        
        output = list(data = list(X = X.raw,
                                  standardize = standardize,
                                  X.mean = X.mean, 
                                  X.sd = X.sd,
                                  X.bal = initialise$X.bal,
                                  eigen = initialise$eigen,
                                  predict = fore_X),
                      params = list(A = A.tilde[1:r,1:r],
                                    Phi = A.tilde[(r+1):k,(r+1):k],
                                    Lambda = Lambda.tilde[,1:r],
                                    Sigma_u = Sigma.u.tilde[1:r,1:r],
                                    Sigma_epsilon = Sigma.u.tilde[(r+1):k,(r+1):k]),
                      state = list(factors = factors.PCA,
                                   errors = errors.PCA,
                                   factors.cov = P0_0[1:r,1:r],
                                   errors.cov = P0_0[(r+1):k,(r+1):k]))
        
        return(output)
        
      }else {
        
        output = list(data = list(X = X.raw,
                                  standardize = standardize,
                                  X.mean = X.mean, 
                                  X.sd = X.sd,
                                  X.bal = initialise$X.bal,
                                  eigen = initialise$eigen,
                                  predict = fore_X),
                      params = list(A = A.tilde,
                                    Lambda = Lambda.tilde,
                                    Sigma_u = Sigma.u.tilde,
                                    Sigma_epsilon = Sigma.eta),
                      state = list(factors = factors.PCA,
                                   factors.cov = P0_0[1:r,1:r]))
        
        return(output)
        
      }
        
  }else if(alg == '2Stage'){          # 2 stage algorithm applied (Doz, 2011)
    
    ## Initialise with PCA and VAR(1) 
    
      initialise <- initPCA(X,r,err)
    
        a0_0 = initialise$a0_0
        P0_0 = initialise$P0_0
        A.tilde = initialise$A.tilde
        Lambda.tilde = initialise$Lambda.tilde
        Sigma.u.tilde = initialise$Sigma.u.tilde
        Sigma.eta = initialise$Sigma.eta
        factors.PCA = initialise$factors.pca
        loadings.PCA = initialise$loadings.pca
  
    ## Kalman Filter and Smoother (Doz (2011) 2-step) 
  
      if(kalman == 'univariate'){
        KFS <- kalmanUnivariate(X, a0_0, P0_0, A.tilde, Lambda.tilde, Sigma.eta, Sigma.u.tilde)
      }else{
        KFS <- kalmanCpp(X, a0_0, P0_0, A.tilde, Lambda.tilde, Sigma.eta, Sigma.u.tilde)
      }
        
        state.KF = t(as.matrix(KFS$at_t))
        covariance.KF = KFS$Pt_t 
        state.KS = t(as.matrix(KFS$at_n))
        covariance.KS = KFS$Pt_n
        
    ## Fill in missing data in X - KF and KS 
        
        fore_X_KF = state.KF %*% t(Lambda.tilde)
        if(standardize){
          fore_X_KF = kronecker(t(X.sd),rep(1,n))*fore_X_KF + kronecker(t(X.mean),rep(1,n))
        }
        
        fore_X_KS = state.KS %*% t(Lambda.tilde)
        if(standardize){
          fore_X_KS = kronecker(t(X.sd),rep(1,n))*fore_X_KS + kronecker(t(X.mean),rep(1,n))
        }
        
        
    ## Output for 2Stage - depends on if err = 'AR1' or 'IID'
    
      if(err == 'AR1'){
        
        output = list(data = list(X = X.raw,
                                  standardize = standardize,
                                  X.mean = X.mean, 
                                  X.sd = X.sd,
                                  X.bal = initialise$X.bal,
                                  eigen = initialise$eigen,
                                  predict.KF = fore_X_KF,
                                  predict.KS = fore_X_KS),
                      params = list(A = A.tilde[1:r,1:r],
                                    Phi = A.tilde[(r+1):k,(r+1):k],
                                    Lambda = Lambda.tilde[,1:r],
                                    Sigma_u = Sigma.u.tilde[1:r,1:r],
                                    Sigma_epsilon = Sigma.u.tilde[(r+1):k,(r+1):k]),
                      state = list(factors.KF = state.KF[,1:r],
                                   errors.KF = state.KF[,(r+1):k],
                                   factors.KS = state.KS[,1:r],
                                   errors.KS = state.KS[,(r+1):k],
                                   factors.KF.cov = covariance.KF[1:r,1:r,],
                                   errors.KF.cov = covariance.KF[(r+1):k,(r+1):k,],
                                   factors.KS.cov = covariance.KS[1:r,1:r,],
                                   errors.KS.cov = covariance.KS[(r+1):k,(r+1):k,]))
        
        return(output)
        
      }else {
        
        output = list(data = list(X = X.raw,
                                  standardize = standardize,
                                  X.mean = X.mean, 
                                  X.sd = X.sd,
                                  X.bal = initialise$X.bal,
                                  eigen = initialise$eigen,
                                  predict.KF = fore_X_KF,
                                  predict.KS = fore_X_KS),
                      params = list(A = A.tilde,
                                    Lambda = Lambda.tilde,
                                    Sigma_u = Sigma.u.tilde,
                                    Sigma_epsilon = Sigma.eta),
                      state = list(factors.KF = state.KF,
                                   factors.KS = state.KS,
                                   factors.KF.cov = covariance.KF,
                                   factors.KS.cov = covariance.KS))
        
        return(output)
        
      }
        
        
      
  }else if(alg == 'EM'){             # EM algorithm applied (Banbura and Modugno, 2014)
    
    ## Initialise with PCA and VAR(1) 
    
      initialise <- initPCA(X,r,err)
    
        a0_0 = initialise$a0_0
        P0_0 = initialise$P0_0
        A.tilde = initialise$A.tilde
        Lambda.tilde = initialise$Lambda.tilde
        Sigma.u.tilde = initialise$Sigma.u.tilde
        Sigma.eta = initialise$Sigma.eta
        factors.PCA = initialise$factors.pca
        loadings.PCA = initialise$loadings.pca
        
        
    ## EM Algorithm
          
      # EM iterations function
        EM.fit <- EM(X, a0_0, P0_0, A.tilde, Lambda.tilde, Sigma.eta, Sigma.u.tilde, 
                   err = err, kalman = kalman, sparse = FALSE, max_iter = max_iter, threshold = threshold)
                   
      
      # Optimal parameters
        
        a0_0 = EM.fit$a0_0
        P0_0 = EM.fit$P0_0
        A.tilde = EM.fit$A.tilde
        Lambda.tilde = EM.fit$Lambda.tilde
        Sigma.u.tilde = EM.fit$Sigma.u.tilde
        Sigma.eta = EM.fit$Sigma.eta
        loglik.store = EM.fit$loglik.store
        converged = EM.fit$converged 
        num_iter = EM.fit$num_iter 
        
      # run KFS on final parameter estimates 
      
        if(kalman == 'univariate'){
          best.KFS <- kalmanUnivariate(X, a0_0, P0_0, A.tilde, Lambda.tilde, Sigma.eta, Sigma.u.tilde)
        }else{
          best.KFS <- kalmanCpp(X, a0_0, P0_0, A.tilde, Lambda.tilde, Sigma.eta, Sigma.u.tilde)
        }
      
      state.EM = t(best.KFS$at_n)
      covariance.EM = best.KFS$Pt_n
      
      
      # fill in missing data in X
      
      fore_X = state.EM %*% t(Lambda.tilde)
      if(standardize){
        fore_X = kronecker(t(X.sd),rep(1,n))*fore_X + kronecker(t(X.mean),rep(1,n))
      }
      
      ## Output for EM - depends on if err = 'AR1' or 'IID'
      
      if(err == 'AR1'){
        
        output = list(data = list(X = X.raw,
                                  standardize = standardize,
                                  X.mean = X.mean, 
                                  X.sd = X.sd,
                                  X.bal = initialise$X.bal,
                                  eigen = initialise$eigen,
                                  predict = fore_X),
                      params = list(A = A.tilde[1:r,1:r],
                                    Phi = A.tilde[(r+1):k,(r+1):k],
                                    Lambda = Lambda.tilde[,1:r],
                                    Sigma_u = Sigma.u.tilde[1:r,1:r],
                                    Sigma_epsilon = Sigma.u.tilde[(r+1):k,(r+1):k]),
                      state = list(factors = state.EM[,1:r],
                                   errors = state.EM[,(r+1):k],
                                   factors.cov = covariance.EM[1:r,1:r,],
                                   errors.cov = covariance.EM[(r+1):k,(r+1):k,]),
                      converged = list(converged = converged,
                                       loglik = loglik.store,
                                       num_iter = num_iter,
                                       tol = threshold,
                                       max_iter = max_iter))
        
        return(output)
        
      }else {
        
        output = list(data = list(X = X.raw,
                                  standardize = standardize,
                                  X.mean = X.mean, 
                                  X.sd = X.sd,
                                  X.bal = initialise$X.bal,
                                  eigen = initialise$eigen,
                                  predict = fore_X),
                      params = list(A = A.tilde,
                                    Lambda = Lambda.tilde,
                                    Sigma_u = Sigma.u.tilde,
                                    Sigma_epsilon = Sigma.eta),
                      state = list(factors = state.EM,
                                   factors.cov = covariance.EM),
                      converged = list(converged = converged,
                                       loglik = loglik.store,
                                       num_iter = num_iter,
                                       tol = threshold,
                                       max_iter = max_iter))
        
        return(output)
        
      }

    
    
  }else {         # sparse EM algorithm applied (Mosley et al, 2022)
    
    
    ## Initialise with PCA and VAR(1) 
    
      initialise <- initPCA(X,r,err)
      
        a0_0 = initialise$a0_0
        P0_0 = initialise$P0_0
        A.tilde = initialise$A.tilde
        Lambda.tilde = initialise$Lambda.tilde
        Sigma.u.tilde = initialise$Sigma.u.tilde
        Sigma.eta = initialise$Sigma.eta
        factors.PCA = initialise$factors.pca
        loadings.PCA = initialise$loadings.pca
    
    
    ## Sparsified EM Algorithm
      
      ## Loop over alphas and calculate BIC until column of Lambda becomes 0 
        
        bic <- c()
        num_iter = c()
        best.bic <- .Machine$double.xmax
        
        for(alphas.index in 1:length(alphas)) {
          

          # lasso regularisation parameter 
            
            alpha.value = alphas[alphas.index]
          
          # make into a matrix 
           
            alpha.value = matrix(alpha.value, nrow = p, ncol = r)
            
          # alpha parameter set to 0 for no regularisation on first q series
            
           if(q > 0){
             alpha.value[1:q,] = 0
           }
          
          # EM iterations function
            
            EM.fit <- EM(X, a0_0, P0_0, A.tilde, Lambda.tilde, Sigma.eta, Sigma.u.tilde, 
                       alpha.lasso = alpha.value, err = err, kalman = kalman, sparse = TRUE, 
                       max_iter = max_iter, threshold = threshold)
          
          # store number of iterations for each alpha
            
            num_iter[alphas.index] = EM.fit$num_iter
          
          # check if a column of Lambda has been set entirely to 0
            
            if(any(colSums(EM.fit$Lambda.tilde[(q+1):p,1:r]) == 0)){
              break
            }
          
          # Update parameters - used for warm start of the EM algorithm  
            
            a0_0 = EM.fit$a0_0
            P0_0 = EM.fit$P0_0
            A.tilde = EM.fit$A.tilde
            Lambda.tilde = EM.fit$Lambda.tilde 
            Sigma.u.tilde = EM.fit$Sigma.u.tilde
            Sigma.eta = EM.fit$Sigma.eta
          
          
          # run KFS on final parameter estimates  
            
            if(kalman == 'univariate'){
              KFS <- kalmanUnivariate(X, a0_0, P0_0, A.tilde, Lambda.tilde, Sigma.eta, Sigma.u.tilde)
            }else{
              KFS <- kalmanCpp(X, a0_0, P0_0, A.tilde, Lambda.tilde, Sigma.eta, Sigma.u.tilde)
            }
          

          # calculate BIC 
            
            bic[alphas.index] = bic_function(X, t(KFS$at_n[1:r,]), Lambda.tilde[,1:r])
          
          
          # store estimates if BIC improved 
          
            if(bic[alphas.index] < best.bic){
            
            best.EM = EM.fit
            best.KFS = KFS
            best.bic = bic[alphas.index]
            
          }
          
          
        } 
        
        
      ## Store the optimal outputs  
      
        alphas.used = alphas[1:length(bic)]
        best.alpha = alphas[which.min(bic)]
        loglik.store = best.EM$loglik.store 
        converged = best.EM$converged 

        A.tilde = best.EM$A.tilde
        Lambda.tilde = best.EM$Lambda.tilde
        Sigma.u.tilde = best.EM$Sigma.u.tilde
        Sigma.eta = best.EM$Sigma.eta
        
        state.EM = t(best.KFS$at_n)
        covariance.EM = best.KFS$Pt_n
        
        
      ## Fill in missing data in X
        
        fore_X = state.EM %*% t(Lambda.tilde)
        if(standardize){
          fore_X = kronecker(t(X.sd),rep(1,n))*fore_X + kronecker(t(X.mean),rep(1,n))
        }
        
      ## Output for EM-sparse - depends on if err = 'AR1' or 'IID'
        
        if(err == 'AR1'){
          
          output = list(data = list(X = X.raw,
                                    standardize = standardize,
                                    X.mean = X.mean, 
                                    X.sd = X.sd,
                                    X.bal = initialise$X.bal,
                                    eigen = initialise$eigen,
                                    predict = fore_X),
                        params = list(A = A.tilde[1:r,1:r],
                                      Phi = A.tilde[(r+1):k,(r+1):k],
                                      Lambda = Lambda.tilde[,1:r],
                                      Sigma_u = Sigma.u.tilde[1:r,1:r],
                                      Sigma_epsilon = Sigma.u.tilde[(r+1):k,(r+1):k]),
                        state = list(factors = state.EM[,1:r],
                                     errors = state.EM[,(r+1):k],
                                     factors.cov = covariance.EM[1:r,1:r,],
                                     errors.cov = covariance.EM[(r+1):k,(r+1):k,]),
                        converged = list(converged = converged,
                                         alpha_grid = alphas.used,
                                         alpha_opt = best.alpha,
                                         bic = bic,
                                         loglik = loglik.store,
                                         num_iter = num_iter,
                                         tol = threshold,
                                         max_iter = max_iter))
          
          return(output)
          
        }else {
          
          output = list(data = list(X = X.raw,
                                    standardize = standardize,
                                    X.mean = X.mean, 
                                    X.sd = X.sd,
                                    X.bal = initialise$X.bal,
                                    eigen = initialise$eigen,
                                    predict = fore_X),
                        params = list(A = A.tilde,
                                      Lambda = Lambda.tilde,
                                      Sigma_u = Sigma.u.tilde,
                                      Sigma_epsilon = Sigma.eta),
                        state = list(factors = state.EM,
                                     factors.cov = covariance.EM),
                        converged = list(converged = converged,
                                         alpha_grid = alphas.used,
                                         alpha_opt = best.alpha,
                                         bic = bic,
                                         loglik = loglik.store,
                                         num_iter = num_iter,
                                         tol = threshold,
                                         max_iter = max_iter))
          
          return(output)
          
        }
        
    
  }
  
}
  
