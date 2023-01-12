## Function to get initial parameter estimates for the EM algorithm.
##
## Inputs:
#
# X: n x p matrix of (stationary) time series
# r: number of factors
# err: idiosyncratic error structure, 'AR1' or 'IID'. Default is 'AR1'.


## library(pracma)

initPCA <- function(X,r,err='AR1') {

  # interpolate the missing data in X
  fillx = fillNA(X)
  X.balanced<- fillx$X
  X.NAidx <- fillx$idx.na

  n <- dim(X.balanced)[1] # n observations
  p <- dim(X.balanced)[2] # p time series


  ## PCA on balanced panel

  evd = eigen(cov(X.balanced))
  loadings.pca = as.matrix(evd$vectors[,1:r]) # (p x r)
  loadings.pca[,which(colSums(loadings.pca)<0)] = -loadings.pca[,which(colSums(loadings.pca)<0)] # flip sign
  factors.pca = X.balanced %*% loadings.pca # (n x r)
  factors = factors.pca


  ## Two error structures: err = 'AR1' or err = 'IID'

  if(err == 'AR1'){

    ## Measurement equation parameters
    #
    ## Lambda.tilde = [Lambda I_p] is p x (r+p)
    ## Sigma.eta = \kappa*I_p where \kappa is a small number

      Lambda.tilde = cbind(loadings.pca, diag(p))
      kappa = 1e-4
      Sigma.eta = kappa*diag(p)


    ## State equation parameters
    #
    ## A.tilde = [A 0; 0 Phi] is (r+p) x (r+p)
    ## Sigma.u.tilde = [Sigma_u 0; 0 Sigma_{epsilon}] is (r+p) x (r+p)

      VAR_fit =  VAR(factors, 1)
      A = t(VAR_fit$A)
      Sigma_u = cov(VAR_fit$res)

      e = X.balanced - factors %*% t(loadings.pca)
      e[X.NAidx] <- NA # do not consider missing data when working out these parameters

      phi = sigma = c()
      for(i in 1:p){
        ei = na.omit(e[,i])
        n_ei = length(ei)
        phi[i] = solve(ei[1:(n_ei-1)]%*%ei[1:(n_ei-1)])%*%ei[1:(n_ei-1)]%*%ei[2:(n_ei)]
        sigma[i] = var(ei[2:(n_ei)]-phi[i]*ei[1:(n_ei-1)])
      }

      Phi = diag(phi)
      Sigma_epsilon = diag(sigma)

      A.tilde = blkdiag(A, Phi)

      Sigma.u.tilde = blkdiag(Sigma_u, Sigma_epsilon)


    ## Initial mean and variance of state
    #
    # a0_0 is (r+p)x1 mean of state at t=0
    # P0_0 is (r+p)x(r+p) variance of state at t=0

      a0_0 = as.matrix(rep(0, r+p))
      P_F = matrix(corpcor::pseudoinverse(diag(r*r)- kronecker(A,A)) %*% matrix(Sigma_u, ncol = 1), r, r)
      P_eps = diag(1 / diag(diag(dim(Phi)[1]) - Phi ^ 2)) * Sigma_epsilon
      P0_0 = blkdiag(P_F, P_eps) # covariance of F and e initialised to be 0
  }else { # IID errors

    ## Measurement equation parameters: Lambda and Sig_e

      Lambda.tilde = loadings.pca

      e = X.balanced - factors %*% t(loadings.pca)
      e[X.NAidx] <- NA # do not consider missing data when working out these parameters

      sigma = c()
      for(i in 1:p){
        ei = na.omit(e[,i])
        sigma[i] = var(ei)
      }

      Sigma.eta = diag(sigma)

    ## State equation parameters: A and Sig_u

      VAR_fit =  VAR(factors, 1)
      A.tilde = t(VAR_fit$A)
      Sigma.u.tilde = cov(VAR_fit$res)

    ## Initial mean and variance of state
    #
    # a0_0 is rx1 mean of state at t=0
    # P0_0 is rxr variance of state at t=0

      a0_0 = as.matrix(rep(0, r))
      P0_0 = matrix(corpcor::pseudoinverse(diag(r*r)- kronecker(A.tilde,A.tilde)) %*% matrix(Sigma.u.tilde, ncol = 1), r, r)

  }


  output = list('a0_0' = a0_0, 'P0_0' = P0_0, 'factors.pca' = factors.pca,
                'loadings.pca' = loadings.pca, 'A.tilde' = A.tilde, 'Sigma.u.tilde' = Sigma.u.tilde,
                'Lambda.tilde' = Lambda.tilde, 'Sigma.eta' = Sigma.eta, 'X.bal' = X.balanced,
                'eigen' = evd)

  return(output)

}
