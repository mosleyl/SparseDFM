##  Sparsified EM algorithm for IID or AR(1) idio errors (adapted from Banbura and Modugno (2014))
##
##  E-step: Kalman filter and smoother using parameters estimated in M-step to obtain
##          state mean, covariance and lagged-covariance. Log-likelihood calculated here.
##
##  M-step: Maximisation of expected log-likelihood equations using E-step output.
##          With LASSO regularisation applied to the loadings matrix.
##
##  Convergence: Uses the log-likelihood convergence rule from from Doz (2012)

##  Input:
##
##  X: n x p, matrix of (stationary) time series
##  a0_0: 1 x k, initial state mean vector
##  P0_0: k x k, initial state covariance matrix
##  A.tilde: k x k, initial state transition matrix
##  Lambda.tilde: p x k, initial measurement matrix
##  Sigma.eta: p x p, initial measurement equation residuals covariance matrix (diagonal)
##  Sigma.u.tilde: k x k, initial state equation residuals covariance matrix
##  alpha.lasso: lasso tuning parameter value (>= 0). Default set to 0.1.
##  err: idiosyncratic error structure, 'AR1' or 'IID'. Default is 'AR1'.
##  kalman: KFS equations, 'multivariate' or 'univariate'. Default is 'univariate'.
##  sparse: if TRUE, lasso regularisation applied to the loadings. Default is TRUE.
##  max_iter: maximum number of iterations for the EM algorithm. Default is 500.
##  threshold: threshold for log-likelihood convergence check. Default is 1e-4.
##
##  NOTE: For the DFM with AR(1) errors we have that k = r + p where r is
##        the number of factors and p is the number of variables.
##        Otherwise, k = r.
##
##  Output:
##
##  a0_0: 1 x k, optimal initial state mean vector
##  P0_0: k x k, optimal initial state covariance matrix
##  A.tilde: k x k, optimal state transition matrix
##  Lambda.tilde: p x k, optimal measurement matrix
##  Sigma.eta: p x p, optimal measurement equation residuals covariance matrix (diagonal)
##  Sigma.u.tilde: k x k, optimal state equation residuals covariance matrix
##  converged : TRUE or FALSE, did the EM algorithm converge?
##  num_iter : how many iterations did the EM algorithm take to converge?
##  loglik.store: 1 x num_iter, log-likelihood for each EM iteration


## library(pracma)


EM <- function(X, a0_0, P0_0, A.tilde, Lambda.tilde, Sigma.eta, Sigma.u.tilde, alpha.lasso = 0.1, err = 'AR1', kalman = 'univariate', sparse = TRUE, max_iter=100, threshold=1e-4) {


  if(err == 'AR1'){

    ## Initialise

    n = dim(X)[1]
    p = dim(X)[2]
    k = dim(A.tilde)[1]
    r = k-p

    previous_loglik = -.Machine$double.xmax   # at least 2 iterations
    num_iter = 0        # counter begins at 0
    converged = 0       # convergence initialised to FALSE
    loglik.store = c()  # store log likelihoods

    W = 1*!is.na(X)   # data availability - n x p - 0 for missing, 1 for observed


    ## EM iterations

    while ((num_iter < max_iter) & !converged) {

      ## E-step
      if(kalman == 'univariate'){
        KFS <- SparseDFM::kalmanUnivariate(X, a0_0, P0_0, A.tilde, Lambda.tilde, Sigma.eta, Sigma.u.tilde)
      }else{
        KFS <- SparseDFM::kalmanCpp(X, a0_0, P0_0, A.tilde, Lambda.tilde, Sigma.eta, Sigma.u.tilde)
      }

      at_n = KFS$at_n             # state mean: k x n matrix (t=1,...,n)
      Pt_n = KFS$Pt_n             # state covariance: k x k x n array (t=1,...,n)
      Pt_tlag_n = KFS$Pt_tlag_n   # state covariance with lag: k x k x n array (t=1,...,n)

      loglik = KFS$logl           # log-likelihood value

      num_iter = num_iter + 1     # new iteration


      ## M-step

      ## Initial state mean and covariance update

      a0_0 = at_n[,1]
      P_F = as.matrix(Pt_n[1:r,1:r,1])
      P_E = diag(diag(as.matrix(Pt_n[(r+1):k,(r+1):k,1])))
      P0_0 = blkdiag(P_F, P_E)

      ## State transition equation parameter updates: A.tilde and Sigma.u.tilde

      ## required sums for the update equations

      E_Ft_Ftlag = at_n[1:r,2:n] %*% t(at_n[1:r,1:(n-1)]) + apply(Pt_tlag_n[1:r,1:r,],c(1,2),sum)
      E_Ftlag_Ftlag = at_n[1:r,1:(n-1)] %*% t(at_n[1:r,1:(n-1)]) + apply(Pt_n[1:r,1:r,1:(n-1)],c(1,2),sum)
      E_Ft_Ft = at_n[1:r,2:n] %*% t(at_n[1:r,2:n]) + apply(Pt_n[1:r,1:r,2:n],c(1,2),sum)
      E_Ftlag_Ft = at_n[1:r,1:(n-1)] %*% t(at_n[1:r,2:n]) + apply(Pt_tlag_n[1:r,1:r,],c(1,2),sum)

      E_et_etlag = diag(diag(at_n[(r+1):k,2:n] %*% t(at_n[(r+1):k,1:(n-1)])) + diag(apply(Pt_tlag_n[(r+1):k,(r+1):k,],c(1,2),sum)))
      E_etlag_etlag = diag(diag(at_n[(r+1):k,1:(n-1)] %*% t(at_n[(r+1):k,1:(n-1)])) + diag(apply(Pt_n[(r+1):k,(r+1):k,1:(n-1)],c(1,2),sum)))
      E_et_et = diag(diag(at_n[(r+1):k,2:n] %*% t(at_n[(r+1):k,2:n])) + diag(apply(Pt_n[(r+1):k,(r+1):k,2:n],c(1,2),sum)))
      E_etlag_et = diag(diag(at_n[(r+1):k,1:(n-1)] %*% t(at_n[(r+1):k,2:n])) + diag(apply(Pt_tlag_n[(r+1):k,(r+1):k,],c(1,2),sum)))

      ## A update

      A = E_Ft_Ftlag %*% solve(E_Ftlag_Ftlag)

      ## Sig_u update

      Sig_u = (E_Ft_Ft - A %*% E_Ftlag_Ft)/n

      ## Phi update

      Phi = E_et_etlag %*% diag(1 / diag(E_etlag_etlag))

      ## Sig_epsilon update

      Sig_epsilon = (E_et_et - Phi %*% t(E_etlag_et)) / n

      ## A.tilde construction

      A.tilde = blkdiag(A, Phi)

      ## Sigma.u.tilde construction

      Sigma.u.tilde = blkdiag(Sig_u, Sig_epsilon)


      ## Measurement equation parameter updates: Lambda.tilde and Sigma.eta

      ## Sigma.eta update (just the same)

      kappa = 1e-4
      Sigma.eta = kappa*diag(p)

      ## Lambda.tilde update (ADMM algorithm for lasso regularisation - see paper)

      y = t(X)
      y[is.na(y)] = 0

      ## Calculate components A, B and C in the fast algorithm for Lambda computation (see paper)
      ##
      ## At: r x r x n array
      ## Bt: p x p matrix
      ## Ct: p x r matrix

      At = array(NA, dim=c(r, r, n))
      Bt = t(W)
      Ct = 0

      for(t in 1:n){

        At[,,t] = at_n[1:r,t] %*% t(at_n[1:r,t]) + Pt_n[1:r,1:r,t]
        Ct = Ct + (as.matrix(Bt[,t]*y[,t]) %*% at_n[1:r,t]) - (diag(Bt[,t]) %*% (at_n[(r+1):k,t] %*% t(at_n[1:r,t]) + Pt_n[(r+1):k,1:r,t]))

      }

      ## Estimating loadings: with or without sparsity

      if(sparse){

        ## Solve fast algorithm for Lambda

        D_cube = SparseDFM::solveCube(At, Bt, nu = kappa)
        sol = solveLambda(D_cube, Ct, nu = kappa, alpha = alpha.lasso)
        Lambda = sol$Z

        ## Lambda.tilde construction

        Lambda.tilde = cbind(Lambda, diag(p))

      }else{

        D_cube = SparseDFM::solveCube(At, Bt, nu = 0)
        Lambda = SparseDFM::fastLambda(D_cube, Ct)

        Lambda.tilde = cbind(Lambda, diag(p))

      }


      ## Check convergence

      converged <- emConverged(loglik, previous_loglik, threshold)
      previous_loglik <- loglik
      loglik.store[num_iter] = loglik


    }


    output <- list('a0_0' = a0_0, 'P0_0' = P0_0, 'A.tilde' = A.tilde, 'Sigma.u.tilde' = Sigma.u.tilde,
                   'Lambda.tilde' = Lambda.tilde, 'Sigma.eta' = Sigma.eta, 'loglik.store' = loglik.store,
                   'converged' = converged, 'num_iter' = num_iter)

    return(output)

  }else {

    ## Initialise

    n = dim(X)[1]
    p = dim(X)[2]
    k = dim(A.tilde)[1]
    r = k

    previous_loglik = -.Machine$double.xmax   # at least 2 iterations
    num_iter = 0        # counter begins at 0
    converged = 0       # convergence initialised to FALSE
    loglik.store = c()  # store log likelihoods

    W = 1*!is.na(X)   # data availability - n x p - 0 for missing, 1 for observed


    ## EM iterations

    while ((num_iter < max_iter) & !converged) {

      ## E-step
      if(kalman == 'univariate'){
        KFS <- SparseDFM::kalmanUnivariate(X, a0_0, P0_0, A.tilde, Lambda.tilde, Sigma.eta, Sigma.u.tilde)
      }else{
        KFS <- SparseDFM::kalmanCpp(X, a0_0, P0_0, A.tilde, Lambda.tilde, Sigma.eta, Sigma.u.tilde)
      }

      at_n = KFS$at_n             # state mean: k x n matrix (t=1,...,n)
      Pt_n = KFS$Pt_n             # state covariance: k x k x n array (t=1,...,n)
      Pt_tlag_n = KFS$Pt_tlag_n   # state covariance with lag: k x k x n array (t=1,...,n)

      loglik = KFS$logl           # log-likelihood value


      num_iter = num_iter + 1     # new iteration

      ## M-step

      ## Initial state mean and covariance update

      a0_0 = at_n[,1]
      P0_0 = as.matrix(Pt_n[,,1])

      ## State transition equation parameter updates: A.tilde and Sigma.u.tilde

      ## required sums for the update equations

      E_Ft_Ftlag = at_n[1:r,2:n] %*% t(at_n[1:r,1:(n-1)]) + apply(Pt_tlag_n[1:r,1:r,],c(1,2),sum)
      E_Ftlag_Ftlag = at_n[1:r,1:(n-1)] %*% t(at_n[1:r,1:(n-1)]) + apply(Pt_n[1:r,1:r,1:(n-1)],c(1,2),sum)
      E_Ft_Ft = at_n[1:r,2:n] %*% t(at_n[1:r,2:n]) + apply(Pt_n[1:r,1:r,2:n],c(1,2),sum)
      E_Ftlag_Ft = at_n[1:r,1:(n-1)] %*% t(at_n[1:r,2:n]) + apply(Pt_tlag_n[1:r,1:r,],c(1,2),sum)

      ## A update

      A = E_Ft_Ftlag %*% solve(E_Ftlag_Ftlag)
      A.tilde = A

      ## Sig_u update

      Sig_u = (E_Ft_Ft - A %*% E_Ftlag_Ft)/n
      Sig.u.tilde = Sig_u


      ## Measurement equation parameter updates: Lambda.tilde and Sigma.eta


      ## Lambda.tilde update (ADMM algorithm for lasso regularisation - see paper)

      y = t(X)
      y[is.na(y)] = 0

      ## Calculate components A, B and C in the fast algorithm for Lambda computation (see paper)
      ##
      ## At: r x r x n array
      ## Bt: p x p matrix
      ## Ct: p x r matrix

      At = array(NA, dim=c(r, r, n))
      Bt = matrix(NA, nrow = p, ncol = n)
      Ct = 0
      Sigma.eta.inv = 1/diag(Sigma.eta)

      for(t in 1:n){

        At[,,t] = at_n[1:r,t] %*% t(at_n[1:r,t]) + Pt_n[1:r,1:r,t]
        Bt[,t] = W[t,]*Sigma.eta.inv*W[t,]
        Ct = Ct + as.matrix(Bt[,t]*y[,t]) %*% at_n[1:r,t]

      }

      ## Estimating loadings: with or without sparsity

      if(sparse){

        ## Solve fast algorithm for Lambda

        D_cube = SparseDFM::solveCube(At, Bt, nu = 1)
        sol = solveLambda(D_cube, Ct, nu = 1, alpha = alpha.lasso)
        Lambda = sol$Z

        ## Lambda.tilde construction

        Lambda.tilde = Lambda

      }else{

        D_cube = SparseDFM::solveCube(At, Bt, nu = 0)
        Lambda = SparseDFM::fastLambda(D_cube, Ct)

        Lambda.tilde = Lambda

      }



      ## Sigma.eta update

      Sig_e_new = 0
      I = rep(1,p)
      for(t in 1:n){

        Sig_e_new = Sig_e_new + (tcrossprod(W[t,]*y[,t]) - tcrossprod(W[t,]*y[,t],W[t,]*Lambda%*%at_n[1:r,t])
                                 - tcrossprod(W[t,]*Lambda%*%at_n[1:r,t],W[t,]*y[,t]) + tcrossprod(W[t,]*Lambda%*%(At[,,t]),W[t,]*Lambda)
                                 + diag((I - W[t,])*diag(Sigma.eta)*(I - W[t,])))
      }
      Sig_e = Sig_e_new/n
      Sig_e = diag(diag(Sig_e))

      Sigma.eta = Sig_e

      if(any(diag(Sigma.eta)<0.001)){
        converged = TRUE
        previous_loglik <- loglik
        loglik.store[num_iter] = loglik
      }else{
        ## Check convergence

        converged <- emConverged(loglik, previous_loglik, threshold)
        previous_loglik <- loglik
        loglik.store[num_iter] = loglik
      }

    }

    output <- list('a0_0' = a0_0, 'P0_0' = P0_0, 'A.tilde' = A.tilde, 'Sigma.u.tilde' = Sigma.u.tilde,
                   'Lambda.tilde' = Lambda.tilde, 'Sigma.eta' = Sigma.eta, 'loglik.store' = loglik.store,
                   'converged' = converged, 'num_iter' = num_iter)

    return(output)


  }

}

