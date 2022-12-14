#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List kalmanUnivariate(const mat& X, const mat& a0_0, const mat& P0_0,
                      const mat& A, const mat& Lambda, const mat& Sig_e,
                      const mat& Sig_u) {

  //  Univariate treatment of multivariate series
  //  Kalman filter and smoother equations from Durbin and Koopman (2012)
  //  Smoothed state covariance with lag from de Jong and Mackinnon (1988)

  //  Works for a general state space model of the form:
  //
  //  X_t = Lambda*F_t + e_t,  e_t ~ N(0,Sig_e)
  //  F_t = A*F_{t-1} + u_t,   u_t ~ N(0,Sig_u)

  //  Inputs:
  //
  //  X: n x p, matrix of (stationary) time series
  //  a0_0: k x 1, initial state mean vector at t=0
  //  P0_0: k x k, initial state covariance matrix at t=0
  //  A: k x k, state matrix
  //  Lambda: p x k, measurement matrix
  //  Sig_e: p x p, measurement equation residuals covariance matrix (diagonal)
  //  Sig_u: k x k, state equation residuals covariance matrix
  //
  //  Outputs:
  //
  //  logl: log-likelihood required for convergence check in EM
  //  at_tlag: k x n, predicted state mean vectors
  //  Pt_tlag: k x k x n, predicted state covariance matrices
  //  at_t: k x n, filtered state mean vectors
  //  Pt_t: k x k x n, filtered state covariance matrices
  //  at_n: k x n, smoothed state mean vectors
  //  Pt_n: k x k x n, smoothed state covariance matrices
  //  Pt_tlag_n: k x k x n, smoothed state covariance with lag
  //
  //  NOTE: For the DFM with AR(1) errors we have that k = r + p where r is
  //        the number of factors and p is the number of variables.

  // Initialise
  const uword n = X.n_rows;
  const uword p = X.n_cols;
  const uword k = A.n_rows;

  mat at_tlag(k, n, fill::none);        // predicted state mean
  cube Pt_tlag(k, k, n, fill::none);    // predicted state covariance
  mat at_t(k, n, fill::none);           // filtered state mean
  cube Pt_t(k, k, n, fill::none);       // filtered state covariance
  mat at_n(k, n, fill::none);           // smoothed state mean
  cube Pt_n(k, k, n, fill::none);       // smoothed state covariance
  cube Pt_tlag_n(k, k, n, fill::none);  // smoothed state covariance with lag

  vec at_i = vectorise(a0_0);           // initial state mean
  mat Pt_i = P0_0;                      // initial state covariance

  double logl = 0.0;                    // log-likelihood

  mat vt(p, n, fill::none);             // innovation
  mat inv_Ft(p, n, fill::zeros);        // inverse innovation variance
  cube Kt(k, p, n, fill::none);         // Kalman gain

  // Kalman filter loop for t = 1,...,n
  const vec diag_Sig_e = diagvec(Sig_e);
  const double log2pi = log(2.0 * datum::pi);

  for (uword t = 0; t < n; ++t) {
    // Prediction equations
    at_i = A * at_i;
    Pt_i = A * Pt_i * A.t() + Sig_u;
    at_tlag.col(t) = at_i;
    Pt_tlag.slice(t) = Pt_i;

    // Update equations
    for (uword i = 0; i < p; ++i) {
      const double yt_i = X(t, i);      // y = t(X)
      if (std::isnan(yt_i)) continue;   // omit NA

      const rowvec Zt_i = Lambda.row(i);
      const double vt_i = yt_i - dot(Zt_i, at_i);
      const vec PZt_i = Pt_i * Zt_i.t();
      const double Ft_i = dot(Zt_i, PZt_i) + diag_Sig_e(i);

      // Skip yt_i if Ft_i is zero
      if (abs(Ft_i) <= datum::eps) continue;

      const double inv_Ft_i = 1.0 / Ft_i;
      const vec Kt_i = PZt_i * inv_Ft_i;
      at_i += Kt_i * vt_i;
      Pt_i -= Kt_i * PZt_i.t();

      // Enforce symmetry
      Pt_i = 0.5 * (Pt_i + Pt_i.t());

      // Calculate log-likelihood
      logl -= 0.5 * (log2pi + log(Ft_i) + vt_i * inv_Ft_i * vt_i);

      vt(i, t) = vt_i;
      inv_Ft(i, t) = inv_Ft_i;
      Kt.slice(t).col(i) = Kt_i;
    }
    at_t.col(t) = at_i;
    Pt_t.slice(t) = Pt_i;
  }

  // Initialise Kalman smoother with zeros
  vec rt_i(k, fill::zeros);
  mat Nt_i(k, k, fill::zeros);

  // Kalman smoother loop for t = n,...,1
  const mat Ik = eye(k, k);
  for (uword t = n - 1; t != uword(-1); --t) {
    for (uword i = p - 1; i != uword(-1); --i) {
      const rowvec Zt_i = Lambda.row(i);
      const double inv_Ft_i = inv_Ft(i, t);
      const vec Kt_i = Kt.slice(t).col(i);
      const mat Lt_i = Ik - Kt_i * Zt_i;
      if (inv_Ft_i) {
        // yt_i is not missing and Ft_i is not zero
        const vec ZFt_i = Zt_i.t() * inv_Ft_i;
        rt_i = ZFt_i * vt(i, t) + Lt_i.t() * rt_i;
        Nt_i = ZFt_i * Zt_i + Lt_i.t() * Nt_i * Lt_i;
      } else {
        rt_i = Lt_i.t() * rt_i;
        Nt_i = Lt_i.t() * Nt_i * Lt_i;
      }
    }
    const mat Pt_1 = Pt_tlag.slice(t);
    at_n.col(t) = at_tlag.col(t) + Pt_1 * rt_i;
    Pt_n.slice(t) = Pt_1 - Pt_1 * Nt_i * Pt_1;
    if (t > 0) {
      rt_i = A.t() * rt_i;
      Nt_i = A.t() * Nt_i * A;
    }
  }

  // Calculate smoothed state covariance with lag
  Pt_tlag_n.slice(0) = Pt_n.slice(0) * solve(Pt_tlag.slice(0), A) * P0_0;
  for (uword t = 1; t < n; ++t) {
    Pt_tlag_n.slice(t) = Pt_n.slice(t) * solve(Pt_tlag.slice(t), A) * Pt_t.slice(t - 1);
  }

  // Output
  return List::create(Named("logl")       = logl,
                      Named("at_t")       = at_t,
                      Named("Pt_t")       = Pt_t,
                      Named("at_n")       = at_n,
                      Named("Pt_n")       = Pt_n,
                      Named("Pt_tlag_n")  = Pt_tlag_n);
}
