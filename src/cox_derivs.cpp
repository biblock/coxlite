// log partial likelihood, score and information (C++ version of engine_r.R)
// rows sorted by stratum then decreasing time, so the risk set is just every
// row seen so far in the stratum
//
// S0, S1, S2: sums of r, r*x, r*x*x' over the risk set, r = exp(x'beta)
// D0, D1, D2: same sums over the events at the current time
// efron uses S - (l/d) D for the l-th of d ties, breslow just uses S

#include <Rcpp.h>
#include <cmath>
#include <vector>

using namespace Rcpp;

// [[Rcpp::export]]
List cox_derivs_cpp(const NumericVector& time,
                    const IntegerVector& status,
                    const NumericMatrix& x,
                    const IntegerVector& strata,
                    const NumericVector& beta,
                    const bool efron) {
  const int n = x.nrow();
  const int p = x.ncol();

  std::vector<double> eta(n, 0.0), risk(n);
  for (int k = 0; k < p; ++k)
    for (int i = 0; i < n; ++i) eta[i] += x(i, k) * beta[k];
  for (int i = 0; i < n; ++i) risk[i] = std::exp(eta[i]);

  double loglik = 0.0;
  NumericVector score(p);
  NumericMatrix info(p, p);

  // p x p sums stored flat, only the upper triangle
  double s0 = 0.0, d0 = 0.0;
  std::vector<double> s1(p), s2(p * p), d1(p), d2(p * p), a1(p);

  int i = 0;
  while (i < n) {
    // reset at each new stratum
    if (i == 0 || strata[i] != strata[i - 1]) {
      s0 = 0.0;
      std::fill(s1.begin(), s1.end(), 0.0);
      std::fill(s2.begin(), s2.end(), 0.0);
    }

    // add rows i..j-1 (tied times) to the sums
    d0 = 0.0;
    std::fill(d1.begin(), d1.end(), 0.0);
    std::fill(d2.begin(), d2.end(), 0.0);
    int ndead = 0;
    int j = i;
    while (j < n && strata[j] == strata[i] && time[j] == time[i]) {
      const double r = risk[j];
      const bool dead = status[j] == 1;
      s0 += r;
      if (dead) {
        ++ndead;
        loglik += eta[j];
        d0 += r;
      }
      for (int k = 0; k < p; ++k) {
        const double rxk = r * x(j, k);
        s1[k] += rxk;
        if (dead) {
          score[k] += x(j, k);
          d1[k] += rxk;
        }
        for (int m = k; m < p; ++m) {
          s2[k * p + m] += rxk * x(j, m);
          if (dead) d2[k * p + m] += rxk * x(j, m);
        }
      }
      ++j;
    }

    // event contributions
    for (int l = 0; l < ndead; ++l) {
      const double frac = efron ? static_cast<double>(l) / ndead : 0.0;
      const double den = s0 - frac * d0;
      loglik -= std::log(den);
      for (int k = 0; k < p; ++k) {
        a1[k] = (s1[k] - frac * d1[k]) / den;
        score[k] -= a1[k];
      }
      for (int k = 0; k < p; ++k) {
        for (int m = k; m < p; ++m) {
          const double a2 = (s2[k * p + m] - frac * d2[k * p + m]) / den;
          info(k, m) += a2 - a1[k] * a1[m];
        }
      }
    }
    i = j;
  }

  // fill in lower triangle
  for (int k = 0; k < p; ++k)
    for (int m = 0; m < k; ++m) info(k, m) = info(m, k);

  return List::create(Named("loglik") = loglik,
                      Named("score") = score,
                      Named("info") = info);
}
