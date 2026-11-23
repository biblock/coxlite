# coxlite

<!-- badges: start -->
[![R-CMD-check](https://github.com/biblock/coxlite/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/biblock/coxlite/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Cox proportional hazards model written from scratch, for Biostat 625 HW4.
It fits the same models as `survival::coxph()` and gets the same answers,
but the partial likelihood, score, information and Newton-Raphson are all
done here, in C++ (with [Rcpp](https://www.rcpp.org/)) and in plain R.

## What's in it

- `cox_ph()`: formula interface, works with factors and `strata()`
- Efron or Breslow ties
- C++ (`engine = "cpp"`, default) or R (`engine = "R"`) for the likelihood
  calculations
- `print()`, `summary()`, `coef()`, `vcov()`, `confint()`, `logLik()`,
  `AIC()`, `predict()`
- `cox_basehaz()` for the baseline hazard / survival curves
- `sim_cox()` for simulating data, and an example data set `sim_trial`

## Installation

```r
# install.packages("remotes")
remotes::install_github("biblock/coxlite", build_vignettes = TRUE)
```

You need a C++ compiler (`xcode-select --install` on a Mac, or
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) on Windows).

## Example

```r
library(coxlite)

fit <- cox_ph(Surv(time, status) ~ trt + age + biomarker + sex + strata(center),
              data = sim_trial)
fit
#> Call:
#> cox_ph(formula = Surv(time, status) ~ trt + age + biomarker + 
#>     sex + strata(center), data = sim_trial)
#> 
#>                  coef exp(coef) se(coef)      z  Pr(>|z|)
#> trttreatment -0.46420    0.6286 0.126451 -3.671 2.416e-04
#> age           0.03391    1.0345 0.007449  4.552 5.305e-06
#> biomarker     0.37815    1.4596 0.064573  5.856 4.736e-09
#> sexmale       0.29949    1.3492 0.125200  2.392 1.675e-02
#> 
#> Likelihood ratio test = 65.79 on 4 df, p = 1.753e-13
#> n = 400, number of events = 264
```

So treatment cuts the hazard by about 37% (HR 0.63). Same as
`survival::coxph()`:

```r
ref <- survival::coxph(Surv(time, status) ~ trt + age + biomarker + sex +
                         strata(center), data = sim_trial)
all.equal(coef(fit), coef(ref))
#> [1] TRUE
all.equal(vcov(fit), vcov(ref), check.attributes = FALSE)
#> [1] TRUE
all.equal(fit$loglik, ref$loglik)
#> [1] TRUE
```

Other stuff:

```r
summary(fit)                           # HRs, CIs and global tests
predict(fit, newdata, type = "risk")   # relative risks for new subjects
cox_basehaz(fit)                       # cumulative baseline hazard
cox_ph(..., ties = "breslow")          # Breslow ties
cox_ph(..., engine = "R")              # pure-R engine
```

## How it works

Sorting by decreasing time means the risk set just grows as you go down the
rows, so the log partial likelihood, score and information can all be
computed in one pass. `cox_fit()` runs Newton-Raphson on that, halving the
step if the likelihood goes down. More detail in `vignette("coxlite")`.

## Correctness and speed

`vignette("benchmarks", package = "coxlite")` has the comparisons with
`coxph()` on real and simulated data, plus timings from `bench::mark()`.
Everything matches up to rounding error. The C++ version runs about as fast
as `coxph()` and is a lot faster than the R version.

## Tests

Tests use **testthat** and run on GitHub Actions:

```r
devtools::test()
```

## Help

See `?cox_ph`, `?cox_fit`, `?predict.coxlite`, `?cox_basehaz`, `?sim_cox`
and the vignettes:

```r
browseVignettes("coxlite")
```
