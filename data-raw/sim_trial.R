# makes data/sim_trial.rda (run from the package root)

set.seed(625)
n <- 400

trt <- factor(rep(c("control", "treatment"), each = n / 2),
              levels = c("control", "treatment"))
age <- round(stats::rnorm(n, mean = 62, sd = 9))
biomarker <- round(stats::rnorm(n), 2)
sex <- factor(sample(c("female", "male"), n, replace = TRUE))
center <- factor(sample(LETTERS[1:4], n, replace = TRUE))

# true log hazard ratios
lp <- -0.5 * (trt == "treatment") + 0.03 * (age - 62) + 0.4 * biomarker +
  0.2 * (sex == "male")

# weibull baseline, different scale for each centre
shape <- 1.3
scale <- c(A = 0.010, B = 0.015, C = 0.020, D = 0.030)[as.character(center)]
event_time <- (-log(stats::runif(n)) / (scale * exp(lp)))^(1 / shape)

# accrual over 24 months, analysis at 48 months, plus some dropout
cens_time <- pmin(48 - stats::runif(n, 0, 24), stats::rexp(n, 0.01))

time <- pmin(event_time, cens_time)
status <- as.integer(event_time <= cens_time)
time <- ceiling(time * 2) / 2   # half months, gives ties

sim_trial <- data.frame(time, status, trt, age, biomarker, sex, center)

save(sim_trial, file = "data/sim_trial.rda", compress = "xz")
