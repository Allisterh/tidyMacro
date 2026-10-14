# Independent R reference for the cluster meats in src/se_cluster.cpp:
# one-way S'S with S the within-cluster score sums, and Cameron-Gelbach-Miller
# two-way M_a + M_b - M_(a,b).
cl_meat_ref <- function(g, cl) crossprod(rowsum(g, cl))

test_that("fLPPanel time, unit and two-way cluster SEs match the reference", {
  d <- sim_panel_df()
  X <- resid(lm(cbind(shock, size) ~ factor(unit) + factor(tt), d))
  u <- resid(lm(y ~ shock + size + factor(unit) + factor(tt), d))
  g <- X * u
  B <- solve(crossprod(X))
  meats <- list(
    cl_meat_ref(g, d$tt),
    cl_meat_ref(g, d$unit),
    cl_meat_ref(g, d$unit) + cl_meat_ref(g, d$tt) -
      cl_meat_ref(g, paste(d$unit, d$tt))
  )
  cls <- list(~tt, ~unit, ~unit + tt)
  for (k in seq_along(cls)) {
    f <- fLPPanel(y ~ shock + size | unit + tt, data = d,
                  panel_id = c("unit", "tt"), shock = "shock", horizons = 0,
                  cluster = cls[[k]], n_threads = NT)
    V <- B %*% meats[[k]] %*% B
    expect_equal(f$irfs_se[1, 1], sqrt(V[1, 1]), tolerance = 1e-10)
  }
})

test_that("fLPDID CR1 SE matches a reghdfe-style reference at event time 0", {
  d <- sim_did_df()
  d <- d[order(d$id, d$year), ]
  lagv <- function(v) ave(v, d$id, FUN = function(z) c(NA, head(z, -1)))
  dD <- d$treat - lagv(d$treat)
  dy <- d$y - lagv(d$y)
  keep <- !is.na(dD) & (dD == 1 | (dD == 0 & d$treat == 0))
  s <- data.frame(dy = dy, dD = dD, year = d$year, id = d$id)[keep, ]
  s <- s[ave(s$dy, s$year, FUN = length) > 1, ]     # singleton years
  x <- resid(lm(dD ~ factor(year), s))
  e <- resid(lm(dy ~ dD + factor(year), s))
  G <- length(unique(s$id)); m <- nrow(s); nY <- length(unique(s$year))
  cadj <- G / (G - 1) * (m - 1) / (m - 1 - nY)
  se <- sqrt(cadj * sum(rowsum(x * e, s$id)^2)) / sum(x^2)

  r <- fLPDID(y ~ treat, data = d, panel_id = c("id", "year"), treat = "treat",
              post = 2, pre = 2, n_threads = NT)
  r0 <- r[r$event_time == 0, ]
  expect_equal(r0$se, se, tolerance = 1e-10)
  expect_equal(r0$nobs, m)
  expect_equal(r0$nclust, G)
})
