# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016 Trustees of Columbia University
# Copyright (C) 2017 Sam Brilleman
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

# tests can be run using devtools::test() or manually by loading testthat 
# package and then running the code below possibly with options(mc.cores = 4).

library(rstanarm)
stopifnot(require(lme4))
stopifnot(require(survival))
ITER <- 1000
CHAINS <- if (interactive()) 1 else 1
SEED <- 12345
REFRESH <- ITER
set.seed(SEED)
if (interactive()) options(mc.cores = parallel::detectCores())
TOLSCALE_lmer  <- 0.2 # how many SEs can stan_jm fixefs be from lmer fixefs
TOLSCALE_glmer <- 0.3 # how many SEs can stan_jm fixefs be from glmer fixefs
TOLSCALE_event <- 0.2 # how many SEs can stan_jm fixefs be from coxph fixefs
TOLSCALE_jms <- 0.2 # how many SEs can stan_jm fixefs/ses be from JM fixefs
FIXEF_tol <- 0.02
RANEF_tol <- 0.05
EVENT_tol <- 0.05
PP_tol <- 0.3
PPMEANS_tol <- 0.05

expect_stanreg <- function(x) expect_s3_class(x, "stanreg")
expect_stanjm  <- function(x) expect_s3_class(x, "stanjm")
SW <- function(expr) capture.output(suppressWarnings(expr))
colMedians <- function(x) apply(x, 2L, median)
colSds <- function(x) apply(x, 2L, sd)
colLower <- function(x) apply(x, 2L, quantile, 0.025)
colUpper <- function(x) apply(x, 2L, quantile, 0.975)

context("stan_jm")

#--------  Models

examplejm1 <- 
  stan_jm(
    formulaLong = logBili ~ year + (1 | id), 
    dataLong = pbcLong,
    formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
    dataEvent = pbcSurv,
    time_var = "year", iter = 1,
    chains = CHAINS, seed = SEED)

examplejm2 <- 
  stan_jm(
    formulaLong = logBili ~ year + (year | id), 
    dataLong = pbcLong,
    formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
    dataEvent = pbcSurv,
    time_var = "year", iter = 1,
    chains = CHAINS, seed = SEED)

examplejm3 <- 
  stan_jm(
    formulaLong = list(
      logBili ~ year + (year | id),
      albumin ~ year + (1 | id)),
    dataLong = pbcLong,
    formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
    dataEvent = pbcSurv,
    time_var = "year", iter = 1,
    chains = CHAINS, seed = SEED)

examplejm4 <- 
  stan_jm(
    formulaLong = list(
      logBili ~ year + (year | id),
      albumin ~ year + (year | id)),
    dataLong = pbcLong,
    formulaEvent = Surv(futimeYears, death) ~ sex + trt, 
    dataEvent = pbcSurv,
    time_var = "year", iter = 1,
    chains = CHAINS, seed = SEED)

#--------  Tests

fmL1 <- logBili ~ year + (1 | id)
fmL2 <- exp(logBili) ~ year + (1 | id)
fmL3 <- logBili ~ I(year / 2) + (1 | id)
fmL4 <- logBili ~ poly(year, degree = 2) + (1 | id)
fmS1  <- Surv(futimeYears, death) ~ sex + trt
fmS2  <- Surv(futimeYears, death) ~ sex + I(trt / 2)
test_that("formula argument works, including evaluation of transformations", {
  expect_stanjm(stan_jm(logBili ~ year + (1 | id), pbcLong, 
                        Surv(futimeYears, death) ~ sex + trt, pbcSurv, 
                        time_var = "year", iter = 1, chains = CHAINS, seed = SEED))
  expect_stanjm(stan_jm(fmL1, pbcLong, fmS1, pbcSurv, time_var = "year",
                        iter = 1, chains = CHAINS, seed = SEED)) 
  expect_stanjm(stan_jm(fmL2, pbcLong, fmS1, pbcSurv, time_var = "year",
                        iter = 1, chains = CHAINS, seed = SEED)) 
  expect_stanjm(stan_jm(fmL3, pbcLong, fmS1, pbcSurv, time_var = "year",
                        iter = 1, chains = CHAINS, seed = SEED)) 
  expect_stanjm(stan_jm(fmL4, pbcLong, fmS1, pbcSurv, time_var = "year",
                        iter = 1, chains = CHAINS, seed = SEED))   
  expect_stanjm(stan_jm(fmL1, pbcLong, fmS2, pbcSurv, time_var = "year",
                        iter = 1, chains = CHAINS, seed = SEED))  
})


test_that("data argument works", {
  
})


test_that("id_var argument works", {
  
  # Models with a single grouping factor
  expect_output(update(examplejm1, id_var = "id"))
  expect_output(expect_warning(update(examplejm1, id_var = "year"), 
                               "are not the same; 'id_var' will be ignored"))
  
  # Models with more than one grouping factor
  tmpdat <- pbcLong
  tmpdat$practice <- cut(pbcLong$id, c(0,10,20,30,40))
  tmpfm <- logBili ~ year + (1 | id) + (1 | practice)
  expect_error(update(examplejm1, formulaLong. = tmpfm, dataLong = tmpdat), 
               "'id_var' must be specified")
  expect_error(update(examplejm1, formulaLong. = tmpfm, dataLong = tmpdat, id_var = "year"), 
               "'id_var' must be included as a grouping factor")
  expect_output(update(examplejm1, formulaLong. = tmpfm, dataLong = tmpdat, id_var = "id", init = 0))
})


test_that("family argument works", {
  # not checked here, since can be checked when comparing stan_jm and stan_glmer below
})


test_that("assoc argument works", {
  
  # Univariate joint models
  expect_output(ret <- update(examplejm2, assoc = NULL))
  expect_output(ret <- update(examplejm2, assoc = "null"))
  expect_output(ret <- update(examplejm2, assoc = "etavalue"))
  expect_output(ret <- update(examplejm2, assoc = "etaslope"))
  expect_output(ret <- update(examplejm2, assoc = "muvalue"))
  expect_output(ret <- update(examplejm2, assoc = "muslope"))
  expect_output(ret <- update(examplejm2, assoc = c("etavalue", "etaslope"))) 
  expect_output(ret <- update(examplejm2, assoc = c("etavalue", "muslope"))) 
  expect_output(ret <- update(examplejm2, assoc = c("muvalue", "etaslope"))) 
  expect_output(ret <- update(examplejm2, assoc = c("muvalue", "muslope")))
  
  expect_error(ret <- update(examplejm2, assoc = c("etavalue", "muvalue")), "cannot be specified together")
  expect_error(ret <- update(examplejm2, assoc = c("etaslope", "muslope")), "cannot be specified together")
  
  expect_output(ret <- update(examplejm2, assoc = "shared_b"))
  expect_output(ret <- update(examplejm2, assoc = "shared_b(1)"))
  expect_output(ret <- update(examplejm2, assoc = "shared_b(2)"))
  expect_output(ret <- update(examplejm2, assoc = "shared_b(1:2)"))
  expect_output(ret <- update(examplejm2, assoc = "shared_b(1,2)"))
  
  expect_output(ret <- update(examplejm2, assoc = "shared_coef"))
  expect_output(ret <- update(examplejm2, assoc = "shared_coef(1)"))
  expect_output(ret <- update(examplejm2, assoc = "shared_coef(2)"))
  expect_output(ret <- update(examplejm2, assoc = "shared_coef(1:2)"))
  expect_output(ret <- update(examplejm2, assoc = "shared_coef(1,2)"))
  
  expect_error(ret <- update(examplejm2, assoc = "shared_b(10)"), "greater than the number of")
  expect_error(ret <- update(examplejm2, assoc = "shared_coef(10)"), "greater than the number of")
  expect_error(ret <- update(examplejm2, assoc = c("shared_b(1)", "shared_coef(1)")), "should not be specified in both")
  expect_error(ret <- update(examplejm2, assoc = c("shared_b", "shared_coef")), "should not be specified in both")
  
  expect_output(ret <- update(examplejm2, assoc = list(NULL)))
  expect_output(ret <- update(examplejm2, assoc = list("null")))
  expect_output(ret <- update(examplejm2, assoc = list("etavalue")))
  expect_output(ret <- update(examplejm2, assoc = list("etaslope")))
  expect_output(ret <- update(examplejm2, assoc = list("muvalue")))
  expect_output(ret <- update(examplejm2, assoc = list("muslope")))
  expect_output(ret <- update(examplejm2, assoc = list(c("etavalue", "etaslope")))) 
  expect_output(ret <- update(examplejm2, assoc = list(c("etavalue", "muslope")))) 
  expect_output(ret <- update(examplejm2, assoc = list(c("muvalue", "etaslope")))) 
  expect_output(ret <- update(examplejm2, assoc = list(c("muvalue", "muslope"))))  
  
  expect_error(ret <- update(examplejm2, assoc = NA), "'assoc' should be") 
  expect_error(ret <- update(examplejm2, assoc = 123), "'assoc' should be") 
  expect_error(ret <- update(examplejm2, assoc = c(1,2,3)), "'assoc' should be") 
  
  expect_error(ret <- update(examplejm2, assoc = c("wrong")), "unsupported association type") 
  expect_error(ret <- update(examplejm2, assoc = list("wrong")), "unsupported association type") 
  
  expect_error(ret <- update(examplejm2, assoc = list(NULL, NULL)), "length equal to the number of") 
  expect_error(ret <- update(examplejm2, assoc = list("etavalue", "etavalue")), "length equal to the number of") 
  expect_error(ret <- update(examplejm2, assoc = list(c("etavalue", "etaslope"), "etavalue")), "length equal to the number of") 
  
  # Multivariate joint models
  
  expect_output(ret <- update(examplejm3, assoc = "etavalue"))
  expect_output(ret <- update(examplejm3, assoc = "etaslope"))
  expect_output(ret <- update(examplejm3, assoc = "muvalue"))
  expect_output(ret <- update(examplejm3, assoc = "muslope"))
  expect_output(ret <- update(examplejm3, assoc = c("etavalue", "etaslope"))) 
  expect_output(ret <- update(examplejm3, assoc = c("etavalue", "muslope"))) 
  expect_output(ret <- update(examplejm3, assoc = c("muvalue", "etaslope"))) 
  expect_output(ret <- update(examplejm3, assoc = c("muvalue", "muslope")))
  
  expect_output(ret <- update(examplejm3, assoc = list("etavalue")))
  expect_output(ret <- update(examplejm3, assoc = list("etavalue", "etavalue")))
  expect_output(ret <- update(examplejm3, assoc = list(c("etavalue", "etaslope"), "etavalue")))
  expect_output(ret <- update(examplejm3, assoc = list("etavalue", c("etavalue", "etaslope"))))
  expect_output(ret <- update(examplejm3, assoc = list(c("etavalue", "etaslope"), c("muvalue", "muslope"))))
  
  expect_error(ret <- update(examplejm3, assoc = list("wrong", "etavalue")), "unsupported association type")
  expect_error(ret <- update(examplejm3, assoc = list("null", "etavalue", "etaslope")), "length equal to the number of")
  expect_error(ret <- update(examplejm3, assoc = data.frame("etavalue", "etaslope")), "'assoc' should be") 
  
  
})


test_that("basehaz argument works", {
  
  expect_output(update(examplejm1, basehaz = "weibull"))
  expect_output(update(examplejm1, basehaz = "bs"))
  expect_output(update(examplejm1, basehaz = "piecewise"))
  
  expect_output(update(examplejm1, basehaz = "bs", basehaz_ops = list(df = 5)))
  expect_output(update(examplejm1, basehaz = "bs", basehaz_ops = list(knots = c(1,3,5))))
  expect_output(update(examplejm1, basehaz = "piecewise", basehaz_ops = list(df = 5)))
  expect_output(update(examplejm1, basehaz = "piecewise", basehaz_ops = list(knots = c(1,3,5))))
  
  expect_output(expect_warning(update(examplejm1, basehaz = "weibull", basehaz_ops = list(df = 1)), "'df' will be ignored"))
  expect_output(expect_warning(update(examplejm1, basehaz = "weibull", basehaz_ops = list(knots = 1)), "'knots' will be ignored"))
  
  expect_output(update(examplejm1, basehaz = "piecewise", basehaz_ops = list(knots = c(1,3,5))))
  
  expect_error(update(examplejm1, basehaz = "bs", basehaz_ops = list(df = 1)), "must be atleast 3")
  expect_error(update(examplejm1, basehaz = "bs", basehaz_ops = list(knots = -1)), "'knots' must be non-negative")
  expect_error(update(examplejm1, basehaz = "piecewise", basehaz_ops = list(knots = -1)), "'knots' must be non-negative")
  expect_error(update(examplejm1, basehaz = "piecewise", basehaz_ops = list(knots = c(1,2,50))), "cannot be greater than the largest event time")
  
})


test_that("quadnodes argument works", {
  
  expect_output(update(examplejm1, quadnodes = 7))
  expect_output(update(examplejm1, quadnodes = 11))
  expect_output(update(examplejm1, quadnodes = 15))
  
  expect_error(update(examplejm1, quadnodes = 1), "'quadnodes' must be either 7, 11 or 15")
  expect_error(update(examplejm1, quadnodes = c(1,2)), "should be a numeric vector of length 1")
  expect_error(update(examplejm1, quadnodes = "wrong"), "should be a numeric vector of length 1")
  
})


test_that("weights argument works", {
  
  idvec0 <- pbcSurv[["id"]]
  idvec1 <- head(idvec0)            # missing IDs
  idvec2 <- rep(idvec0, each = 2)   # repeated IDs
  idvec3 <- c(idvec0, 9998, 9999)    # extra IDs not in model
  
  wts0 <- data.frame(id = idvec0, weights = rep_len(c(1,2), length(idvec0)))
  wts1 <- data.frame(id = idvec1, weights = rep_len(c(1,2), length(idvec1)))
  wts2 <- data.frame(id = idvec2, weights = rep_len(c(1,2), length(idvec2)))
  wts3 <- data.frame(id = idvec0, weights = rep_len(c(1,2), length(idvec0)),
                     junkcol = idvec0)
  wts4 <- data.frame(id = idvec0, weights = rep_len(c("word"), length(idvec0)))
  wts5 <- data.frame(id = idvec0, weights = rep_len(c(NA), length(idvec0)))
  wts6 <- data.frame(id = idvec0, weights = rep_len(c(-1, 1), length(idvec0)))
  wts7 <- data.frame(id = idvec3, weights = rep_len(c(1,2), length(idvec3)))
  
  expect_output(update(examplejm1, weights = wts0, iter = 5))
  expect_output(update(examplejm1, weights = wts7, iter = 5)) # ok to supply extra IDs in weights
  
  expect_error(update(examplejm1, weights = as.matrix(wts0)), "should be a data frame")
  expect_error(update(examplejm1, weights = wts1), "do not have weights supplied")
  expect_error(update(examplejm1, weights = wts2), "should only have one row")
  expect_error(update(examplejm1, weights = wts3), "should be a data frame with two columns")
  expect_error(update(examplejm1, weights = wts4), "weights supplied must be numeric")
  expect_error(update(examplejm1, weights = wts5), "weights supplied must be numeric")
  expect_error(update(examplejm1, weights = wts6), "Negative weights are not allowed")
  
})

test_that("init argument works", {
  expect_output(update(examplejm1, init = "model_based"))
  expect_output(update(examplejm1, init = "0"))
  expect_output(update(examplejm1, init = 0))
  expect_output(update(examplejm1, init = "random"))
})

test_that("prior_PD argument works", {
  expect_output(update(examplejm1, prior_PD = TRUE))
})

test_that("adapt_delta argument works", {
  expect_output(update(examplejm1, adapt_delta = NULL))
  expect_output(update(examplejm1, adapt_delta = 0.8))
})

test_that("max_treedepth argument works", {
  expect_output(update(examplejm1, max_treedepth = NULL))
  expect_output(update(examplejm1, max_treedepth = 5))
  expect_output(update(examplejm1, max_treedepth = 5L))
})

test_that("QR argument works", {
  expect_error(update(examplejm1, QR = TRUE), "not yet implemented")
})

#--------  Interactive tests

if (interactive()) {

  pbcLong$ybino_successes <- as.integer(cut(pbcLong$logBili, 6))
  pbcLong$ybino_failures  <- 6 - as.integer(cut(pbcLong$logBili, 6)) 
  pbcLong$ybino_trials    <- rep_len(6, length(pbcLong$ybino_successes))
  pbcLong$ybino_prop      <- pbcLong$ybino_successes / pbcLong$ybino_trials
  pbcLong$ybern <- as.integer(pbcLong$logBili >= mean(pbcLong$logBili))
  pbcLong$ypois <- as.integer(pbcLong$albumin)
  pbcLong$ygamm <- as.integer(pbcLong$platelet)
  pbcLong$xbern <- as.numeric(pbcLong$platelet / 100)
  pbcLong$xpois <- as.numeric(pbcLong$platelet / 100)
  pbcLong$xgamm <- as.numeric(pbcLong$logBili)
  
  compare_lmer <- function(fm) {
    y1 <- stan_lmer(fm, pbcLong, iter = 2000, chains = CHAINS, seed = SEED)
    s1 <- coxph(Surv(futimeYears, death) ~ sex + trt, data = pbcSurv)
    j1 <- stan_jm(fm, pbcLong, Surv(futimeYears, death) ~ sex + trt, pbcSurv,
                  time_var = "year", assoc = NULL, 
                  iter = 2000, chains = CHAINS, seed = SEED)        
    expect_equal(ranef(y1), ranef(j1)$Long1, tol = RANEF_tol)
    tols <- TOLSCALE_lmer * y1$ses
    for (i in 1:length(fixef(y1)))
      expect_equal(fixef(y1)[[i]], fixef(j1)$Long1[[i]], tol = tols[[i]])
    tols <- TOLSCALE_event * summary(s1)$coefficients[, "se(coef)"]
    for (i in 1:length(coef(s1)))
      expect_equal(coef(s1)[[i]], fixef(j1)$Event[[i+1]], tol = tols[[i]])        
  }  
    
  compare_glmer <- function(fm, fam) {
    SW(y1 <- stan_glmer(fm, pbcLong, fam, iter = 2000, chains = CHAINS, seed = SEED))
    SW(s1 <- coxph(Surv(futimeYears, death) ~ sex + trt, data = pbcSurv))
    SW(j1 <- stan_jm(fm, pbcLong, Surv(futimeYears, death) ~ sex + trt, pbcSurv,
                     time_var = "year", assoc = NULL, family = fam,
                     iter = 2000, chains = CHAINS, seed = SEED))
    
    # check parameter estimates
    expect_equal(ranef(y1), ranef(j1)$Long1, tol = RANEF_tol)
    tols <- TOLSCALE_glmer * y1$ses
    for (i in 1:length(fixef(y1)))
      expect_equal(fixef(y1)[[i]], fixef(j1)$Long1[[i]], tol = tols[[i]])
    tols <- TOLSCALE_event * summary(s1)$coefficients[, "se(coef)"]
    for (i in 1:length(coef(s1)))
      expect_equal(coef(s1)[[i]], fixef(j1)$Event[[i+1]], tol = tols[[i]])
    # check posterior_predict - in sample
    y1pp1 <- posterior_predict(y1, seed = SEED)
    j1pp1 <- posterior_predict(j1, seed = SEED, m = 1)
    expect_equal(y1pp1, j1pp1, tol = PP_tol)
    expect_equal(colMeans(y1pp1), colMeans(j1pp1), tol = PPMEANS_tol)
    expect_equal(colSds(y1pp1), colSds(j1pp1), tol = PPMEANS_tol)
    # check posterior_predict - new data
    y1pp2 <- posterior_predict(y1, newdata = pbcLong[1:10,], seed = SEED)
    j1pp2 <- posterior_predict(j1, newdata = pbcLong[1:10,], seed = SEED, m = 1)
    expect_equal(y1pp2, j1pp2, tol = PP_tol)
    expect_equal(colMeans(y1pp2), colMeans(j1pp2), tol = PPMEANS_tol)
    expect_equal(colSds(y1pp2), colSds(j1pp2), tol = PPMEANS_tol)
    # check posterior_predict - new data with new levels
    nd <- pbcLong[1:10,]
    nd$id <- paste0("new", nd$id)
    y1pp3 <- posterior_predict(y1, newdata = nd, seed = SEED)
    j1pp3 <- posterior_predict(j1, newdata = nd, seed = SEED, m = 1)   
    expect_equal(colMeans(y1pp3), colMeans(j1pp3), tol = PPMEANS_tol)
    expect_equal(colSds(y1pp3), colSds(j1pp3), tol = PPMEANS_tol)
    # check posterior_traj without interpolation
    nd <- pbcLong[pbcLong$id == 2, ]
    j1pp4 <- posterior_predict(j1, newdata = nd, seed = SEED, m = 1)
    j1pt4 <- posterior_traj(j1, newdata = nd, seed = SEED, return_matrix = TRUE,
                           interpolate = FALSE, extrapolate = FALSE)
    expect_identical(j1pp4, j1pt4)
    # check posterior_traj with interpolation
    nd <- pbcLong[pbcLong$id == 2, ]
    nd$year <- (0:8 / 8) * 10 # create exact interpolation times for posterior_predict
    j1pp5 <- posterior_predict(j1, newdata = nd, seed = SEED, m = 1)
    y1pt5 <- posterior_traj(j1, newdata = nd, seed = SEED, control = list(ipoints = 9))
    expect_equal(colMedians(j1pp5), y1pt5$yfit, tol = PPMEANS_tol) # not exactly equal because yfit is median of predicted mu, not predicted y
    expect_equal(colLower(j1pp5), y1pt5$pi_lb) # should be exactly equal, since pi_lb based on predicted y
    expect_equal(colUpper(j1pp5), y1pt5$pi_ub) # should be exactly equal, since pi_ub based on predicted y
  } 

  test_that("coefs same for stan_jm and stan_lmer/coxph", {
    compare_lmer(logBili ~ year + (1 | id))})  
  
  test_that("stan_jm/stan_glmer parameter estimates and predictions are the same, gaussian", {
    compare_glmer(logBili ~ year + (1 | id), gaussian)})  
  test_that("stan_jm/stan_glmer parameter estimates and predictions are the same, binomial as cbind(success, failure)", {
    compare_glmer(cbind(ybino_successes, ybino_failures) ~ year + xbern + (1 | id), binomial)})
  test_that("stan_jm/stan_glmer parameter estimates and predictions are the same, binomial as cbind(success, trials-success)", {
    compare_glmer(cbind(ybino_successes, ybino_trials - ybino_successes) ~ year + xbern + (1 | id), binomial)})
  test_that("stan_jm/stan_glmer parameter estimates and predictions are the same, bernoulli", {
    compare_glmer(ybern ~ year + xbern + (1 | id), binomial)})
  test_that("stan_jm/stan_glmer parameter estimates and predictions are the same, poisson", {
    compare_glmer(ypois ~ year + xpois + (1 | id), poisson)})
  test_that("stan_jm/stan_glmer parameter estimates and predictions are the same, negative binomial", {
    compare_glmer(ypois ~ year + xpois + (1 | id), neg_binomial_2)})
  #test_that("stan_jm/stan_glmer parameter estimates and predictions are the same, Gamma", {
  #compare_glmer(ygamm ~ year + xgamm + (1 | id), Gamma)})
  #test_that("stan_jm/stan_glmer parameter estimates and predictions are the same, inverse gaussian", {
  #compare_glmer(ygamm ~ year + xgamm + (1 | id), inverse.gaussian)})  

  test_that("coefs same for stan_jm and stan_glmer, binomial as prop as outcome and trials as weights", {
    jm_weights <- data.frame(id = pbcLong$id, weights = pbcLong$ybino_trials)
    y1 <- stan_glmer(ybino_prop ~ year + (1 | id), pbcLong, weights = ybino_trials, 
                     iter = 2000, chains = CHAINS, seed = SEED)
    s1 <- coxph(Surv(futimeYears, death) ~ sex + trt, data = pbcSurv)
    j1 <- stan_jm(ybino_prop ~ year + (1 | id), data = pbcLong, weights = jm_weights,
                  formulaEvent = Surv(futimeYears, death) ~ sex + trt,
                  dataLong = pbcLong, dataEvent = pbcSurv,
                  time_var = "year", assoc = NULL,
                  iter = 2000, chains = CHAINS, seed = SEED)        
    expect_equal(ranef(y1), ranef(j1)$Long1, tol = RANEF_tol)
    tols <- TOLSCALE_lmer * y1$ses
    for (i in 1:length(fixef(y1)))
      expect_equal(fixef(y1)[[i]], fixef(j1)$Long1[[i]], tol = tols[[i]])
    tols <- TOLSCALE_event * summary(s1)$coefficients[, "se(coef)"]
    for (i in 1:length(coef(s1)))
      expect_equal(coef(s1)[[i]], fixef(j1)$Event[[i+1]], tol = tols[[i]])       
  })  

  test_that("coefs same for stan_jm and JM, linear trajectory", {
    fmL1 <- logBili ~ year + (year | id)
    fmS1 <- Surv(futimeYears, death) ~ 1
    j1 <- SW(stan_jm(fmL1, pbcLong, fmS1, pbcSurv, time_var = "year", assoc = "etavalue", 
                  iter = 1000, chains = CHAINS, seed = SEED))
    l1 <- nlme::lme(logBili ~ year, random = ~ year | id, pbcLong)
    s1 <- coxph(fmS1, pbcSurv, x = TRUE, model = TRUE)  
    m1 <- JM::jointModel(l1, s1, timeVar = "year")
    # compare estimates and SEs for longitudinal submodel
    estL <- summary(m1)[["CoefTable-Long"]][, "Value"]
    sesL <- summary(m1)[["CoefTable-Long"]][, "Std.Err"]
    tols <- TOLSCALE_jms * sesL
    expect_equal(estL, fixef(j1)$Long1, check.attributes = FALSE, tol = tols)
    # compare estimates and SEs for event submodel
    estE <- summary(m1)[["CoefTable-Event"]][, "Value"]
    sesE <- summary(m1)[["CoefTable-Event"]][, "Std.Err"]
    tols <- TOLSCALE_jms * sesE
    expect_equal(estE[-length(estE)], fixef(j1)$Event, check.attributes = FALSE,
                 tol = tols[-length(estE)]) # the -length(estE) removes log(shape)
    expect_equal(sesE[-length(estE)], se(j1)$Event, check.attributes = FALSE,
                 tol = tols[-length(estE)]) # the -length(estE) removes log(shape)
    # compare estimate for log(shape)
    stan_logshape <- log(summary(j1, pars = "Event|weibull-shape")[,1])
    expect_equal(as.vector(estE[length(estE)]), stan_logshape, tol = tols[length(estE)])
  })
  
  test_that("coefs same for stan_jm and JM, raw polynomial trajectory", {
    fmL1 <- logBili ~ poly(year, degree = 2, raw = TRUE) + (poly(year, degree = 2, raw = TRUE) | id)
    fmS1 <- Surv(futimeYears, death) ~ 1
    SW(j1 <- stan_jm(fmL1, pbcLong, fmS1, pbcSurv, time_var = "year", assoc = "etavalue", 
                  iter = 1000, chains = CHAINS, seed = SEED))
    l1 <- nlme::lme(logBili ~ poly(year, degree = 2, raw = TRUE), 
                    random = ~ poly(year, degree = 2, raw = TRUE) | id, pbcLong)
    s1 <- coxph(fmS1, pbcSurv, x = TRUE, model = TRUE)  
    m1 <- JM::jointModel(l1, s1, timeVar = "year")
    # compare estimates and SEs for longitudinal submodel
    estL <- summary(m1)[["CoefTable-Long"]][, "Value"]
    sesL <- summary(m1)[["CoefTable-Long"]][, "Std.Err"]
    tols <- TOLSCALE_jms * sesL
    expect_equal(estL, fixef(j1)$Long1, check.attributes = FALSE, tol = tols)
    # compare estimates and SEs for event submodel
    estE <- summary(m1)[["CoefTable-Event"]][, "Value"]
    sesE <- summary(m1)[["CoefTable-Event"]][, "Std.Err"]
    tols <- TOLSCALE_jms * sesE
    expect_equal(estE[-length(estE)], fixef(j1)$Event, check.attributes = FALSE,
                 tol = tols[-length(estE)]) # the -length(estE) removes log(shape)
    expect_equal(sesE[-length(estE)], se(j1)$Event, check.attributes = FALSE,
                 tol = tols[-length(estE)]) # the -length(estE) removes log(shape)
    # compare estimate for log(shape)
    stan_logshape <- log(summary(j1, pars = "Event|weibull-shape")[,1])
    expect_equal(as.vector(estE[length(estE)]), stan_logshape, tol = tols[length(estE)])
  })  
 
  test_that("coefs same for stan_jm and JM, orthogonalised polynomial trajectory", {
    fmL1 <- logBili ~ poly(year, degree = 2) + (poly(year, degree = 2) | id)
    fmS1 <- Surv(futimeYears, death) ~ 1
    j1 <- SW(stan_jm(fmL1, pbcLong, fmS1, pbcSurv, time_var = "year", assoc = "etavalue", 
                  iter = 1000, chains = CHAINS, seed = SEED))
    l1 <- nlme::lme(logBili ~ poly(year, degree = 2), 
                    random = ~ poly(year, degree = 2) | id, pbcLong)
    s1 <- coxph(fmS1, pbcSurv, x = TRUE, model = TRUE)  
    m1 <- JM::jointModel(l1, s1, timeVar = "year")
    # compare estimates and SEs for longitudinal submodel
    estL <- summary(m1)[["CoefTable-Long"]][, "Value"]
    sesL <- summary(m1)[["CoefTable-Long"]][, "Std.Err"]
    tols <- TOLSCALE_jms * sesL
    expect_equal(estL, fixef(j1)$Long1, check.attributes = FALSE, tol = tols)
    # compare estimates and SEs for event submodel
    estE <- summary(m1)[["CoefTable-Event"]][, "Value"]
    sesE <- summary(m1)[["CoefTable-Event"]][, "Std.Err"]
    tols <- TOLSCALE_jms * sesE
    expect_equal(estE[-length(estE)], fixef(j1)$Event, check.attributes = FALSE,
                 tol = tols[-length(estE)]) # the -length(estE) removes log(shape)
    expect_equal(sesE[-length(estE)], se(j1)$Event, check.attributes = FALSE,
                 tol = tols[-length(estE)]) # the -length(estE) removes log(shape)
    # compare estimate for log(shape)
    stan_logshape <- log(summary(j1, pars = "Event|weibull-shape")[,1])
    expect_equal(as.vector(estE[length(estE)]), stan_logshape, tol = tols[length(estE)])
  })    
  
  
  
   
 compare_jm <- function(fm, fmFixed, fmRandom) {
   fmS1 <- Surv(futimeYears, death) ~ 1
   SW(j1 <- stan_jm(fm, pbcLong, fmS1, pbcSurv, time_var = "year", assoc = "etavalue", 
                 iter = 1000, chains = CHAINS, seed = SEED))
   l1 <- nlme::lme(fmFixed, random = fmRandom, pbcLong)
   s1 <- coxph(fmS1, pbcSurv, x = TRUE, model = TRUE)  
   m1 <- JM::jointModel(l1, s1, timeVar = "year")
   # compare estimates and SEs for longitudinal submodel
   estL <- summary(m1)[["CoefTable-Long"]][, "Value"]
   sesL <- summary(m1)[["CoefTable-Long"]][, "Std.Err"]
   tols <- TOLSCALE_jms * sesL
   expect_equal(estL, fixef(j1)$Long1, check.attributes = FALSE, tol = tols)
   # compare estimates and SEs for event submodel
   estE <- summary(m1)[["CoefTable-Event"]][, "Value"]
   sesE <- summary(m1)[["CoefTable-Event"]][, "Std.Err"]
   tols <- TOLSCALE_jms * sesE
   expect_equal(estE[-length(estE)], fixef(j1)$Event, check.attributes = FALSE,
                tol = tols[-length(estE)]) # the -length(estE) removes log(shape)
   expect_equal(sesE[-length(estE)], se(j1)$Event, check.attributes = FALSE,
                tol = tols[-length(estE)]) # the -length(estE) removes log(shape)
   # compare estimate for log(shape)
   stan_logshape <- log(summary(j1, pars = "Event|weibull-shape")[,1])
   expect_equal(as.vector(estE[length(estE)]), stan_logshape, tol = tols[length(estE)])   
 }

 test_that("coefs same for stan_jm and JM, linear trajectory", {
   compare_jm(logBili ~ year + (1 | id), logBili ~ year, ~ 1 | id)})
 test_that("coefs same for stan_jm and JM, raw polynomial trajectory", {
   compare_jm(logBili ~ poly(year, degree = 2, raw = TRUE) + (poly(year, degree = 2, raw = TRUE) | id), 
              logBili ~ poly(year, degree = 2, raw = TRUE), ~ poly(year, degree = 2, raw = TRUE) | id)})
 test_that("coefs same for stan_jm and JM, orthogonalised polynomial trajectory", {
   compare_jm(logBili ~ year + (1 | id), logBili ~ year, ~ 1 | id)})
 
  
   
}

