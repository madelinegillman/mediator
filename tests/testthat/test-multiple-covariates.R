context("multiple covariates")

# import data for use in tests - roughly based off data from SAS macro --------
load("test_dat.RData")

# set uniform values to test across all types ---------------------------------
data = test_dat
a = 1
a_star = 0
m = 0

test_that("continuous outcome and mediator match SAS macro",{

  # setting options (specificed by user) --------------------------------------
  out.model = glm(y ~ x + c + m + cens + x*c, data = test_dat)
  med.model = glm(c ~ x + m + cens, data = test_dat)
  treat = "x"
  mediator = "c"
  out.reg = "linear"
  med.reg = "linear"

  # calculating values to use later on ----------------------------------------
  out_vars <- if (out.reg=="coxph") names(attr(out.model$terms,"dataClasses"))[-1] else
    names(attr(out.model$terms,"dataClasses"))

  var_set <- unique(c(out_vars,
                      names(attr(med.model$terms,"dataClasses"))))
  data <- data %>% dplyr::select(all_of(var_set))

  betas <- stats::coef(med.model) # coefficients from mediation model
  beta_info <- cov_pred(treat, mediator, med.model, data)
  betasum <- sum(beta_info$betasum, na.rm=TRUE)
  betameans <- beta_info$betamean

  # get covariate names
  cnames <- names(betameans)

  # Covariance matrix for standar errors
  sigmaV <- stats::sigma(med.model)^2
  Sigma <- comb_sigma(med.model, out.model, treat, mediator,
                      out.reg, cnames, med.reg)

  # setting coefficients for no interaction = 0 -------------------------------
  if(is.na(out.model$coefficients[paste0(treat, ":", mediator)])){
    out.model$coefficients[paste0(treat, ":", mediator)] <- 0
  } else {
    out.model$coefficients[paste0(treat, ":", mediator)] <-
      out.model$coefficients[paste0(treat, ":", mediator)]
  }

  # pulling coefficients from models
  theta1 <- out.model$coefficients[treat]
  theta2 <- out.model$coefficients[mediator]
  theta3 <- out.model$coefficients[paste0(treat, ":", mediator)]

  beta0 <- med.model$coefficients["(Intercept)"]
  beta1 <- med.model$coefficients[treat]

  arg_list <- list(theta1 = theta1, theta2 = theta2, theta3 = theta3,
                   beta0 = beta0, beta1 = beta1,
                   betasum = betasum, betameans = betameans,
                   a = a, a_star = a_star, m = m, out.reg = out.reg,
                   med.reg = med.reg,
                   sigmaV = sigmaV)

  # calculate effect estimates ------------------------------------------------

  ## controlled direct effect
  CDE <- do.call(controlled_direct_effect, arg_list)
  expect_equal(round(CDE, 6), 0.091774)

  ## natural direct effect
  NDE <- do.call(natural_direct_effect, arg_list)
  expect_equal(round(as.numeric(NDE), 6), 0.084069)

  ## natural indirect effect
  NIE <- do.call(natural_indirect_effect, arg_list)
  expect_equal(round(as.numeric(NIE), 6), -0.009772)

  ## total effect
  TE <- total_effect(NDE, NIE, out.reg)
  expect_equal(round(as.numeric(TE), 6), 0.074297)

  ## proportion mediated
  PM <- prop_mediated(NDE, NIE, out.reg, TE)
  expect_equal(round(as.numeric(PM), 5), -0.13153)

})

test_that("binary outcome and continuous mediator match SAS macro",{

  # setting options (specificed by user) --------------------------------------
  out.model = glm(cens ~ x + y + m + c + x*y, data = test_dat,
                  family="binomial")
  med.model = glm(y ~ x + m + c, data = test_dat)
  treat = "x"
  mediator = "y"
  out.reg = "logistic"
  med.reg = "linear"

  # calculating values to use later on ----------------------------------------
  out_vars <- if (out.reg=="coxph") names(attr(out.model$terms,"dataClasses"))[-1] else
    names(attr(out.model$terms,"dataClasses"))

  var_set <- unique(c(out_vars,
                      names(attr(med.model$terms,"dataClasses"))))
  data <- data %>% dplyr::select(all_of(var_set))

  betas <- stats::coef(med.model) # coefficients from mediation model
  beta_info <- cov_pred(treat, mediator, med.model, data)
  betasum <- sum(beta_info$betasum, na.rm=TRUE)
  betameans <- beta_info$betamean

  # get covariate names
  cnames <- names(betameans)

  # Covariance matrix for standar errors
  sigmaV <- stats::sigma(med.model)^2
  Sigma <- comb_sigma(med.model, out.model, treat, mediator,
                      out.reg, cnames, med.reg)

  # setting coefficients for no interaction = 0 -------------------------------
  if(is.na(out.model$coefficients[paste0(treat, ":", mediator)])){
    out.model$coefficients[paste0(treat, ":", mediator)] <- 0
  } else {
    out.model$coefficients[paste0(treat, ":", mediator)] <-
      out.model$coefficients[paste0(treat, ":", mediator)]
  }

  # pulling coefficients from models
  theta1 <- out.model$coefficients[treat]
  theta2 <- out.model$coefficients[mediator]
  theta3 <- out.model$coefficients[paste0(treat, ":", mediator)]

  beta0 <- med.model$coefficients["(Intercept)"]
  beta1 <- med.model$coefficients[treat]

  arg_list <- list(theta1 = theta1, theta2 = theta2, theta3 = theta3,
                   beta0 = beta0, beta1 = beta1,
                   betasum = betasum, betameans = betameans,
                   a = a, a_star = a_star, m = m, out.reg = out.reg,
                   med.reg = med.reg,
                   sigmaV = sigmaV)
  # calculate effect estimates ------------------------------------------------

  ## controlled direct effect
  CDE <-  do.call(controlled_direct_effect, arg_list)
  expect_equal(round(CDE, 5), 0.42995)

  ## natural direct effect
  NDE <- do.call(natural_direct_effect, arg_list)
  expect_equal(round(as.numeric(NDE), 5), 0.72323)

  ## natural indirect effect
  NIE <- do.call(natural_indirect_effect, arg_list)
  expect_equal(round(as.numeric(NIE), 5), 1.05791, tolerance = 0.00001)

  ## total effect
  TE <- total_effect(NDE, NIE, out.reg)
  expect_equal(round(as.numeric(TE), 5), 0.76512)

  ## proportion mediated
  PM <- prop_mediated(NDE, NIE, out.reg, TE)
  expect_equal(round(as.numeric(PM), 5), -0.17833)

})

test_that("continuous outcome and binary mediator match SAS macro",{

  # setting options (specificed by user) --------------------------------------
  out.model = glm(y ~ x + m + c + cens + x*m, data = test_dat)
  med.model = glm(m ~ x + c + cens, data = test_dat, family = "binomial")
  treat = "x"
  mediator = "m"
  out.reg = "linear"
  med.reg = "logistic"

  # calculating values to use later on ----------------------------------------
  out_vars <- if (out.reg=="coxph") names(attr(out.model$terms,"dataClasses"))[-1] else
    names(attr(out.model$terms,"dataClasses"))

  var_set <- unique(c(out_vars,
                      names(attr(med.model$terms,"dataClasses"))))
  data <- data %>% dplyr::select(all_of(var_set))

  betas <- stats::coef(med.model) # coefficients from mediation model
  # beta_info <- cov_pred(cmeans, cmodes, treat, mediator, med.model, data)
  beta_info <- cov_pred(treat, mediator, med.model, data)
  betasum <- sum(beta_info$betasum, na.rm=TRUE)
  betameans <- beta_info$betamean

  # get covariate names
  cnames <- names(betameans)

  # Covariance matrix for standar errors
  sigmaV <- stats::sigma(med.model)^2
  Sigma <- comb_sigma(med.model, out.model, treat, mediator,
                      out.reg, cnames, med.reg)

  # setting coefficients for no interaction = 0 -------------------------------
  if(is.na(out.model$coefficients[paste0(treat, ":", mediator)])){
    out.model$coefficients[paste0(treat, ":", mediator)] <- 0
  } else {
    out.model$coefficients[paste0(treat, ":", mediator)] <-
      out.model$coefficients[paste0(treat, ":", mediator)]
  }

  # pulling coefficients from models
  theta1 <- out.model$coefficients[treat]
  theta2 <- out.model$coefficients[mediator]
  theta3 <- out.model$coefficients[paste0(treat, ":", mediator)]

  beta0 <- med.model$coefficients["(Intercept)"]
  beta1 <- med.model$coefficients[treat]

  arg_list <- list(theta1 = theta1, theta2 = theta2, theta3 = theta3,
                   beta0 = beta0, beta1 = beta1,
                   betasum = betasum, betameans = betameans,
                   a = a, a_star = a_star, m = m, out.reg = out.reg,
                   med.reg = med.reg,
                   sigmaV = sigmaV)

  # calculate effect estimates ------------------------------------------------

  ## controlled direct effect
  CDE <- do.call(controlled_direct_effect, arg_list)
  expect_equal(round(CDE, 6), 0.26042, tolerance = 0.00001)

  ## natural direct effect
  NDE <- do.call(natural_direct_effect, arg_list)
  expect_equal(round(as.numeric(NDE), 6), 0.09575, tolerance = 0.00001)

  ## natural indirect effect
  NIE <- do.call(natural_indirect_effect, arg_list)
  expect_equal(round(as.numeric(NIE), 6), -0.00369, tolerance = 0.00001)

  ## total effect
  TE <- total_effect(NDE, NIE, out.reg)
  expect_equal(round(as.numeric(TE), 6), 0.09206, tolerance = 0.00001)

  ## proportion mediated
  PM <- prop_mediated(NDE, NIE, out.reg, TE)
  expect_equal(round(as.numeric(PM), 5), -0.040117, tolerance = 0.00001)
})

test_that("binary outcome and mediator match SAS macro",{

  # setting options (specificed by user) --------------------------------------
  out.model = glm(cens ~ x + m + c + y + x*m, data = test_dat, family = "binomial")
  med.model = glm(m ~ x + c + y,data = test_dat, family = "binomial")
  treat = "x"
  mediator = "m"
  out.reg = "logistic"
  med.reg = "logistic"

  # calculating values to use later on ----------------------------------------
  out_vars <- if (out.reg=="coxph") names(attr(out.model$terms,"dataClasses"))[-1] else
    names(attr(out.model$terms,"dataClasses"))

  var_set <- unique(c(out_vars,
                      names(attr(med.model$terms,"dataClasses"))))
  data <- data %>% dplyr::select(all_of(var_set))

  betas <- stats::coef(med.model) # coefficients from mediation model
  # beta_info <- cov_pred(cmeans, cmodes, treat, mediator, med.model, data)
  beta_info <- cov_pred(treat, mediator, med.model, data)
  betasum <- sum(beta_info$betasum, na.rm=TRUE)
  betameans <- beta_info$betamean

  # get covariate names
  cnames <- names(betameans)

  # Covariance matrix for standar errors
  sigmaV <- stats::sigma(med.model)^2
  Sigma <- comb_sigma(med.model, out.model, treat, mediator,
                      out.reg, cnames, med.reg)

  # setting coefficients for no interaction = 0 -------------------------------
  if(is.na(out.model$coefficients[paste0(treat, ":", mediator)])){
    out.model$coefficients[paste0(treat, ":", mediator)] <- 0
  } else {
    out.model$coefficients[paste0(treat, ":", mediator)] <-
      out.model$coefficients[paste0(treat, ":", mediator)]
  }

  # pulling coefficients from models
  theta1 <- out.model$coefficients[treat]
  theta2 <- out.model$coefficients[mediator]
  theta3 <- out.model$coefficients[paste0(treat, ":", mediator)]

  beta0 <- med.model$coefficients["(Intercept)"]
  beta1 <- med.model$coefficients[treat]

  arg_list <- list(theta1 = theta1, theta2 = theta2, theta3 = theta3,
                   beta0 = beta0, beta1 = beta1,
                   betasum = betasum, betameans = betameans,
                   a = a, a_star = a_star, m = m, out.reg = out.reg,
                   med.reg = med.reg,
                   sigmaV = sigmaV)

  # calculate effect estimates ------------------------------------------------

  ## controlled direct effect
  CDE <- do.call(controlled_direct_effect, arg_list)
  expect_equal(round(CDE, 5), 0.54184, tolerance = 0.00001)

  ## natural direct effect
  NDE <- do.call(natural_direct_effect, arg_list)
  expect_equal(round(as.numeric(NDE), 5), 0.61768, tolerance = 0.00001)

  ## natural indirect effect
  NIE <- do.call(natural_indirect_effect, arg_list)
  expect_equal(round(as.numeric(NIE), 5), 1.00201, tolerance = 0.00001)

  ## total effect
  TE <- total_effect(NDE, NIE, out.reg)
  expect_equal(round(as.numeric(TE), 5), 0.61892, tolerance = 0.00001)

  ## proportion mediated
  PM <- prop_mediated(NDE, NIE, out.reg, TE)
  expect_equal(round(as.numeric(PM), 6), -0.003255, tolerance = 0.00001)
})
