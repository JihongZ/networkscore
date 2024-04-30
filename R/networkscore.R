# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'
#' @export

transform_centrality <- function(centrality){
  centrality = centrality / (length(centrality) - 1)
  (centrality - min(centrality)) / (sd(centrality))
}

hb_centrality <- function(omega){
  omega[omega == 0] <- .0001
  diag(omega) = 0
  NetworkToolbox::hybrid(omega)
}

networkscore <- function(dat) {
  out <- list()
  # Model 0: saturated model
  mod0 <- psychonetrics::ggm(dat, estimator = "ML", standardize = 'z')
  fit0 <- mod0 |> psychonetrics::runmodel()
  Omega = psychonetrics::getmatrix(fit0, "omega")
  # Delta = psychonetrics::getmatrix(fit0, "delta")
  # Kappa = psychonetrics::getmatrix(fit0, "kappa")
  # Scaled_Y = Kappa %*% t(dat)

  ## Model 0 Score: Strength based Network Score
  # S_mod0  = rowSums(abs(Omega))
  S_mod0  = transform_centrality(rowSums(abs(Omega)))
  EI_mod0  = transform_centrality(rowSums(Omega))
  HB_mod0 = hb_centrality(Omega)

  out$NS_S <- as.numeric(S_mod0 %*% t(dat)) # NS based on raw strength
  out$NS_RS <- as.numeric(sqrt(S_mod0) %*% t(dat)) # NS based on square root of strength
  out$NS_EI <- as.numeric(EI_mod0 %*% t(dat)) # NS based on raw strength
  out$NS_REI <- as.numeric(sqrt(EI_mod0) %*% t(dat)) # NS based on square root of strength
  out$NS_H <- as.numeric(HB_mod0 %*% t(dat)) ## Hybrid based Network Score

  # Model 1: regularized model
  fit1 <- mod0 |> psychonetrics::prune() |> psychonetrics::runmodel()
  Omega1 = psychonetrics::getmatrix(fit1, "omega")
  # Delta1 = psychonetrics::getmatrix(fit1, "delta")
  # Kappa1 = psychonetrics::getmatrix(fit1, "kappa")
  # Scaled_Y1 = Kappa1 %*% t(dat)

  ## Model 1 Score: Strength-based Network Score for regularized model
  # S_mod1  = rowSums(abs(Omega1))
  S_mod1  = transform_centrality(rowSums(abs(Omega1)))
  EI_mod1  = transform_centrality(rowSums(Omega1))
  HB_mod1 = hb_centrality(Omega1)

  out$RegNS_S <- as.numeric(S_mod1 %*% t(dat)) # Regularized NS based on strength
  out$RegNS_RS <- as.numeric(sqrt(S_mod1) %*% t(dat)) # Regularized NS based on square roots of strength
  out$RegNS_EI <- as.numeric(EI_mod1 %*% t(dat)) # Regularized NS based on strength
  out$RegNS_REI <- as.numeric(sqrt(EI_mod1) %*% t(dat)) # Regularized NS based on square roots of strength
  out$RegNS_H <- as.numeric(HB_mod1 %*% t(dat))

  ## Hybrid based Network Score for regularized model
  # out$RegNS_RH <- as.numeric(sqrt(Hybrid_mod1) %*% Scaled_Y1)

  ## Linear Combination Christensen (2018)
  # out$NS_H_0 <- as.numeric(HB_mod0 %*% t(dat))
  # out$NS_S_0 <- as.numeric(S_mod0 %*% t(dat))
  # out$NS_RS_0 <- as.numeric(sqrt(S_mod0) %*% t(dat))
  # out$NS_EI_0 <- as.numeric(EI_mod0 %*% t(dat))
  # out$NS_REI_0 <- as.numeric(sqrt(EI_mod0) %*% t(dat))

  out <- lapply(out, \(x) as.numeric(scale(x))) # standarized
  out
}
