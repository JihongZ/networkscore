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

networkscore <- function(dat) {
  out <- list()

  # Model 0: saturated model
  mod0 <- psychonetrics::ggm(dat)
  fit0 <- mod0 |> psychonetrics::runmodel()
  Delta = psychonetrics::getmatrix(fit0, "delta")
  Omega = psychonetrics::getmatrix(fit0, "omega")
  Scaled_Y = solve(Delta) %*% (diag(1, ncol(dat)) - Omega) %*% solve(Delta) %*% t(dat)

  ## Strength based Network Score
  Node_Strength_mod0  = rowSums(abs(Omega))
  out$NS_S <- as.numeric(Node_Strength_mod0 %*% Scaled_Y) # NS based on raw strength
  out$NS_RS <- as.numeric(sqrt(Node_Strength_mod0) %*% Scaled_Y) # NS based on square root of strength

  ## Hybrid based Network Score
  Hybrid_mod0 = NetworkToolbox::hybrid(Omega)
  out$NS_H <- as.numeric(Hybrid_mod0 %*% Scaled_Y)
  out$NS_RH <- as.numeric(sqrt(Hybrid_mod0) %*% Scaled_Y)

  # Model 1: regularized model
  fit1 <- mod0 |> psychonetrics::prune() |> psychonetrics::runmodel()
  Delta1 = psychonetrics::getmatrix(fit1, "delta")
  Omega1 = psychonetrics::getmatrix(fit1, "omega")
  Scaled_Y1 = solve(Delta1) %*% (diag(1, ncol(dat)) - Omega1) %*% solve(Delta1) %*% t(dat)

  ## Strength-based Network Score for regularized model
  Node_Strength_mod1  = rowSums(abs(Omega1))
  out$RegNS_S <- as.numeric(Node_Strength_mod1 %*% Scaled_Y1) # Regularized NS based on strength
  out$RegNS_RS <- as.numeric(sqrt(Node_Strength_mod1) %*% Scaled_Y1) # Regularized NS based on square roots of strength

  ## Hybrid based Network Score for regularized model
  Omega1[Omega1 == 0] <- .001
  diag(Omega1) = 0
  Hybrid_mod1 = NetworkToolbox::hybrid(Omega1)
  out$RegNS_H <- as.numeric(Hybrid_mod1 %*% Scaled_Y1)
  out$RegNS_RH <- as.numeric(sqrt(Hybrid_mod1) %*% Scaled_Y1)

  ## Linear Combination
  out$NS_Origin <- as.numeric(Hybrid_mod0 %*% t(dat))
  out$NS_Origin_S <- as.numeric(Node_Strength_mod1 %*% t(dat))

  out
}
