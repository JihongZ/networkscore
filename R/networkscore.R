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

centrality_rescale <- function(centrality){
  (centrality - min(centrality)) / (sd(centrality))
  # (rank(centrality) - 1)/(length(centrality)-1)
}
centrality_rank <- function(centrality){
  #(centrality - min(centrality)) / (max(centrality) - min(centrality))
  (rank(centrality) - 1)/(length(centrality)-1)
}

(centrality <- runif(10, .8, .9))
(centrality2 <- centrality + .1)
rank <- centrality_rank(centrality)
unique(rank)[order(unique(rank))]

rank2 <- centrality_rank(centrality2)
plot(rank, rank2)

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
  Delta = psychonetrics::getmatrix(fit0, "delta")
  Kappa = psychonetrics::getmatrix(fit0, "kappa")
  Scaled_Y = dat%*%Kappa

  ## Model 0 Score: Strength based Network Score
  S_mod0  = centrality_rescale(rowSums(abs(Omega)))
  EI_mod0  = centrality_rescale(rowSums(Omega))
  BC_mod0  = centrality_rescale(NetworkToolbox::betweenness(Omega)) # betweenness centrality (BC);
  LC_mod0  = centrality_rescale(NetworkToolbox::closeness(Omega)) # closeness centrality;
  S_mod0_rank  = centrality_rank(rowSums(abs(Omega)))
  EI_mod0_rank  = centrality_rank(rowSums(Omega))
  BC_mod0_rank  = centrality_rank(NetworkToolbox::betweenness(Omega)) # betweenness centrality (BC);
  LC_mod0_rank  = centrality_rank(NetworkToolbox::closeness(Omega)) # closeness centrality;
  HB_mod0 = hb_centrality(Omega)

  out$NS_S <- as.numeric(S_mod0 %*% t(dat)) # NS based on raw strength
  out$NS_EI <- as.numeric(EI_mod0 %*% t(dat)) # NS based on raw strength
  out$NS_BC <- as.numeric(BC_mod0 %*% t(dat)) # BC based on raw strength
  out$NS_LC <- as.numeric(LC_mod0 %*% t(dat)) # LC based on raw strength
  out$NS_S_rank <- as.numeric(S_mod0_rank %*% t(dat)) # NS based on raw strength
  out$NS_EI_rank <- as.numeric(EI_mod0_rank %*% t(dat)) # NS based on raw strength
  out$NS_BC_rank <- as.numeric(BC_mod0_rank %*% t(dat)) # BC based on raw strength
  out$NS_LC_rank <- as.numeric(LC_mod0_rank %*% t(dat)) # LC based on raw strength
  out$NS_H <- as.numeric(HB_mod0 %*% t(dat)) ## Hybrid based Network Score

  # Model 1: regularized model
  fit1 <- mod0 |> psychonetrics::prune() |> psychonetrics::runmodel()
  Omega1 = psychonetrics::getmatrix(fit1, "omega")
  Delta1 = psychonetrics::getmatrix(fit1, "delta")
  Kappa1 = psychonetrics::getmatrix(fit1, "kappa")
  Scaled_Y1 = dat%*%Kappa1

  ## Model 1 Score: Strength-based Network Score for regularized model
  # S_mod1  = rowSums(abs(Omega1))
  S_mod1  = centrality_rank(rowSums(abs(Omega1)))
  EI_mod1  = centrality_rank(rowSums(Omega1))
  BC_mod1  = centrality_rank(NetworkToolbox::betweenness(Omega1)) # betweenness centrality (BC);
  LC_mod1  = centrality_rank(NetworkToolbox::closeness(Omega1)) # closeness centrality;
  HB_mod1 = hb_centrality(Omega1)

  out$RegNS_S <- as.numeric(S_mod1 %*% t(dat)) # Regularized NS based on strength
  out$RegNS_EI <- as.numeric(EI_mod1 %*% t(dat)) # Regularized NS based on strength
  out$RegNS_BC <- as.numeric(BC_mod1 %*% t(dat)) # BC based on raw strength
  out$RegNS_LC <- as.numeric(LC_mod1 %*% t(dat)) # LC based on raw strength
  out$RegNS_H <- as.numeric(HB_mod1 %*% t(dat))

  ## Linear Combination Christensen (2018)
  # out$NS_S_0 <- as.numeric(S_mod0 %*% t(dat))
  # out$NS_EI_0 <- as.numeric(EI_mod0 %*% t(dat))
  # out$NS_BC_0 <- as.numeric(BC_mod0 %*% t(dat)) # NS based on raw strength
  # out$NS_LC_0 <- as.numeric(LC_mod0 %*% t(dat)) # NS based on raw strength
  # out$NS_H_0 <- as.numeric(HB_mod0 %*% t(dat))

  out <- lapply(out, \(x) as.numeric(scale(x))) # standarized
  out
}
