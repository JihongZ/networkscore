library(lavaan)
## Centrality: Strength
### x: precision matrix with diagonal value is 0
centralityTableNew <- function(x) {
  # Strength is sum of absolute values of weights
  raw_strength <- rowSums(abs(x))
  std_strength <- scale(raw_strength)
  centralityTable(x) |>
    filter(measure == "Strength") |>
    select(node:value) |>
    mutate(selcalc.std.value = std_strength,
           selcalc.raw.value = raw_strength)
}


# Make data into factor score and Lambda---------------------------------------------
PredFactorScore <- function(dat){
  colnames(dat) <- paste0("X", 1:ncol(dat))
  mod <- paste0("F1=~", paste0(colnames(dat), collapse = "+"))
  fit <- cfa(model = mod, data = dat, std.lv = FALSE)
  parameters = parameterEstimates(fit, standardized = TRUE)
  Lambda = parameters[parameters$op == "=~", "std.all"]
  Psi = parameters[(parameters$op == "~~" & parameters$lhs != "F1"), "std.all"]
  Psi = diag(Psi)
  # browser()
  fs <- as.numeric(lavaan::lavPredict(fit, method = "regression"))

  list(
    fit = fit, # lavaan object
    factorScore = fs, # estimates of factor scores
    Lambda = Lambda, # estimates of factor loadings
    Cor = t(Lambda) %*% solve(Lambda %*% t(Lambda) + Psi) %*% Lambda # unscaled correlation
  )
}


# Make correlation vector into Network/Factor Score ------------------------------
cor_vecTomat <- function(cor_vec) {
  cor_mat <- diag(1, length(cor_vec))
  cor_mat[upper.tri(cor_mat)] <- cor_vec
  cor_mat[lower.tri(cor_mat)] <- t(cor_mat)[lower.tri(t(cor_mat))]
  return(cor_mat)
}

## input original data, fit a unidimensional CFA and a saturated network model
datToScore <- function(dat){
  Y = dat
  ## Calculate network scores
  Results <- networkscore::networkscore(Y)
  ## Fit a GGM model
  FS <- scale(PredFactorScore(Y)$factorScore) # calculate Factor scores
  cbind(
    as.data.frame(Results),
    FS = FS
  )
}

# (1) the sample size (n = 500, n = 1000, n = 2000);
# (2) the test length (J = 6, 12 or 24) was selected to reflect a range of potential applications spanning small to large;
# (3) the measurement quality, as indicated by the loadings (all $\lambda$s = .7 or all $\lambda s$ = .4);
# (4) the residual covariances between the first indicator and the second indicator ($\boldsymbol\Psi_{12}$ = 0 or $\boldsymbol\Psi_{12}$ = .2)
data.generation <- function(N, J, Lambda, Psi12){
  if (length(Lambda) == 1) {
    Lambda_mat <- matrix(rep(Lambda, J), ncol = 1)
  }else{
    Lambda_mat <- matrix(Lambda, ncol = 1)
  }
  Mu_mat <- matrix(0, nrow = N, ncol = J) # item intercept
  Residual_Sigma <- diag(1-Lambda^2, J) # residual variance covariance matrix

  Residual_Sigma[1,2] = Residual_Sigma[2,1] = Psi12
  Residual <- mvtnorm::rmvnorm(N, mean = rep(0, J), sigma = Residual_Sigma)
  FactorScore <- rnorm(N, mean = 0, sd = 1)

  myData <- FactorScore %*% t(Lambda_mat) + Mu_mat + Residual
  list(
    data= myData,
    true_factor_score = FactorScore
  )
}
