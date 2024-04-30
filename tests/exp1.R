## Centrality measures is a secondary analysis transformed from precision matrix.
library(greekLetters)
library(here)
library(lavaan)
library(mvtnorm)
library(networkscore) # for network score , install the package via `devtools::install_github("JihongZ/networkscore")`
library(parSim) # devtools::install_github("sachaepskamp/parSim")
library(psych)
library(psychonetrics)
library(qgraph)
library(tidyverse)
source("R/Helper.R")

J = 12
OneRep <- data.generation(N = 300, J = J, Lambda = runif(J, .3, .7), Psi12 = 0.2) # generate data and true factor scores
sapply(as.data.frame(networkscore(OneRep$data)), \(x) cor(x, OneRep$true_factor_score))
