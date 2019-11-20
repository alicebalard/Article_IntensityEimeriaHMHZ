source("DataPreparation.R")

# D. Compare old transect / new transect

# To compare both dataset hybrid effect (alpha) we run 2 models with our data:
#
#   * model A with fixed alpha = 1.39 (Bairds et al. 2012 value for pinworms)
#
#   * model B with variable alpha (model B being fit_BAL_abundance$H1)
#
#   then we compare via G-test.

FitAdvancedNoAlphaNegbin <- function(data, response, hybridIndex, paramBounds, config){
  print("Fitting model advanced without alpha")
  data$response <- data[[response]]
  HI <- data[[hybridIndex]]
  start <-  list(L1 = paramBounds[["L1start"]],
                 L2 = paramBounds[["L2start"]],
                 A1 = paramBounds[["A1start"]],
                 A2 = paramBounds[["A2start"]],
                 Z = paramBounds[["Zstart"]])
  fit <- bbmle::mle2(
    response ~ dnbinom(mu = MeanLoad(L1, L2, 0, HI),
                       size = SizeNegBin(A1, A2, Z, HI)),
    data = data,
    start = start,
    lower = c(L1 = paramBounds[["L1LB"]],
              L2 = paramBounds[["L2LB"]],
              A1 = paramBounds[["A1LB"]],
              A2 = paramBounds[["A2LB"]],
              Z = paramBounds[["ZLB"]]),
    upper = c(L1 = paramBounds[["L1UB"]],
              L2 = paramBounds[["L2UB"]],
              A1 = paramBounds[["A1UB"]],
              A2 = paramBounds[["A2UB"]],
              Z = paramBounds[["ZUB"]]),
    optimizer = config$optimizer,
    method = config$method,
    control = config$control)
  printConvergence(fit)
  return(fit)
}

FitAdvancedAlphaNegbin <- function(data, response, hybridIndex, paramBounds, config){
  print("Fitting model advanced with alpha")
  data$response <- data[[response]]
  HI <- data[[hybridIndex]]
  start <-  list(L1 = paramBounds[["L1start"]],
                 L2 = paramBounds[["L2start"]],
                 A1 = paramBounds[["A1start"]],
                 A2 = paramBounds[["A2start"]],
                 alpha = paramBounds[["alphaStart"]],
                 Z = paramBounds[["Zstart"]])
  fit <- bbmle::mle2(
    response ~ dnbinom(mu = MeanLoad(L1, L2, alpha, HI),
                       size = SizeNegBin(A1, A2, Z, HI)),
    data = data,
    start = start,
    lower = c(L1 = paramBounds[["L1LB"]],
              L2 = paramBounds[["L2LB"]],
              A1 = paramBounds[["A1LB"]],
              A2 = paramBounds[["A2LB"]],
              alpha = paramBounds[["alphaLB"]],
              Z = paramBounds[["ZLB"]]),
    upper = c(L1 = paramBounds[["L1UB"]],
              L2 = paramBounds[["L2UB"]],
              A1 = paramBounds[["A1UB"]],
              A2 = paramBounds[["A2UB"]],
              alpha = paramBounds[["alphaUB"]],
              Z = paramBounds[["ZUB"]]),
    optimizer = config$optimizer,
    method = config$method,
    control = config$control)
  printConvergence(fit)
  return(fit)
}

defaultConfig <- list(optimizer = "optimx",
                      method = c("L-BFGS-B", "bobyqa"),
                      control = list(follow.on = TRUE))

model = "negbin"
data = pinwormsdata_bal
response = "Aspiculuris.Syphacia+1"
fixedAlphaParamBounds <- parasiteLoad::getParamBounds(model, data, response)
fixedAlphaParamBounds[["alphaStart"]] <- 1.39
fixedAlphaParamBounds[["alphaLB"]] <- 1.38
fixedAlphaParamBounds[["alphaUB"]] <- 1.40

fit_alphafixed <- FitAdvancedAlphaNegbin(data, response, hybridIndex = "HI",
                               paramBounds = fixedAlphaParamBounds, 
                               config = defaultConfig)

fit_alphavaries <- FitAdvancedAlphaNegbin(data, response, hybridIndex = "HI",
                                          paramBounds = parasiteLoad::getParamBounds(model, data, response),
                                          config = defaultConfig)

fit_noalpha <- FitAdvancedNoAlphaNegbin(data, response, hybridIndex = "HI",
                                          paramBounds = parasiteLoad::getParamBounds(model, data, response),
                                          config = defaultConfig)

Gtest <- function(fitA, fitB){
  LL1 <- logLik(fitA)
  LL2 <- logLik(fitB)
  dLL <- LL1 - LL2
  dDF <- 1 # alpha fixed in one case
  pvalue <- 1 - pchisq(2*dLL, df=dDF)
  chisqvalue <- qchisq(p = pvalue, df = dDF)
  data.frame(dLL = round(dLL, 2), dDF = dDF, pvalue = pvalue, chisqvalue = chisqvalue)
}

Gtest(fit_alphavaries, fit_noalpha)
Gtest(fit_alphavaries, fit_alphafixed)
