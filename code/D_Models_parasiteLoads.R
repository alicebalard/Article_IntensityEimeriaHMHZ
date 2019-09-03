source("DataPreparation.R")

fit_WATWM_intensity <- parasiteLoad::analyse(
  data = pinwormsdata_watwm[pinwormsdata_watwm$Aspiculuris.Syphacia >=1,],
  response = "Aspiculuris.Syphacia", group = "Sex", model = "negbin")

coef(fit_WATWM_intensity$H1)

# Pinworms (A. tetraptera and S. obvelata) : BAL

fit_BAL_abundance <- parasiteLoad::analyse(data = pinwormsdata_bal,
                             response = "Aspiculuris.Syphacia+1",
                             model = "negbin", group = "Sex")

coef(fit_BAL_abundance$H0)
coef(fit_BAL_abundance$H1)
coef(fit_BAL_abundance$H2$groupA)
coef(fit_BAL_abundance$H2$groupB)
coef(fit_BAL_abundance$H3$groupA)
coef(fit_BAL_abundance$H3$groupB)

plotAllWorms <- parasiteLoad::bananaPlot(mod = fit_BAL_abundance$H3,
                           data = pinwormsdata_bal,
                           response = "Aspiculuris.Syphacia+1",
                           islog10 = TRUE, group = "Sex",
                           cols = c("#E69F00", "#009E73"))
plotAllWorms

fit_BAL_intensity <- parasiteLoad::analyse(
  data = pinwormsdata_bal[pinwormsdata_bal$Aspiculuris_Syphacia >=1,],
  response = "Aspiculuris_Syphacia",
  model = "negbin", group = "Sex")

coef(fit_BAL_intensity$H0)
coef(fit_BAL_intensity$H1)
coef(fit_BAL_intensity$H2$groupA)
coef(fit_BAL_intensity$H2$groupB)
coef(fit_BAL_intensity$H3$groupA)
coef(fit_BAL_intensity$H3$groupB)

plotInfWorms <- parasiteLoad::bananaPlot(
  mod = fit_BAL_intensity$H3,
  data = pinwormsdata_bal[pinwormsdata_bal$Aspiculuris_Syphacia >=1,],
  response = "Aspiculuris.Syphacia+1",
  islog10 = TRUE, group = "Sex",
  cols = c("#E69F00", "#009E73")) +
  theme_bw() 

breaks = c(2,11,31,101,301)
labels = as.character(breaks -1)

plotInfWorms <- plotInfWorms +
  scale_y_log10(
    breaks = breaks,
    labels = labels,
    paste("Pinworm count"))

# qPCR Eimeria

# intensity : delta_ct_max_MminusE+5 to have it positive (in DataPReparation)
# qpcr_intensity_data <- qpcrdata[qpcrdata$`delta_ct_max_MminusE+5` > 0,]

# start with 6, and find out which value of shift maximize the likelihood
fit_qpcr_intensity_findShift <- parasiteLoad::analyse(
  data = qpcr_intensity_data,
  response = "delta_ct_max_MminusE+6",
  model = "weibullshifted",
  group = "Sex")

# Get the shift optimal for H0
shift <- coef(fit_qpcr_intensity_findShift$H0)["SHIFT"]
shift
# fit_qpcr_intensity_findShift$H0 : 1.142989

# Now run while adding this 1.14 value
qpcr_intensity_data$shiftedDeltaCt <-
  qpcr_intensity_data$`delta_ct_max_MminusE+6` + shift
    
fit_qpcr_intensity <- parasiteLoad::analyse(
  data = qpcr_intensity_data,
  response = "shiftedDeltaCt",
  model = "weibull",
  group = "Sex")

fit_qpcr_intensity
coef(fit_qpcr_intensity$H0)
# L1     alpha   myshape
# 6.4376774 0.7391913 2.3332110

# Get actual values of L1 and L2 for raw deltaCt
getUnshiftedData <- function(vector){
  vector[c("L1", "L2")] <- vector[c("L1", "L2")] - 6 - shift
  return(vector)
}
getUnshiftedData(coef(fit_qpcr_intensity$H0))
getUnshiftedData(coef(fit_qpcr_intensity$H1))
getUnshiftedData(coef(fit_qpcr_intensity$H2$groupA))
getUnshiftedData(coef(fit_qpcr_intensity$H2$groupB))
getUnshiftedData(coef(fit_qpcr_intensity$H3$groupA))
getUnshiftedData(coef(fit_qpcr_intensity$H3$groupB))

plotInfEimeria <- parasiteLoad::bananaPlot(
  mod = fit_qpcr_intensity$H0,
  data = qpcr_intensity_data,
  response = "shiftedDeltaCt",
  islog10 = FALSE, group = "Sex",
  cols = c("white", "white")) +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        legend.position = "none") +
  scale_y_continuous(
    breaks = c(2:12 + shift), # to remove the shift
    labels = as.character(round(c(2:12 + shift) - (6 + shift))),
    expression(
      paste("Eimeria intensity (", Delta, "Ct Mouse - Eimeria)")))
plotInfEimeria
    
