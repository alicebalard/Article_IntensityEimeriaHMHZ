source("DataPreparation.R")

#### Investigate coinfection pinworms-Eimeria

# qpcrdata contains mice that have been investigated for both worms and eimeria
coinfDf <- qpcrdata[!is.na(rowSums(qpcrdata[c(l istWorms)])),] # N = 383
coinfDf$presence_pinworms <- as.numeric(coinfDf$Aspiculuris_Syphacia > 0)
coinfDf$presence_eimeria_tissues <- as.numeric(as.character(coinfDf$presence_eimeria_tissues))

table(pinworms = coinfDf$presence_pinworms)
211 / (172+211) # 55% - 52.5% in the complete dataset
table(Eimeria = coinfDf$presence_eimeria_tissues)
70/(70+313) # 18%

table(Eimeria = coinfDf$presence_eimeria_tissues, pinworms = coinfDf$presence_pinworms)

sum(table(Eimeria = coinfDf$presence_eimeria_tissues, pinworms = coinfDf$presence_pinworms))

chisq.test(coinfDf$presence_eimeria_tissues, coinfDf$presence_pinworms)

summary(lm(presence_eimeria_tissues ~ presence_pinworms, data = coinfDf))



###### Other stats, fail cause low power ###### 
## Analysis 1. Richness along HI
coinfDf$numberOfWormsGroup <- rowSums(coinfDf[c(listWorms)] > 0)

descdist(coinfDf$numberOfWormsGroup)
# we log transform to normalise
coinfDf$numberOfWormsGroupLOG <- log10(coinfDf$numberOfWormsGroup + 1)

ggplot(coinfDf, aes(HI,numberOfWormsGroup))+
  geom_point(col = "black", pch = 21, size = 4) + theme_bw() +
  geom_smooth()

fitNWorms <- parasiteLoad::analyse(data = coinfDf,
                                   response = "numberOfWormsGroupLOG",
                                   model = "normal", group = "Sex")
# [1] "Testing H0 no alpha vs alpha"
# dLL dDF    pvalue
# 1   0   1 0.9527644
# [1] "Testing H1 no alpha vs alpha"
# dLL dDF   pvalue
# 1 0.16   1 0.567228
# [1] "Testing H2 groupA no alpha vs alpha"
# dLL dDF    pvalue
# 1 0.11   1 0.6380954
# [1] "Testing H2 groupB no alpha vs alpha"
# dLL dDF    pvalue
# 1 0.21   1 0.5197232
# [1] "Testing H3 groupA no alpha vs alpha"
# dLL dDF    pvalue
# 1 0.81   1 0.2034286
# [1] "Testing H3 groupB no alpha vs alpha"
# dLL dDF    pvalue
# 1 0.11   1 0.6397515
# [1] "Testing H1 vs H0"
# dLL dDF       pvalue
# 1 5.84   1 0.0006293671
# [1] "Testing H2 vs H0"
# dLL dDF    pvalue
# 1 0.6   3 0.7530853
# [1] "Testing H3 vs H1"
# dLL dDF   pvalue
# 1 1.59   4 0.528784
# [1] "Testing H3 vs H2"
# dLL dDF      pvalue
# 1 6.83   2 0.001078924

fitNWorms$H1
# Coefficients:
#   L1        L2     alpha      mysd 
# 0.1854120 0.2811688 0.1573972 0.1854662
log10(4+1)

10^(fitNWorms$H1@coef["L1"]) -1 
# Mmd richness: 0.53
10^(fitNWorms$H1@coef["L2"]) -1
# Mmm richness: 0.91

plotNWorms <- parasiteLoad::bananaPlot(mod = fitNWorms$H1,
                                       data = coinfDf,
                                       response = "numberOfWormsGroupLOG",
                                       islog10 = F, group = "Sex",
                                       cols = c("#E69F00", "#009E73"))+ 
  scale_y_continuous(breaks = log10((0:5) +1), labels = as.character(0:5), 
                     name= "Worms richness")

plotNWorms 

####### Part 2: Eimeria and pinworms

coinfDf$presence_eimeria_tissues <- as.numeric(as.character(coinfDf$presence_eimeria_tissues))
coinfDf$presence_pinworms <- as.numeric(coinfDf$Aspiculuris_Syphacia > 0)

coinfDf$RichnessEimeriaPinworms <- coinfDf$presence_eimeria_tissues + coinfDf$presence_pinworms

hist(coinfDf$RichnessEimeriaPinworms)

ggplot(coinfDf, aes(HI,RichnessEimeriaPinworms))+
  geom_point(col = "black", pch = 21, size = 4) + theme_bw() +
  geom_smooth()

fitNEimeriaPinworms <- parasiteLoad::analyse(data = coinfDf,
                                   response = "RichnessEimeriaPinworms",
                                   model = "normal", group = "Sex")

# [1] "Testing H0 no alpha vs alpha"
# dLL dDF    pvalue
# 1 0.97   1 0.1638067
# [1] "Testing H1 no alpha vs alpha"
# dLL dDF    pvalue
# 1 1.11   1 0.1367972
# [1] "Testing H2 groupA no alpha vs alpha"
# dLL dDF    pvalue
# 1 0.69   1 0.2386343
# [1] "Testing H2 groupB no alpha vs alpha"
# dLL dDF    pvalue
# 1 0.28   1 0.4572745
# [1] "Testing H3 groupA no alpha vs alpha"
# dLL dDF    pvalue
# 1 1.12   1 0.1344645
# [1] "Testing H3 groupB no alpha vs alpha"
# dLL dDF    pvalue
# 1 0.21   1 0.5211931
# [1] "Testing H1 vs H0"
# dLL dDF    pvalue
# 1 0.17   1 0.5594216
# [1] "Testing H2 vs H0"
# dLL dDF    pvalue
# 1 0.26   3 0.9131108
# [1] "Testing H3 vs H1"
# dLL dDF    pvalue
# 1 1.12   4 0.6904749
# [1] "Testing H3 vs H2"
# dLL dDF   pvalue
# 1 1.03   2 0.356768
fitNEimeriaPinworms$H0

# Coefficients:
#   L1      mysd     alpha 
# 0.8080063 0.6503932 0.3726448 

# Mice richness: 0.8 on both sides

plotNEimeriaPinworms <- parasiteLoad::bananaPlot(mod = fitNEimeriaPinworms$H0,
                                       data = coinfDf,
                                       response = "RichnessEimeriaPinworms",
                                       islog10 = F, group = "Sex",
                                       cols = c("#E69F00", "#009E73")) 
  # scale_y_continuous(breaks = log10((0:5) +1), labels = as.character(0:5), 
                     # name= "Worms richness")

plotNEimeriaPinworms 

## Analysis 2. Effect of presence of Eimeria on presence of pinworms

table(Eimeria = coinfDf$presence_eimeria_tissues, Pinworms = coinfDf$presence_pinworms)

ggplot(coinfDf, aes(presence_eimeria_tissues, presence_pinworms)) +
  geom_point() + geom_smooth(method = "lm")

modE <- lm(presence_eimeria_tissues ~ presence_pinworms, data = coinfDf)
summary(modE)

## pinworms + eimeria
coinfDf$RichnessPandE <- coinfDf$presence_eimeria_tissues + coinfDf$presence_pinworms

hist(coinfDf$RichnessPandE)
descdist(coinfDf$RichnessPandE)

ggplot(coinfDf, aes(HI,RichnessPandE))+
  geom_point(col = "black", pch = 21, size = 4) + theme_bw() +
  geom_smooth()

fitRichnessPandE <- parasiteLoad::analyse(data = coinfDf,
                                     response = "RichnessPandE",
                                     model = "normal", group = "Sex")
fitRichnessPandE$H0
# Coefficients:
#   L1      mysd     alpha 
# 0.8080063 0.6503932 0.3726448
plotRichnessPandE <- parasiteLoad::bananaPlot(mod = fitRichnessPandE$H0,
                                         data = coinfDf,
                                         response = "RichnessPandE",
                                         islog10 = F, group = "Sex",
                                         cols = c("#E69F00", "#009E73")) 
plotRichnessPandE

### All together, all worms + eimeria?
coinfDf$fullRichness <- coinfDf$presence_eimeria_tissues + coinfDf$numberOfWormsGroup

hist(coinfDf$fullRichness)
descdist(coinfDf$fullRichness)

# we log transform to normalise
coinfDf$fullRichnessLOG <- log10(coinfDf$fullRichness + 1)
descdist(coinfDf$fullRichnessLOG)

ggplot(coinfDf, aes(HI,fullRichnessLOG))+
  geom_point(col = "black", pch = 21, size = 4) + theme_bw() +
  geom_smooth()

fitRichness <- parasiteLoad::analyse(data = coinfDf,
                                   response = "fullRichnessLOG",
                                   model = "normal", group = "Sex")
fitRichness
# Coefficients:
#   L1        L2     alpha      mysd 
# 0.2555224 0.3028713 0.2301157 0.1910755 
plotRichness <- parasiteLoad::bananaPlot(mod = fitRichness$H0,
                                                 data = coinfDf,
                                                 response = "fullRichnessLOG",
                                                 islog10 = F, group = "Sex",
                                                 cols = c("#E69F00", "#009E73")) 
plotRichness

#################################################
#### Eimeria in pinworms model, and vice versa###
#################################################
coinfDf$`Aspiculuris.Syphacia+1` <- coinfDf$Aspiculuris_Syphacia + 1
coinfDf$presence_eimeria_tissues <- as.factor(coinfDf$presence_eimeria_tissues)
coinfDf$presence_pinworms <- as.factor(coinfDf$presence_pinworms)

fit_pin_intensity_coinfEim <- parasiteLoad::analyse(
  data = coinfDf[coinfDf$Aspiculuris_Syphacia >=1,],
  response = "Aspiculuris_Syphacia",
  model = "negbin", group = "presence_eimeria_tissues")
# [1] "Testing H0 no alpha vs alpha"
# dLL dDF   pvalue
# 1 0.06   1 0.721007
# [1] "Testing H1 no alpha vs alpha"
# dLL dDF    pvalue
# 1 0.91   1 0.1783199
# [1] "Testing H2 groupA no alpha vs alpha"
# dLL dDF     pvalue
# 1 2.04   1 0.04320654
# [1] "Testing H2 groupB no alpha vs alpha"
# dLL dDF     pvalue
# 1 3.99   1 0.00473545
# [1] "Testing H3 groupA no alpha vs alpha"
# dLL dDF     pvalue
# 1 3.28   1 0.01040524
# [1] "Testing H3 groupB no alpha vs alpha"
# dLL dDF     pvalue
# 1 2.24   1 0.03427029
# [1] "Testing H1 vs H0"
# dLL dDF       pvalue
# 1 19.44   2 3.601098e-09
# [1] "Testing H2 vs H0"
# dLL dDF     pvalue
# 1 6.33   4 0.01305669
# [1] "Testing H3 vs H1"
# dLL dDF    pvalue
# 1 6.76   6 0.0355307
# [1] "Testing H3 vs H2"
# dLL dDF       pvalue
# 1 19.87   4 4.900146e-08

breaks = c(2,11,31,101,301)
labels = as.character(breaks -1)

plotpin_intensity_coinfEim <- parasiteLoad::bananaPlot(
  mod = fit_pin_intensity_coinfEim$H3,
  data = coinfDf[coinfDf$Aspiculuris_Syphacia >=1,],
  response = "Aspiculuris.Syphacia+1",
  islog10 = F, group = "presence_eimeria_tissues",
  cols = c("#E69F00", "#009E73")) +
  theme_bw() +
  scale_y_log10(
    breaks = breaks,
    labels = labels,
    paste("Pinworm count"))

plotpin_intensity_coinfEim

# qPCR Eimeria

# intensity : delta_ct_max_MminusE+5 to have it positive (in DataPReparation)
coinfDf_intensity <- coinfDf[coinfDf$`delta_ct_max_MminusE+5` > 0,]

# start with 6, and find out which value of shift maximize the likelihood
fit_qpcr_intensity_findShift_coinfPW <- parasiteLoad::analyse(
  data = coinfDf_intensity,
  response = "delta_ct_max_MminusE+6",
  model = "weibullshifted",
  group = "presence_pinworms")

# Get the shift optimal for H3 : both are very low, so we can accept 6
shift <- coef(fit_qpcr_intensity_findShift_coinfPW$H3$groupB)["SHIFT"]
shift

fit_qpcr_intensity_coinfPW <- parasiteLoad::analyse(
  data = coinfDf_intensity,
  response = "delta_ct_max_MminusE+6",
  model = "weibull",
  group = "presence_pinworms")
# [1] "Testing H0 no alpha vs alpha"
# dLL dDF     pvalue
# 1 2.47   1 0.02633725
# [1] "Testing H1 no alpha vs alpha"
# dLL dDF     pvalue
# 1 2.96   1 0.01492237
# [1] "Testing H2 groupA no alpha vs alpha"
# dLL dDF   pvalue
# 1 0.31   1 0.431386
# [1] "Testing H2 groupB no alpha vs alpha"
# dLL dDF      pvalue
# 1 4.57   1 0.002501114
# [1] "Testing H3 groupA no alpha vs alpha"
# dLL dDF    pvalue
# 1 0.27   1 0.4656953
# [1] "Testing H3 groupB no alpha vs alpha"
# dLL dDF      pvalue
# 1 4.93   1 0.001687092
# [1] "Testing H1 vs H0"
# dLL dDF    pvalue
# 1 0.5   1 0.3186214
# [1] "Testing H2 vs H0"
# dLL dDF     pvalue
# 1 3.38   3 0.07995699
# [1] "Testing H3 vs H1"
# dLL dDF    pvalue
# 1 3.29   4 0.1592259
# [1] "Testing H3 vs H2"
# dLL dDF    pvalue
# 1 0.41   2 0.6621981

fit_qpcr_intensity_coinfPW
coef(fit_qpcr_intensity_coinfPW$H2$groupA)
# L1      alpha    myshape 
# 3.5234749 -0.6754201  2.0245766 
coef(fit_qpcr_intensity_coinfPW$H2$groupB)
# L1    alpha  myshape 
# 6.253862 1.183998 1.866617 

plotInfEimeria_coinfPW <- parasiteLoad::bananaPlot(
  mod = fit_qpcr_intensity_coinfPW$H3,
  data = coinfDf_intensity,
  response = "delta_ct_max_MminusE+6",
  islog10 = FALSE, group = "presence_pinworms",
  cols = c("black", "red")) +
  theme_bw() #+
  theme(axis.title.x=element_blank(),
        legend.position = "none") #+
  # scale_y_continuous(
  #   breaks = c(2:12 + shift), # to remove the shift
  #   labels = as.character(round(c(2:12 + shift) - (6 + shift))),
  #   expression(
  #     paste("Eimeria intensity (", Delta, "Ct Mouse - Eimeria)")))
plotInfEimeria_coinfPW

### NB. Arbitrary cut. Just for validation to review. Not meaningful!
data2cut_eimeria <- qpcr_intensity_data

data2cut_eimeria$shiftedDeltaCt == data2cut_eimeria$delta_ct_max_MminusE + 6 + shift

mean(data2cut_eimeria$shiftedDeltaCt)
mean(data2cut_eimeria$delta_ct_max_MminusE) + 6 + shift

#arbitrary catagories HI: 0-0.2 / 0.2-0.6 / 0.6-1
data2cut_eimeria$hybridCat <- "hybrids"
data2cut_eimeria$hybridCat[data2cut_eimeria$HI <= 0.2] <- "Mmd"
data2cut_eimeria$hybridCat[data2cut_eimeria$HI >= 0.6] <- "Mmm"
data2cut_eimeria$hybridCat <- factor(data2cut_eimeria$hybridCat, levels = c("Mmd", "hybrids", "Mmm"))

table(data2cut_eimeria$hybridCat)

data2cut_eimeria %>% 
  dplyr::group_by(hybridCat) %>%
  dplyr::summarize(meanInt = mean(delta_ct_max_MminusE))

data2cut_eimeria %>% 
  dplyr::group_by(hybridCat) %>%
  dplyr::summarize(meanInt = mean(shiftedDeltaCt) - 6 - shift)

mod2 <- survreg(Surv(shiftedDeltaCt)~hybridCat, data = data2cut_eimeria, dist="weibull")
anova(mod2)

library(ggeffects)
dat <- ggpredict(mod2, terms = "hybridCat")

# to come back to delta Ct original
dat$predicted <- dat$predicted - 6 - shift
dat$conf.low <- dat$conf.low - 6 - shift
dat$conf.high <- dat$conf.high - 6 - shift

ggplot(dat, aes(x = x, y = predicted)) + 
  geom_point() +
  geom_errorbar(aes(ymin =conf.low, ymax= conf.high), width = 0.3) +
  scale_x_continuous(labels = c("Mmd", "", "hybrids", "", "Mmm"), name = "") +
  scale_y_continuous(name = "predicted deltaCt") +
  ggtitle("Eimeria intensity")



t <- table(rowSums(!is.na(cleanData[markersHI])), cleanData$Year)

sum(t)

