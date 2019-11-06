source("DataPreparation.R")
# probability of infection

# For logistic regression, we use 0 -> 0.5 to test if linear decrease towards the center

# Fun1: data preparation
getDataPrep <- function(data, presence){
  data$HIdistToCenter <- data$HI
  data$HIdistToCenter[data$HI > 0.5] <-
    1-data$HIdistToCenter[data$HI > 0.5]
  data$presence <- as.numeric(as.character(data[,presence]))
  data$Sex <- factor(data$Sex)
  return(data)
}

# Fun2: model
getModel <- function(data){
  mylogit <- glm(presence ~ HIdistToCenter * Sex, data = data, family = "binomial")
  return(mylogit)
}

# Fun3: predict and plot
getPlot <- function(data, mylogit){
  newdata1 <- with(data,
                   data.frame(HIdistToCenter = mean(HIdistToCenter), 
                              Sex = factor(c("M","F"))))
  newdata1$SexP <- predict(mylogit, newdata = newdata1, type = "response")
  newdata2 <- with(data,
                   data.frame(HIdistToCenter = rep(seq(from = 0, to = 0.5, length.out = 100), 2),
                              Sex = factor(rep(c("M", "F"), each = 100))))
  newdata3 <- cbind(newdata2, predict(mylogit, newdata = newdata2, type = "link",
                                      se = TRUE))
  newdata3 <- within(newdata3, {
    PredictedProb <- plogis(fit)
    LL <- plogis(fit - (1.96 * se.fit))
    UL <- plogis(fit + (1.96 * se.fit))
  })
  myplot <- ggplot(newdata3, aes(x = HIdistToCenter, y = PredictedProb)) +
    geom_ribbon(aes(ymin = LL,
                    ymax = UL, fill = Sex), alpha = 0.2) +
    geom_line(aes(colour = Sex), size = 1) +
    scale_y_continuous(limits = c(0, 1),
                       breaks = seq(0, 1, by = 0.1),
                       labels = c("0", "10%", "20%", "30%", "40%", "50%",
                                  "60%", "70%", "80%", "90%", "100%")) +
    geom_rug(data = data, aes(x = HIdistToCenter, y = presence), alpha=.5) +
    scale_color_manual(values = c("#E69F00", "#009E73")) +
    scale_fill_manual(values = c("#E69F00", "#009E73")) +
    theme_bw() +
    # for adding average line
    geom_hline(yintercept = mean(data$presence), 
               col = "darkgrey", linetype = 2, size = 1)
  return(myplot)
}

# Full analysis:
prevFun <- function(data, presence){
  # prepare data: add HIdistToCenter
  prepData <- getDataPrep(data, presence)
  # Logistic regression, also called a logit model, is used to model dichotomous
  # outcome variables. In the logit model the log odds of the outcome is modeled
  # # as a linear combination of the predictor variables
  mylogit <- getModel(prepData)
  myplot <- getPlot(prepData, mylogit)
  results <- list(plot = myplot, output = summary(mylogit))
  return(results)
}

## Apply on our datasets
A <- prevFun(qpcrdata, "presence_eimeria_tissues")
prevalenceAlongHIEimeria <- A$plot
A$output

B <- prevFun(pinwormsdata_bal, "presence_oxyurids")
prevalenceAlongHIpinworms <- B$plot
B$output

### For supplementary data: on half transects
A_Mmd <- prevFun(qpcrdata[qpcrdata$HI < 0.5,], "presence_eimeria_tissues")
prevalenceAlongHIEimeria_Mmd <- A_Mmd$plot
A_Mmd$output

A_Mmm <- prevFun(qpcrdata[qpcrdata$HI >= 0.5,], "presence_eimeria_tissues")
prevalenceAlongHIEimeria_Mmm <- A_Mmm$plot
A_Mmm$output

B_Mmd <- prevFun(pinwormsdata_bal[pinwormsdata_bal$HI < 0.5,], "presence_oxyurids")
prevalenceAlongHIpinworms <- B_Mmd$plot
B_Mmd$output

B_Mmm <- prevFun(pinwormsdata_bal[pinwormsdata_bal$HI >= 0.5,], "presence_oxyurids")
prevalenceAlongHIpinworms <- B_Mmm$plot
B_Mmm$output

### Model with a new factor, "side", for Efal and Efer
data2 <- qpcrdata
data2$presence_efalci_identified <- 0
data2[grep("falci", data2$eimeriaSpecies),"presence_efalci_identified"] <- 1

data2$side[data2$HI < 0.5] <- "Dom"
data2$side[data2$HI >= 0.5] <- "Mus"

table(data2$presence_efalci_identified, data2$side)
table(data2$presence_eferrisi_identified, data2$side)

modEfal <- glm(presence ~ HIdistToCenter * Sex * side, data = getDataPrep(data2, "presence_efalci_identified"), family = "binomial")
summary(modEfal)
modEfalNoSex <- glm(presence ~ HIdistToCenter * side, data = getDataPrep(data2, "presence_efalci_identified"), family = "binomial")
summary(modEfalNoSex)

modEfer <- glm(presence ~ HIdistToCenter * Sex * side, data = getDataPrep(data2, "presence_eferrisi_identified"), family = "binomial")
summary(modEfer)
modEferNoSex <- glm(presence ~ HIdistToCenter * side, data = getDataPrep(data2, "presence_eferrisi_identified"), family = "binomial")
summary(modEferNoSex)

## Model simply along HI
table(qpcrdata$eimeriaSpecies)
table(qpcr_intensity_data$eimeriaSpecies) # only > -5 kept

ggplot(qpcrdata, aes(delta_ct_ilwe_MminusE, delta_ct_cewe_MminusE, col = eimeriaSpecies)) +
  geom_point(size = 5, alpha = .7) + theme_bw() +
  geom_hline(aes(yintercept = -5))+
  geom_vline(aes(xintercept = -5))

# prevalence between both sides 
library('WVPlots')
data2plot <- qpcrdata
data2plot$Efer <- grepl("ferrisi", qpcrdata$eimeriaSpecies)
data2plot$Efal <- grepl("falcif", qpcrdata$eimeriaSpecies)
data2plot$Ever <- grepl("vermi", qpcrdata$eimeriaSpecies)
#quick visu
P1 <- WVPlots::BinaryYScatterPlot(data2plot, "HI", "Efer", title="Prevalence") +
  theme_bw() + scale_y_continuous(breaks=seq(0, 1, 0.1))
P2 <- WVPlots::BinaryYScatterPlot(data2plot, "HI", "Efal", title="Prevalence") +
  theme_bw() + scale_y_continuous(breaks=seq(0, 1, 0.1))
P3 <- WVPlots::BinaryYScatterPlot(data2plot, "HI", "Ever", title="Prevalence") +
  theme_bw() + scale_y_continuous(breaks=seq(0, 1, 0.1))
plot_grid(P1, P2, P3)

# test -> nothing significant
modA2 <- glm(grepl("ferrisi", qpcrdata$eimeriaSpecies) ~ HI, data = qpcrdata, family = "binomial")
summary(modA2)
modA2

modA3 <- glm(grepl("falci", qpcrdata$eimeriaSpecies) ~ HI, data = qpcrdata, family = "binomial")
summary(modA3)
modA3

# Final plot of model of Efer to have the visual of the model
newdata1 <- with(qpcrdata, data.frame(HI = mean(HI)))
newdata2 <- with(qpcrdata, data.frame(HI = rep(seq(from = 0, to = 1, length.out = 100), 2)))
newdata3 <- cbind(newdata2, predict(modA3, newdata = newdata2, type = "link", se = TRUE))
newdata3 <- within(newdata3, {
  PredictedProb <- plogis(fit)
  LL <- plogis(fit - (1.96 * se.fit))
  UL <- plogis(fit + (1.96 * se.fit))
})
myplot <- ggplot(newdata3, aes(x = HI, y = PredictedProb)) +
  geom_ribbon(aes(ymin = LL,
                  ymax = UL), alpha = 0.2) +
  geom_line() +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.1),
                     labels = c("0", "10%", "20%", "30%", "40%", "50%",
                                "60%", "70%", "80%", "90%", "100%")) +
  geom_rug(data = qpcrdata, aes(x = HI, y = as.integer(grepl("ferrisi", qpcrdata$eimeriaSpecies))), alpha=.5) +
  theme_bw() +
  # for adding average line
  geom_hline(yintercept = mean(grepl("falci", qpcrdata$eimeriaSpecies)),
             col = "darkgrey", linetype = 2, size = 1)
myplot

  