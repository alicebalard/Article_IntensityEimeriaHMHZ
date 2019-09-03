source("DataPreparation.R")
# probability of infection

# For logistic regression, we use 0 -> 0.5 to test if linear decrease towards the center
prevFun <- function(data, presence){
  data$HIdistToCenter <- data$HI
  data$HIdistToCenter[data$HI > 0.5] <-
    1-data$HIdistToCenter[data$HI > 0.5]
  data$pres <- as.numeric(as.character(data[,presence]))
  # Logistic regression, also called a logit model, is used to model dichotomous
  # outcome variables. In the logit model the log odds of the outcome is modeled
  # as a linear combination of the predictor variables
  mydata <- data.frame(HIdistToCenter = data$HIdistToCenter,
                       presence = data$pres,
                       Sex = data$Sex)

  mydata$Sex <- factor(mydata$Sex)
  mylogit <- glm(presence ~ HIdistToCenter * Sex, data = mydata, family = "binomial")
  newdata1 <- with(mydata,
                   data.frame(HIdistToCenter = mean(HIdistToCenter),
                              Sex = factor(c("M","F"))))
  newdata1$SexP <- predict(mylogit, newdata = newdata1, type = "response")
  newdata2 <- with(mydata,
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
    geom_line(aes(colour = Sex),
              size = 1) +
    scale_y_continuous(limits = c(0, 1),
                       breaks = seq(0, 1, by = 0.1),
                       labels = c("0", "10%", "20%", "30%", "40%", "50%",
                                  "60%", "70%", "80%", "90%", "100%")) +
    geom_rug(data = data, aes(x = HIdistToCenter, y = pres),
             alpha=.5) +
    scale_color_manual(values = c("#E69F00", "#009E73")) +
    scale_fill_manual(values = c("#E69F00", "#009E73")) +
    theme_bw() +
    # for adding average line
    geom_hline(yintercept = mean(data$pres), 
                 col = "darkgrey", linetype = 2, size = 1)
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

