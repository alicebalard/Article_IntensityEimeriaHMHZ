source("DataPreparation.R")

##### Supp. material A ##### 

getAgeCat <- function(d){
  d$ageGroup <- "Old"
  d$ageGroup[d$Body_weight<=21] <- "Adult"
  d$ageGroup[d$Body_weight<=17] <- "Mature"
  d$ageGroup[d$Body_weight<=13] <- "Juvenile"
  d$ageGroup[d$Body_weight<=9] <- "Young"
  d$ageGroup[d$Body_weight<5] <- "Weanlings"
  d$ageGroup <- factor(d$ageGroup, levels = c("Weanlings", "Young","Juvenile", "Mature", "Adult", "Old"))
  return(d)
}

d <- getAgeCat(pinwormsdata_bal)
t <- as.data.frame.matrix(table(d$ageGroup, d$Aspiculuris_Syphacia>0))
t$sum <- rowSums(t)
t$prevalence <- round(t$"TRUE"/t$sum, 2)
t$ageGroup <- rownames(t)
t$x <- 1:6
t$lab <- paste0(as.character(100 * t$prevalence),"% (N=", as.character(t$sum), ")")
t$ageGroup <- factor(t$ageGroup, levels = c("Weanlings", "Young","Juvenile", "Mature", "Adult", "Old"))
t
test <- t[c(1,2)]
# Chi-square to test the hypothesis whether the number of infected is independent of age Cat
chisq.test(test)
# remove Weanlings (to few to test)
chisq.test(test[-1,])

d <- getAgeCat(qpcrdata)
t <- as.data.frame.matrix(table(d$ageGroup,  d$`delta_ct_max_MminusE+5`>0))
t$sum <- rowSums(t)
t$prevalence <- round(t$"TRUE"/t$sum, 2)
t$ageGroup <- rownames(t)
t$x <- 1:6
t$lab <- paste0(as.character(100 * t$prevalence),"% (N=", as.character(t$sum), ")")
t$ageGroup <- factor(t$ageGroup, levels = c("Weanlings", "Young","Juvenile", "Mature", "Adult", "Old"))
ggplot(t) + geom_bar(aes(x = ageGroup, y = prevalence), stat = "identity") + 
  geom_text(aes(x=x,y=prevalence,label=lab),vjust=-.2) +
  ggtitle("Eimeria prevalence")
t
test <- t[c(1,2)]
# Chi-square to test the hypothesis whether the number of infected is independent of age Cat
chisq.test(test)
# remove Weanlings (to few to test)
chisq.test(test[-1,])

# And between sexes?
t <- as.data.frame.matrix(table(pinwormsdata_bal$Sex, pinwormsdata_bal$Aspiculuris_Syphacia>0))
t$sum <- rowSums(t)
t$prevalence <- round(t$"TRUE"/t$sum, 2)
t$sex <- rownames(t)
t$x <- 1:2
t$lab <- paste0(as.character(100 * t$prevalence),"% (N=", as.character(t$sum), ")")
ggplot(t) + geom_bar(aes(x = sex, y = prevalence), stat = "identity") + 
  geom_text(aes(x=x,y=prevalence,label=lab),vjust=-.2) +
  ggtitle("Pinworms prevalence")
test <- t[c(1,2)]
# Chi-square to test the hypothesis whether the number of infected is independent of age Cat
chisq.test(test)
t

t <- as.data.frame.matrix(table(qpcrdata$Sex, qpcrdata$`delta_ct_max_MminusE+5`>0))
t$sum <- rowSums(t)
t$prevalence <- round(t$"TRUE"/t$sum, 2)
t$sex <- rownames(t)
t$x <- 1:2
t$lab <- paste0(as.character(100 * t$prevalence),"% (N=", as.character(t$sum), ")")
t
ggplot(t) + geom_bar(aes(x = sex, y = prevalence), stat = "identity") + 
  geom_text(aes(x=x,y=prevalence,label=lab),vjust=-.2) +
  ggtitle("Eimeria prevalence")
test <- t[c(1,2)]
# Chi-square to test the hypothesis whether the number of infected is independent of age Cat
chisq.test(test)

nrow(cleanData[cleanData$UsedForEimeriaRes %in% "yes",])
nrow(cleanData[cleanData$UsedForEimeriaImpactHealth %in% "yes",])
nrow(cleanData[cleanData$UsedForPinwormsRes %in% "yes",])
nrow(cleanData[cleanData$UsedForPinwormsImpactHealth %in% "yes",])

##### Supp. material B ##### 
##### Hybrid markers
getNloci <- function(df){
  df$HI_Num_loci <- as.numeric(gsub("HI ", "", df$HI_NLoci))
  return(df)
}

df <- getNloci(qpcrdata)
A <- ggplot(df, aes(HI_Num_loci)) + geom_histogram() + theme_bw() +
  ggtitle(paste("Mice tested for Eimeria prevalence, N=", nrow(df)))+ xlab("number of alleles")
A2 <- ggplot(df, aes(HI,HI_Num_loci)) + geom_point() + geom_smooth() + theme_bw() +
  ylab("number of alleles") + ggtitle("Mice tested for Eimeria prevalence")

df <- getNloci(body_data_eimeria)
B <-ggplot(df, aes(HI_Num_loci)) + geom_histogram() + theme_bw() +
  ggtitle(paste("Mice tested for body condition & Eimeria, N=", nrow(df)))+ xlab("number of alleles")
B2 <- ggplot(df, aes(HI,HI_Num_loci)) + geom_point() + geom_smooth() + theme_bw() +
  ylab("number of alleles") + ggtitle("Mice tested for body condition & Eimeria")

df <- getNloci(pinwormsdata_bal)
C <- ggplot(df, aes(HI_Num_loci)) + geom_histogram() + theme_bw() +
  ggtitle(paste("Mice tested for pinworms prevalence, N=", nrow(df))) + xlab("number of alleles")
C2 <- ggplot(df, aes(HI,HI_Num_loci)) + geom_point() + geom_smooth() + theme_bw() +
  ylab("number of alleles") + ggtitle("Mice tested for pinworms prevalence")

df <- getNloci(body_data_pinworms)

D <- ggplot(df, aes(HI_Num_loci)) + geom_histogram() + theme_bw() +
  ggtitle(paste("Mice tested for body condition & pinworms, N=", nrow(df)))+ xlab("number of alleles")
D2 <- ggplot(df, aes(HI,HI_Num_loci)) + geom_point() + geom_smooth() + theme_bw() +
  ylab("number of alleles") + ggtitle("Mice tested for body condition & pinworms")

FigSup1 <- plot_grid(A, B, C, D, nrow = 2)
FigSup1

FigSup2 <- plot_grid(A2, B2, C2, D2, nrow = 2)
FigSup2

supplementaryB <- plot_grid(FigSup1, FigSup2, ncol = 1)

# pdf(file = "../figures/supplementaryB.pdf",
#     width = 10, height = 10)
supplementaryB
# dev.off()

##### Supp. material C ##### 
### See "Models_parasiteLoads.R"
### After reveiw JEB
## Separate Eimeria species (Supplementary material)

## --> E.ferrisi found in cecum
sub_qpcr_intensity_data <- qpcr_intensity_data[
  grep("ferrisi", qpcr_intensity_data$eimeriaSpecies),]

table(qpcrdata$presence_eferrisi_identified, qpcrdata$presence_eimeria_tissues)
nrow(sub_qpcr_intensity_data) # only 44 out of 70 identified at the species level

# again, start with 6, and find out which value of shift maximize the likelihood
fit_sub_qpcr_intensity_findShift <- parasiteLoad::analyse(
  data = sub_qpcr_intensity_data,
  response = "delta_ct_max_MminusE+6",
  model = "weibullshifted",
  group = "Sex")

# Get the shift optimal for H0
shift2 <- coef(fit_sub_qpcr_intensity_findShift$H0)["SHIFT"]

# Now run while adding this 0.93 value
sub_qpcr_intensity_data$shiftedDeltaCt <-
  sub_qpcr_intensity_data$`delta_ct_max_MminusE+6` + shift2

fit_suqpcr_intensity <- parasiteLoad::analyse(
  data = sub_qpcr_intensity_data,
  response = "shiftedDeltaCt",
  model = "weibull",
  group = "Sex")

fit_suqpcr_intensity
coef(fit_suqpcr_intensity$H0)
# L1     alpha   myshape 
# 6.5429342 0.7900411 2.5840015 

# Get actual values of L1 and L2 for raw deltaCt
getUnshiftedData2 <- function(vector){
  vector[c("L1", "L2")] <- vector[c("L1", "L2")] - 6 - shift2
  return(vector)
}
getUnshiftedData2(coef(fit_suqpcr_intensity$H0))
getUnshiftedData2(coef(fit_suqpcr_intensity$H1))
getUnshiftedData2(coef(fit_suqpcr_intensity$H2$groupA))
getUnshiftedData2(coef(fit_suqpcr_intensity$H2$groupB))
getUnshiftedData2(coef(fit_suqpcr_intensity$H3$groupA))
getUnshiftedData2(coef(fit_suqpcr_intensity$H3$groupB))

plotInfEimeria_sub <- parasiteLoad::bananaPlot(
  mod = fit_suqpcr_intensity$H0,
  data = sub_qpcr_intensity_data,
  response = "shiftedDeltaCt",
  islog10 = FALSE, group = "Sex",
  cols = c("white", "white")) +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        legend.position = "none") +
  scale_y_continuous(
    breaks = c(2:12 + 0.14), # to remove the shift
    labels = as.character(round(c(2:12 + shift2) - (6 + shift2))),
    expression(
      paste("Eimeria ferrisi intensity (", Delta, "Ct Mouse - Eimeria)")))
plotInfEimeria_sub 

## prevalence with only Eimeria ferrisi identified

## Apply on our datasets
C <- prevFun(qpcrdata, "presence_eferrisi_identified")
prevalenceAlongHIEimeriaFerrisi <- C$plot + ylab("Predicted probability of infection") + 
   theme(legend.position = "none")
C$output # no signif of HI, sex, or both

# pdf("../figures/supplementaryC.pdf", width = 12, height = 5)
plot_grid(prevalenceAlongHIEimeriaFerrisi, plotInfEimeria_sub, labels = c("a", "b"))
# dev.off() 
##### Supp. material D ##### 

## histograms of all detected nematodes (after removing BALdata embrios)
nrow(pinwormsdata_bal)
apply(pinwormsdata_bal[,listWorms], 2, table)

Ntot <- nrow(pinwormsdata_bal)
makeHistWorm <- function(i){
  Ns <- table(pinwormsdata_bal[listWorms[i]]>0)[2]
  xpos <- max(pinwormsdata_bal[listWorms[i]], na.rm = T)/2
  p <- ggplot(pinwormsdata_bal, aes_string(x = listWorms[i])) + 
    geom_histogram(col = "black") + theme_bw() + 
    annotate("text", x = xpos, y= 30, label = paste0(Ns, " positive (N = ", Ntot, ")")) +
    scale_y_log10()
  return(p)
}

plot_grid(makeHistWorm(1), makeHistWorm(2), makeHistWorm(3),
          makeHistWorm(4), makeHistWorm(5), makeHistWorm(6),
          labels=c("A", "B", "C", "D", "E", "F"), label_size=15, ncol=2)

##### Supp. material E ##### 
## Choice of distribution, see choose_distribution.R

##### Supp. material F ##### 
length(intersect(body_data_eimeria$Mouse_ID, body_data_pinworms$Mouse_ID))
length(body_data_eimeria$Mouse_ID)
length(body_data_pinworms$Mouse_ID)
  # body_data_pinworms include all body_data_eimeria (N=333) + 123 mice (previous years)

x <- body_data_pinworms$Body_weight
length(x)

findGoodDist(x, distribs = c("normal", "t", "weibull"), 
             distribs2 = c("norm", "t", "weibull")) # Normal performs well

fitMortality <- parasiteLoad::analyse(
  data = body_data_pinworms,
  response = "Body_weight",
  model = "normal",
  group = "Sex")

plotMortality <- parasiteLoad::bananaPlot(
  mod = fitMortality$H1,
  data = body_data_pinworms,
  response = "Body_weight",
  islog10 = FALSE, group = "Sex",
  cols = c("white", "white")) 

P2 <- plotMortality +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("HI (hybrid index)") +
  ylab("Body weight (g) as approximation of age")

##### Supp. material G ##### 
## See prevalenceAlongHI.R

##### Supp. material H ##### 

# Eimeria+pinworms vs Eimeria without pinworms
table(qpcr_intensity_data$Aspiculuris_Syphacia > 0)
