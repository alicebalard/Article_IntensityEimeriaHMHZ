source("DataPreparation.R")

# C. Choose distributions for each dataset

# Macro parasite counts -> negbin tested vs poisson and normal     
# qpcr -> explore best distribution


# The density plot and the CDF plot may be considered as the basic classical 
# goodness-of-fits plots. The Q-Q plot emphasizes the lack-of-fit at the 
# distribution tails while the P-P plot emphasizes the lack-of-fit at the 
# distribution center. The nbinom distribution could be prefered for its better
# description of the tail of the empirical distribution.
x <- pinwormsdata_bal$Aspiculuris_Syphacia
x <- x[x>0]
hist(x, breaks = 100)
descdist(x)
findGoodDist(x, distribs = c("normal", "negative binomial", "poisson"), 
             distribs2 = c("norm", "nbinom", "pois"))

x <- qpcrdata$`delta_ct_max_MminusE+5`
x <- x[x>0]
hist(x, breaks = 100)
descdist(x)
findGoodDist(x, distribs = c("normal", "t", "weibull"), 
             distribs2 = c("norm", "t", "weibull"))
