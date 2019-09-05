source("DataPreparation.R")

#  Body residuals
# Hyp : hybrid less affected than pure strains by eimeria infection.

# body_data_eimeria

# improve paramBounds
speparam <- c(L1start = 0,
              L1LB = -20,
              L1UB = 20,
              L2start = 0,
              L2LB = -20,
              L2UB = 20,
              alphaStart = 0, alphaLB = -100, alphaUB = 100,
              mysdStart = 1, mysdLB = 1.e-9, mysdUB = 1000)
min(stats::na.omit(body_data_eimeria[["residuals"]]))
hist(body_data_eimeria[["residuals"]])

##All
fitResiduals_eimeria <-
  parasiteLoad::analyse(data = body_data_eimeria,
                        response = "residuals",
                        model = "normal",
                        group = "presence_eimeria_tissues",
                        myparamBounds = speparam)

# dLL dDF    pvalue
# 1   0   1 0.9999996
# [1] "Testing H1 no alpha vs alpha"
# dLL dDF    pvalue
# 1 0.34   1 0.4086401
# [1] "Testing H2 groupA no alpha vs alpha"
# dLL dDF pvalue
# 1   0   1      1
# [1] "Testing H2 groupB no alpha vs alpha"
# dLL dDF    pvalue
# 1 0.23   1 0.4989435
# [1] "Testing H3 groupA no alpha vs alpha"
# dLL dDF    pvalue
# 1 0.15   1 0.5787693
# [1] "Testing H3 groupB no alpha vs alpha"
# dLL dDF    pvalue
# 1 0.01   1 0.8991108
# [1] "Testing H1 vs H0"
# dLL dDF    pvalue
# 1   3   1 0.0142466
# [1] "Testing H2 vs H0"
# dLL dDF    pvalue
# 1 1.08   3 0.5386173
# [1] "Testing H3 vs H1"
# dLL dDF    pvalue
# 1 0.79   4 0.8139575
# [1] "Testing H3 vs H2"
# dLL dDF     pvalue
# 1 2.71   2 0.06679559

plotResiduals_eimeria <-
  bananaPlot(mod = fitResiduals_eimeria$H3,
             data = body_data_eimeria,
             response = "residuals",
             group = "presence_eimeria_tissues") +
  scale_fill_manual(values = c("grey", "green")) +
  scale_color_manual(values = c("grey", "green")) +
  theme_bw()
plotResiduals_eimeria

# pinworms

# body_data_pinworms


##All
fitResiduals_pinworms <-
    parasiteLoad::analyse(data = body_data_pinworms,
                        response = "residuals",
                        model = "normal",
                        group = "presence_oxyurids",
                        myparamBounds = speparam)

# [1] "Testing H0 no alpha vs alpha"
# dLL dDF pvalue
# 1   0   1      1
# [1] "Testing H1 no alpha vs alpha"
# dLL dDF    pvalue
# 1 0.01   1 0.9062973
# [1] "Testing H2 groupA no alpha vs alpha"
# dLL dDF    pvalue
# 1   0   1 0.9999995
# [1] "Testing H2 groupB no alpha vs alpha"
# dLL dDF    pvalue
# 1   0   1 0.9999997
# [1] "Testing H3 groupA no alpha vs alpha"
# dLL dDF    pvalue
# 1 0.3   1 0.4389196
# [1] "Testing H3 groupB no alpha vs alpha"
# dLL dDF    pvalue
# 1   0   1 0.9575637
# [1] "Testing H1 vs H0"
# dLL dDF       pvalue
# 1 6.74   1 0.0002418417
# [1] "Testing H2 vs H0"
# dLL dDF   pvalue
# 1   0   3 0.999871
# [1] "Testing H3 vs H1"
# dLL dDF    pvalue
# 1 0.64   4 0.8662367
# [1] "Testing H3 vs H2"
# dLL dDF       pvalue
# 1 7.37   2 0.0006300432

plotResiduals_pinworms <-
  bananaPlot(mod = fitResiduals_pinworms$H3,
             data = body_data_pinworms,
             response = "residuals",
             group = "presence_oxyurids") +
  scale_fill_manual(values = c("grey", "green")) +
  scale_color_manual(values = c("grey", "green")) +
  theme_bw()
plotResiduals_pinworms
