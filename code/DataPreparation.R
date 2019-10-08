# Installation

## Packages
list.of.packages <- c("parasiteLoad",
                      "bbmle",
                      "devtools",
                      "optimx", # for bbmle it needs to be required(?)
                      "ggplot2",
                      "VennDiagram",
                      "fitdistrplus", # evaluate distribution
                      "epiR", # Sterne's exact method
                      "simpleboot", # BS
                      "plyr", # revalue and other
                      "ggmap",
                      "gridExtra",# several plots in one panel
                      "wesanderson", # nice colors
                      "cowplot",# several plots in one panel
                      "ggpubr")
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(list.of.packages)

## Reinstall the package in case I updated it
devtools::install_github("alicebalard/parasiteLoad@v2.0")
library(parasiteLoad)

## Install_github case
if(!"legendMap" %in% installed.packages()[,"Package"]){
  devtools::install_github("3wen/legendMap")
}
library(legendMap)

# Define function to be used to test, get the log lik and aic
tryDistrib <- function(x, distrib){
  # deals with fitdistr error:
  fit <- tryCatch(MASS::fitdistr(x, distrib), error=function(err) "fit failed")
  return(list(fit = fit,
              loglik = tryCatch(fit$loglik, error=function(err) "no loglik computed"), 
              AIC = tryCatch(fit$aic, error=function(err) "no aic computed")))
}

findGoodDist <- function(x, distribs, distribs2){
  l =lapply(distribs, function(i) tryDistrib(x, i))
  names(l) <- distribs
  print(l)
  listDistr <- lapply(distribs2, function(i){
    if (i %in% "t"){
      fitdistrplus::fitdist(x, i, start = list(df =2))
    } else {
      fitdistrplus::fitdist(x,i)
    }}
  ) 
  par(mfrow=c(2,2))
  denscomp(listDistr, legendtext=distribs2)
  cdfcomp(listDistr, legendtext=distribs2)
  qqcomp(listDistr, legendtext=distribs2)
  ppcomp(listDistr, legendtext=distribs2)
  par(mfrow=c(1,1))
}

## Define functions used for data analysis
area <- get_map(location = c(12, 51.5, 15, 53.5),
                source = "stamen", maptype = "toner-lite")

plotMap <- function(df){
  ggmap(area) +
    geom_point(data = df, shape = 21, size = 2,
               aes(Longitude, Latitude, fill = HI), alpha = .4) + # set up the points
    scale_fill_gradient("Hybrid\nindex", high="red",low="blue") +
    theme_bw() +
    geom_rect(xmin = 12, xmax = 12.7, ymin = 51.5, ymax = 51.9, fill = "white") +
    scale_bar(lon = 12.1, lat = 51.5, arrow_length = 10, arrow_distance = 20,
              distance_lon = 20, distance_lat = 7, distance_legend = 10,
              dist_unit = "km", orientation = TRUE, legend_size = 2,
              arrow_north_size = 4) +
    theme(legend.position = 'none', axis.ticks=element_blank())
}

myQuantitativeParasitology <- function(x){
  intensity <- round(median(x[x>0]),3)
  abundance <- round(median(x), 3)
  max <- max(x)
  Ni <- length(x)
  NiPos <- length(x[x>0])
  # Confidence intervals for prevalence calculated with Sterne's exact method
  sternetest <- epiR::epi.prev(pos = length(x[x > 0]), tested = length(x),
                               se = 1, sp=1, conf.level = .95, method = "sterne")
  cilow <- sternetest$ap["lower"]
  cihigh <- sternetest$ap["upper"]
  prevalence <- sternetest$ap["est"]
  ## Printout results
  Result <- cat(paste0("Prevalence % [CI 95%] (N infected hosts/ N hosts)\n",
                       round(prevalence,1), " [", round(cilow,1), "-", round(cihigh,1), "]",
                       " (", NiPos, "/", Ni,  ")\n",
                       "Abundance (Max parasite load)\n",
                       round(abundance,1), " (", max, ")\n",
                       "Intensity (Max parasite load)\n",
                       round(intensity,1), " (", max, ")"))
  return(Result)
}

## Prepare datasets for each analysis

# Load datasets from parasiteLoad
WATWMdata <- read.csv("https://raw.githubusercontent.com/alicebalard/parasiteLoad/master/data/WATWMdata.csv", na.strings = c("", " ", NA))
BALdata <- read.csv("https://raw.githubusercontent.com/alicebalard/parasiteLoad/master/data/BALdata.csv", na.strings = c("", " ", NA))

# Keep individuals with hybrid index and sex
WATWMdata <- WATWMdata[!is.na(WATWMdata$HI) & !is.na(WATWMdata$Sex),]
# pinworms "where are the wormy mice"
pinwormsdata_watwm <- WATWMdata[!is.na(WATWMdata$Aspiculuris.Syphacia),]
pinwormsdata_watwm$`Aspiculuris.Syphacia+1` <-
  pinwormsdata_watwm$Aspiculuris.Syphacia + 1
pinwormsdata_watwm$presence_oxyurids <- 1
pinwormsdata_watwm$presence_oxyurids[pinwormsdata_watwm$Aspiculuris.Syphacia == 0] <- 0

BALdata <- BALdata[!is.na(BALdata$HI) & !is.na(BALdata$Sex),]
BALdata$Status[BALdata$Status %in% "BA"] <- NA # error
BALdata$Status[is.na(BALdata$Status)] <- "adult" # NAs in the field data status were adults
BALdata$Status <- droplevels(BALdata$Status)

getinfotab <- function(df){
  return(list(Nmice = nrow(df),
              SexRatio = table(df$Sex),
              tableYear = table(df$Year),
              Nfarms = length(table(df$farm)),# define at 0.0001 degree
              meanAnimalperfarm = mean(table(df$farm)),
              medianAnimalperfarm = median(table(df$farm)),
              sdAnimalperfarm = qnorm(0.975)*sd(table(df$farm))/
                sqrt(sum(table(df$farm))),
              latrange = range(df$Latitude),
              lonrange = range(df$Longitude)))
}

nrow(WATWMdata)
table(WATWMdata$Sex)

## HERE PREPARE THE CLEAN TABLE THAT IS USED FOR EACH ANALYSIS OF THIS ARTICLE (sup table S1)
markersHI <- c("mtBamH", "YNPAR", "X332", "X347", "X65", "Tsx", "Btk", "Syap1", "Es1C", "Gpd1C", "Idh1C", "MpiC", "NpC", "Sod1C")
listWorms <- c("Aspiculuris_Syphacia", "Hymenolepis", "Taenia", "Trichuris", "Heterakis", "Mastophorus")
cleanData <- BALdata[c("Mouse_ID", "Sex", "Longitude", "Latitude", "Year", "farm", "Status",
                       markersHI, "HI_NLoci", "HI", listWorms, 
                       "Body_weight", "Body_length", "Tail_length", "Capture",
                       "delta_ct_ilwe_MminusE", "delta_ct_cewe_MminusE", "eimeriaSpecies")] 

cleanData$eimeriaSpecies <- as.character(cleanData$eimeriaSpecies)
cleanData$eimeriaSpecies[cleanData$eimeriaSpecies %in% c("Other", "Negative")] <- "no sp. identified"

# All 6 double infections were E.ferrisi in cecum and E.falciformis in ileum
cleanData$eimeriaSpecies[grep("Double", cleanData$eimeriaSpecies)] <- "E_ferrisi_cecum_E_vermiformis_ileum"

# Remove embryos N=7 mice used in no part of the study
embryos <- cleanData[grep("E", cleanData$Mouse_ID),"Mouse_ID"]
cleanData <- cleanData[!cleanData$Mouse_ID %in% embryos,]

# Verify the number of HI markers
cleanData$HI_NLoci <- as.numeric(gsub("HI ", "", cleanData$HI_NLoci))
table(cleanData$HI_NLoci == apply(cleanData, 1, function(x) sum(!is.na(x[markersHI]))))
cleanData[is.na(cleanData$HI_NLoci) |
            cleanData$HI_NLoci != apply(cleanData, 1, function(x) sum(!is.na(x[markersHI]))),]

# Remove 3 mice with few markers used in no part of the study (SK_2891, SK_3153-5, Sk3173)
cleanData <- cleanData[!cleanData$Mouse_ID %in% c("SK_2891", "SK_3153-5", "Sk3173"),]

# Correct the 2 wrong HI_NLoci (AA_0164, AA_0171)
cleanData$HI_NLoci <- apply(cleanData, 1, function(x) sum(!is.na(x[markersHI])))

# Indicate which mouse used in which part of the study: see further sections
cleanData$UsedForMap <- "no"
cleanData$UsedForEimeriaRes <- "no"
cleanData$UsedForPinwormsRes <- "no"
cleanData$UsedForEimeriaImpactHealth <- "no"
cleanData$UsedForPinwormsImpactHealth <- "no"

##### Geneland map

diploidMarkers <- c("Es1C", "Gpd1C", "Idh1C", "MpiC", "NpC", "Sod1C")

# use for map all individuals with 6 diploid markers

cleanData$UsedForMap <- rowSums(is.na(cleanData[diploidMarkers]))
cleanData$UsedForMap[cleanData$UsedForMap %in% 0] <- "yes"
cleanData$UsedForMap[cleanData$UsedForMap != "yes"] <- "no"

##### Eimeria qpcr ##### 
qpcrdata <- cleanData[!is.na(cleanData$delta_ct_cewe_MminusE) | !is.na(cleanData$delta_ct_ilwe_MminusE),]
df <- qpcrdata[, c("delta_ct_cewe_MminusE", "delta_ct_ilwe_MminusE")]
qpcrdata$delta_ct_max_MminusE <- apply(df, 1, function(x){max(x, na.rm = T)})
rm(df)
# threshold of detection by qPCR = -5. Then we add -5 to all to have positive values
qpcrdata$delta_ct_max_MminusE[qpcrdata$delta_ct_max_MminusE <= -5] <- -5
# 0 will be non infected :
qpcrdata$`delta_ct_max_MminusE+5` <- qpcrdata$delta_ct_max_MminusE + 5
# 1 will be non infected :
qpcrdata$`delta_ct_max_MminusE+6` <- qpcrdata$delta_ct_max_MminusE + 6
# presence/absence
qpcrdata$presence_eimeria_tissues <- 1
qpcrdata$presence_eimeria_tissues[qpcrdata$delta_ct_max_MminusE == -5] <- 0
qpcrdata$presence_eimeria_tissues <- as.factor(qpcrdata$presence_eimeria_tissues)
table(qpcrdata$presence_eimeria_tissues)

qpcrdata$presence_eferrisi_identified <- 0
qpcrdata$presence_eferrisi_identified[grep("ferrisi", qpcrdata$eimeriaSpecies)] <- 1
table(qpcrdata$presence_eferrisi_identified)

getinfotab(qpcrdata)

# for model intensity
qpcr_intensity_data <- qpcrdata[qpcrdata$`delta_ct_max_MminusE+5` > 0,]

cleanData$UsedForEimeriaRes[cleanData$Mouse_ID %in% qpcrdata$Mouse_ID] <- "yes"

##### All mice that were investigated for pinworms were also investigated for other helminths ##### 
pinwormsdata_bal <- cleanData[!is.na(cleanData$Aspiculuris_Syphacia),]
idToCorrect <- pinwormsdata_bal[rowSums(is.na(pinwormsdata_bal[,listWorms])) > 0, "Mouse_ID"]
cleanData[cleanData$Mouse_ID %in% idToCorrect, listWorms][
  is.na(cleanData[cleanData$Mouse_ID %in% idToCorrect, listWorms])] <- 0
pinwormsdata_bal <- cleanData[!is.na(cleanData$Aspiculuris_Syphacia),]

pinwormsdata_bal$`Aspiculuris.Syphacia+1` <-
  pinwormsdata_bal$Aspiculuris_Syphacia + 1
pinwormsdata_bal$presence_oxyurids <- 1
pinwormsdata_bal$presence_oxyurids[pinwormsdata_bal$Aspiculuris_Syphacia == 0] <- 0
pinwormsdata_bal$presence_oxyurids <-  as.factor(pinwormsdata_bal$presence_oxyurids)

getinfotab(pinwormsdata_bal)

cleanData$UsedForPinwormsRes[cleanData$Mouse_ID %in% pinwormsdata_bal$Mouse_ID] <- "yes"

##### Body condition index in Eimeria qpcr ##### 
getBodyCondition <- function(df){
  df <- df[!is.na(df$Body_length) & !is.na(df$Body_weight) & !is.na(df$Sex),]
  # Remove pregnant/post partum and juveniles
  df <- df[!df$Status %in% c("young", "pregnant"),]
  df <- df[df$Body_length > 50,]
  # Regression of BM/BS. Advantage: independant of size!!
  # Step 1: fit the model
  fitRes <- lm(Body_weight ~ Body_length * Sex, data = df)
  # Step 2: obtain predicted and residual values
  df$predicted <- predict(fitRes)   # Save the predicted values
  df$residuals <- residuals(fitRes) # Save the residual values -> to be used as indices!
  # # plot of residuals by sex
  # Plot the actual and predicted values (supplementary figure)
  myplot <- ggplot2::ggplot(df, ggplot2::aes(x = Body_length, y = Body_weight)) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, color = "lightgrey") +  # Plot regression slope
    ggplot2::geom_segment(ggplot2::aes(xend = Body_length, yend = predicted)) +
    ggplot2::geom_point(size = 4, pch = 21, alpha = .8,
                        aes(fill = HI)) +
    ggplot2::scale_fill_gradient(low = "blue", high = "red")+
    ggplot2::geom_point(ggplot2::aes(y = predicted), shape = 1) +
    ggplot2::facet_grid(~ Sex, scales = "free_x") +  # Split panels here by `iv`
    ggplot2::theme_bw()  # Add theme for cleaner look
  return(list(df, myplot))
}

body_data_eimeria <- getBodyCondition(qpcrdata)[[1]]
figResEimeria <- getBodyCondition(qpcrdata)[[2]]

cleanData$UsedForEimeriaImpactHealth[cleanData$Mouse_ID %in% body_data_eimeria$Mouse_ID] <- "yes"

##### Body condition index in pinworms ##### 
body_data_pinworms <- getBodyCondition(pinwormsdata_bal)[[1]]
figResWorm <- getBodyCondition(pinwormsdata_bal)[[2]]

cleanData$UsedForPinwormsImpactHealth[cleanData$Mouse_ID %in% body_data_pinworms$Mouse_ID] <- "yes"

# clean farms
cleanData$farm <- as.numeric(factor(cleanData$farm))

write.csv(cleanData, "../data/cleanedData.csv", row.names = F)

