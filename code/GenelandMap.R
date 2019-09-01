source("DataPreparation.R") 

# Install Geneland modified by Alice
# install_github("alicebalard/Geneland")
library("Geneland")

# Load data 
mydata <- read.csv("../data/cleanedData.csv", na.strings=c(""," ","NA"))
mydata <- mydata[mydata$UsedForMap %in% "yes",]
    nrow(mydata) # N = 598 mice

# Create a folder "outputGeneland" in the folder where you have stored this R script
mainDir <- "../data"
subDir <- "outputGeneland"
dir.create(file.path(mainDir, subDir))

# Choose the markers to work with as well as Latitude and Longitude columns
# Here, 6 diploid markers
mydata <- mydata[names(mydata) %in% c("Es1C", "Gpd1C", "Idh1C", 
                                      "MpiC", "NpC", "Sod1C", 
                                      "Latitude", "Longitude")]
mydata <- na.omit(mydata)
nrow(mydata) # N = 598 mice, check that it didn't change

# Save the genotypes 
geno <- mydata[names(mydata) %in% c("Es1C", "Gpd1C", "Idh1C", "MpiC", "NpC", "Sod1C")]

# Visualise the type of genotype per markers
sapply(geno, levels)

# Format the genotypes table
geno <- as.data.frame(
  sapply(geno, function(x) revalue(x, c("dd" = "d/d", 
                                        "mm" = "m/m", 
                                        "dm" = "d/m"))))
geno[geno == ""] <- NA

geno <- lapply(geno, function(x) data.frame(do.call('rbind', strsplit(as.character(x), '/', fixed=TRUE))))
geno <- do.call(cbind, geno)

# Geneland needs integers instead of m and d:
geno <- data.frame(sapply(geno, as.numeric))

# Coordinates
coord <- data.frame(Longitude = mydata$Longitude, Latitude = mydata$Latitude)

## 2. Inference
plot(coord,xlab="Eastings",ylab="Northings",asp=1)

Geneland::MCMC(coordinates= coord,
               geno.dip.codom = geno,
               path.mcmc = "../data/outputGeneland/",
               nit=1000000,
               npopmax = 2,
               varnpop = FALSE,
               npopinit = 2,
               spatial=TRUE,
               freq.model="Correlated",
               thinning = 100)

# Post-processing MCMC outputs
Geneland::PostProcessChain(coordinates=coord,
                           path.mcmc="../data/outputGeneland/", 
                           nxdom = 1400,
                           nydom=100,
                           burnin=200)

## Generating graphical outputs "GENELAND results for K=n using the spatial model with correlated allele frequencies"
Geneland::Plotnpop(path.mcmc="../data/outputGeneland/",
                   burnin=200)

Geneland::PosteriorMode(coordinates=coord,
                        path.mcmc="../data/outputGeneland/")

## New plot w colors:
Geneland::PlotTessellation(coordinates = coord,
                                     path.mcmc="../data/outputGeneland/",
                                     colfancy =  colorRampPalette(c("blue", "red"))(100), 
                                     plotCities = T, path = "./figures/")
