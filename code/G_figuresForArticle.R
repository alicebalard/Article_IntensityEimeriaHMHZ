source("DataPreparation.R")
source("B_maps.R")
source("D_Models_parasiteLoads.R")
source("I_prevalenceAlongHI.R")
source("F_ModelBodyCondition.R")

# general map
geogDF <- rbind(data.frame(latmin = min(WATWMdata$Latitude), latmax = max(WATWMdata$Latitude),
                           lonmin = min(WATWMdata$Longtitude), lonmax = max(WATWMdata$Longtitude), transect = "watwm"),
                data.frame(latmin = min(BALdata$Latitude), latmax = max(BALdata$Latitude),
                           lonmin = min(BALdata$Longitude), lonmax = max(BALdata$Longitude), transect = "bal"),
                data.frame(latmin = min(BALdata$Latitude, WATWMdata$Latitude),
                           latmax = max(BALdata$Latitude, WATWMdata$Latitude),
                           lonmin = min(BALdata$Longitude, WATWMdata$Longtitude),
                           lonmax = max(BALdata$Longitude, WATWMdata$Longtitude), transect = "all"))

areaBothStudies <- get_map(location = c(0, 38, 30, 60),
                           source = "stamen", 
                           maptype = "toner-lite", 
                           zoom = 5)

# Create HI bar
HIgradientBar <- ggplot(data.frame(hi = seq(0,1,0.0001)),
                        aes(x=hi, y=1, fill = hi)) +
  geom_tile() +
  theme_void() +
  scale_fill_gradient(low = "blue", high = "red")  + 
  scale_x_continuous(expand=c(.01,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  theme(legend.position = 'none')

# fig 2
figEimeria <- plot_grid(
  annotate_figure(plot_grid(mapAllqPCR, 
                            prevalenceAlongHIEimeria +
                              theme(legend.position = "none"),
                            HIgradientBar,  
                            nrow = 3, 
                            rel_heights = c(1.5, 1.3, 1/8))),
  annotate_figure(plot_grid(mapInfqPCR, 
                            plotInfEimeria,
                            HIgradientBar,  
                            nrow = 3, 
                            rel_heights = c(1.5, 1.3, 1/8)))
)
pdf(file = "../figures/figEimeria.pdf",
    width = 12, height = 10)
figEimeria
dev.off()

# fig 3
figPinworms <- plot_grid(
  annotate_figure(plot_grid(mapAllpinworms, 
                            prevalenceAlongHIpinworms +
                              theme(legend.position = "none"),
                            HIgradientBar,  
                            nrow = 3, 
                            rel_heights = c(1.5, 1.3, 1/8))),
  annotate_figure(plot_grid(mapInfpinworms, 
                            plotInfWorms +
                              theme_bw() +
                              theme(legend.position = "none"),
                            HIgradientBar,  
                            nrow = 3, 
                            rel_heights = c(1.5, 1.3, 1/8)))
)

pdf(file = "../figures/figPinworms.pdf",
    width = 12, height = 10)
figPinworms
dev.off()

# Body residuals 

# eimeria
pdf(file = "../figures/residualsEimeria.pdf",
    width = 10, height = 10)
plot_grid(figResWorm, 
          plotResiduals_pinworms,
          nrow = 2)
dev.off()

# pinworms
pdf(file = "../figures/residualsWorms.pdf",
    width = 10, height = 10)
plot_grid(figResEimeria, 
          plotResiduals_eimeria,
          nrow = 2)
dev.off()

