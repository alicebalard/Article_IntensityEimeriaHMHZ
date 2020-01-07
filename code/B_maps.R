source("DataPreparation.R")

# after update ggmap, change area to keep stamen-lite type
area <- get_stamenmap(bbox = c(12, 51.5, 15, 53.5), maptype = "toner-lite", zoom = 8, format = "svg")

ggsave(filename="../figures/figure 2 and 3 better quality maps/test.svg", device = "svg",
       plot=ggmap(area), width = 15, height = 15, dpi = 1000)


pdf(file = "../figures/figure 2 and 3 better quality maps/test.pdf", height = 5, width = 5, pointsize = 20)
ggmap(area)
    dev.off()




# install.packages("svglite")
library(svglite)

ggsave(filename="../figures/figure 2 and 3 better quality maps/test.svg", device = "svg",
       plot=ggmap(area), width = 15, height = 15, dpi = 1000)



ggsave(filename="../figures/figure 2 and 3 better quality maps/test.svg", device = "svg",
       plot=ggmap(area),width=10,height=8,units="cm")


# Plot all samples : pinworms
mapAllpinworms <- plotMap(df = pinwormsdata_bal)

# Plot infected samples : pinworms
mapInfpinworms <- plotMap(df = pinwormsdata_bal[pinwormsdata_bal$Aspiculuris_Syphacia >=1,])

# Plot all samples : eimeria (tissues)
mapAllqPCR <- plotMap(df = qpcrdata)

# Plot infected samples : eimeria (tissues)
mapInfqPCR <- plotMap(df = qpcrdata[qpcrdata$`delta_ct_max_MminusE+5` > 0,])

## maps in better quality for fig 2 and 3

pdf(file = "../figures/figure 2 and 3 better quality maps/fig2A.pdf", height = 30, width = 30)
mapAllpinworms
dev.off()

# install.packages("mapdata")
# library(mapdata)
# 
# 
# install.packages("rmapzen")
# library(rmapzen)
# 
# 
# nextzen_API_key="../../keynextzen"
options(nextzen_API_key="../../keynextzen")
mz_set_tile_host_nextzen(key = getOption("nextzen_API_key"))

#install.packages("cancensus")
library(cancensus)
victoria <- get_census(dataset='CA16',
                        regions=list(CSD=c("5917034","5917040","5917021","5917030","5917041","5917047")),
                        vectors=c(), labels="detailed", geo_format='sf', level='CT')
# 
# library(ggplot2)
# ggplot(victoria)
