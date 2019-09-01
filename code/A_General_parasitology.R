source("DataPreparation.R")

# A. General parasitology
myQuantitativeParasitology(pinwormsdata_watwm$Aspiculuris.Syphacia)

myQuantitativeParasitology(pinwormsdata_bal$Aspiculuris_Syphacia)

myQuantitativeParasitology(qpcrdata$`delta_ct_max_MminusE+5`)
mean(qpcrdata$`delta_ct_max_MminusE+5`[qpcrdata$`delta_ct_max_MminusE+5` > 0]) - 5
median(qpcrdata$`delta_ct_max_MminusE+5`[qpcrdata$`delta_ct_max_MminusE+5` > 0]) - 5
hist(qpcrdata$`delta_ct_max_MminusE+5`[qpcrdata$`delta_ct_max_MminusE+5` > 0])

## histograms of all detected nematodes (after removing BALdata embrios)

apply(BALdata[,listWorms], 2, table)

makeHistWorm <- function(i){
  Ns <- table(!is.na(BALdata[listWorms[i]]))[2]
  xpos <- max(BALdata[listWorms[i]], na.rm = T)/2
  p <- ggplot(BALdata, aes_string(x = listWorms[i])) + geom_histogram(col = "black") + theme_bw() + 
    geom_text(aes(x=xpos,y=30,label=paste0("N = ", Ns)),vjust=0) +
    scale_y_log10()
  return(p)
}

plot_grid(makeHistWorm(1), makeHistWorm(2), makeHistWorm(3),
          makeHistWorm(4), makeHistWorm(5), makeHistWorm(6),
          labels=c("A", "B", "C", "D", "E", "F"), label_size=15, ncol=2)

#### CHECK WITH JENNY --> are Nas Nas or zeros???
         
