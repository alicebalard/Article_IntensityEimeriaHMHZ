source("DataPreparation.R")

# A. General parasitology
myQuantitativeParasitology(pinwormsdata_watwm$Aspiculuris.Syphacia)

myQuantitativeParasitology(pinwormsdata_bal$Aspiculuris_Syphacia)

myQuantitativeParasitology(qpcrdata$`delta_ct_max_MminusE+5`)
mean(qpcrdata$`delta_ct_max_MminusE+5`[qpcrdata$`delta_ct_max_MminusE+5` > 0]) - 5
median(qpcrdata$`delta_ct_max_MminusE+5`[qpcrdata$`delta_ct_max_MminusE+5` > 0]) - 5
hist(qpcrdata$`delta_ct_max_MminusE+5`[qpcrdata$`delta_ct_max_MminusE+5` > 0])

