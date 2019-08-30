source("DataPreparation.R")

# Plot all samples : pinworms
mapAllpinworms <- plotMap(df = pinwormsdata_bal)

# Plot infected samples : pinworms
mapInfpinworms <- plotMap(df = pinwormsdata_bal[pinwormsdata_bal$Aspiculuris_Syphacia >=1,])

# Plot all samples : eimeria (tissues)
mapAllqPCR <- plotMap(df = qpcrdata)

# Plot infected samples : eimeria (tissues)
mapInfqPCR <- plotMap(df = qpcrdata[qpcrdata$`delta_ct_max_MminusE+5` > 0,])
