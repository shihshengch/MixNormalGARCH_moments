# 0_data_download.R

# clear
rm(list = ls())

# library
if (!require("quantmod")) install.packages("quantmod")
library(quantmod)

# get data
getSymbols(Symbols = "^GSPC", from = "1990-01-01", to = "2017-02-28")
GSPC = na.omit(GSPC$GSPC.Close)

# compute return
r = na.omit(diff(log(GSPC))) * 100

# compute epsilon
c = mean(r)
eps = r - c

# keep necessary data
rm(list = setdiff(ls(), c("eps", "r", "c")))

# save data
save.image("data.RData")
