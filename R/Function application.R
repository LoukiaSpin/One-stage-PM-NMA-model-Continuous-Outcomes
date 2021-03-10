#######################################################################################
#                                                                                     #
#                                                                                     #
#                Applying the 'run.model' function to two NMA examples                #
#                                                                                     #
#                                                                                     #
#######################################################################################



## Load necessary library
library(dplyr); library(R2jags)



## Load the NMA examples
stowe <- read.table("./data/21370258_Stowe(2011).txt", header = T, na.strings = "NA")
schwin <- read.table("./data/24996616_Schwingshackl(2014).txt", header = T, na.strings = "NA")



## Load the 'run.model' R function
source("./R/run.model_function.R")



## Run model
# Stowe
run.model(data = stowe, measure = "MD", assumption = "HIE-TRIAL", mean.misspar = 0, var.misspar = 1, D = 0, n.chains = 2, n.iter = 10000, n.burnin = 1000, n.thin = 1)

# Schwingshackl
run.model(data = schwin, measure = "MD", assumption = "HIE-TRIAL", mean.misspar = 0, var.misspar = 1, D = 0, n.chains = 2, n.iter = 10000, n.burnin = 1000, n.thin = 1)

