######################################
## Run twin model in JAGS
######################################

#Install and load packages needed:
install.packages("rjags");library(rjags)

setwd("C:/Users/schwabei/Dropbox/MCMC/Part 2 JAGS") #set working directory

#Simulate twin data to be able to run JAGS analysis: 
install.packages("devtools")
library(devtools)
install_github("ingaschwabe/BayesTwin")
library(BayesTwin)
?simulate_twin_data

#We simulate N = 30 item  data for 500 MZ twin pairs and 700 DZ twin pairs. 
#varA is 0.5, varC is 0.3 and varE (residual variance) is 0.2. No covariates
twin_data = simulate_twin_data(nmz = 500, ndz = 700, var_a = 0.5, 
                               var_c = 0.3, var_e = 0.2, model = "ACE", 
                               n_items = 30, n_var = 0)


#col. 1-30 = item data from twin 1, col. 31-60 = item data from twin 2. 
head(twin_data)
names(twin_data)
twin_data$y_mz #MZ twin pairs data
twin_data$y_dz #DZ twin pairs data

#==========================================================
# Analyse data with JAGS 
#==========================================================
jagsdata = list("n_mz" = 500, "n_dz" = 700, "n_items" = 30,
                "Ymz" = twin_data$y_mz,
                "Ydz" = twin_data$y_dz)

#Het model staat in jags.txt, we gebruiken geen eigen start waarden (inits = NULL) en 1 markov chain:
jags <- jags.model("twin_model_jags.txt", jagsdata, inits = NULL, n.chains = 1, quiet=FALSE)
update(jags, 12000)
out <- jags.samples(jags, c("b", "tau_c", "tau_e", "tau_a"), 17000)  

#Note that we now estimate beta parameters (item parameters)
#although we know the true values because we simulated the item data. 
#we could also have used the simulated beta parameters as input for JAGS and set a prior
#distribution on mu instead of fixing mu to zero in order to identify the scale. 

## Convergence?
plot(1/out$tau_a[,,1], type = "l")
plot(1/out$tau_c[,,1], type = "l")
plot(1/out$tau_e[,,1], type = "l")
plot(out$b[1,,1], type = "l") #item parameter beta, item 1 
plot(out$b[12,,1], type = "l") #item parameter beta, item 12

## Inference
1/out$tau_a; sd(1/out$tau_a)
1/out$tau_c; sd(1/out$tau_c)
1/out$tau_e; sd(1/out$tau_e)

HPD(1/out$tau_a, 0.95)
HPD(1/out$tau_c, 0.95)
HPD(1/out$tau_e, 0.95)
