##############################################
## Run JAGS and OpenBugs from R 
## Using simple regression analysis with one 
## covariate as example
##############################################

#Install and load packages needed:
install.packages("rjags"); install.packages("coda"); install.packages("R2OpenBUGS")
install.packages("devtools")
install.packages('R2OpenBUGS',type='source') #for mac users
library(R2OpenBUGS)
library(rjags)
library(coda)
library(devtools)

#Installeren van libraries hoeft maar 1 keer, laden moet elke keer als je R start weer opnieuw

#Set working directory: 
setwd("/Users/inga/Dropbox/MCMC/2 JAGS and OpenBugs") #hier het pad v.d. bestanden aangeven.

#Simulate some data for the regression analysis: 
set.seed(1030)
N = 100 
mu = rnorm(1, 0, 10) 
beta_coef = runif(1, -3, 3) 

#Simulate data:
x <- rnorm(100, 0, 10)
y <- mu + (beta_coef * x) + rnorm(N, 0, .1)

#==========================================================
# Analyse data with JAGS 
#==========================================================
#1: We maken een "list" aan waarin wij de variabelen/data die we nodig hebben opslaan:
jagsdata = list("N" = N, #"N" = naam, N = objekt/getal 
                "x" = x,
                "y" = y)

#Het model staat in jags.txt, we gebruiken geen eigen start waarden (inits = NULL) en 1 markov chain:
jags <- jags.model("jags.txt", jagsdata, inits = NULL, n.chains = 1, quiet=FALSE)

#10000 (in dit geval een beetje overdreven) samples als burn-in:
update(jags, 10000)

#20000 samples (na de burn in), en we willen de schattingen van de 
#parameters b, mu en tau_r: 
out <- jags.samples(jags, c("b", "mu", "tau_r"), 20000)  

## Convergence (zelf plotten)
plot(out$mu, type = "l", xlab = "Iteration", ylab = "Estimate", 
     main = expression(mu))
plot(out$b, type = "l")
plot(1/out$tau_r, type = "l") #1/tau_r omdat tau_r precisie en niet variantie is!

## Convergence with coda library: 
#hiervoor moeten we de analyse nog een keer runnen (samples nog een keer trekken) met 
#de functie coda.samples ipv jags.samples: 
out_coda_jags = coda.samples(jags, c("b", "mu", "tau_r"), 20000) 
xyplot(out_coda_jags)

#Gelman & Ruben statistic with two chains
jags <- jags.model("jags.txt", jagsdata, inits = NULL, n.chains = 2, quiet=FALSE)
update(jags, 10000)
out_2chains <- coda.samples(jags, c("b", "mu", "tau_r"), 20000)  
gelman.diag(out_2chains) #upper CI must be < 1

## Results: 
out$mu #punt schattingen
out$b
1/out$tau_r

plot(hist(out$mu))
plot(hist(out$b))
plot(density(out$mu, adjust = 5)) #density posterior distribution
plot(density(out$b, adjust = 5)) 
densityplot(out_coda_jags) #density plot met coda package

### HPD: 
HPD <- function(sample1, rel.int) {
    rel.int <- (1 - rel.int)/2 #calculate range outside of credibility region (both sides 2.5 in this case)
    lower <- round(length(sample1) * rel.int, 0) 
    upper <- round(length(sample1) * (1 - rel.int), 0)
    diff.int <- upper - lower
    HPDo <- sample1[order(sample1)][1:lower]
    HPDb <- sample1[order(sample1)][(diff.int + 1):upper]
    HPDI <- round(c(HPDo[order(HPDb - HPDo)[1]], HPDb[order(HPDb - HPDo)[
        1]]), 5)
    #CI <- round(c(sample1[order(sample1)][lower], sample1[order(sample1)][
    #  upper]), 3)
    return(HPDI)
}

HPD(out$mu, 0.95)
HPD(out$b, 0.95)
HPD(1/out$tau_r, 0.95)


#==========================================================
# Analyse data with OpenBugs 
#==========================================================

#for the mac user follow these steps: 
#http://www.davideagle.org/r-2/bayesian-modeling-using-winbugs-and-openbugs/running-openbugs-on-mac-using-wine

bugsdata = list("N" = N, 
                "x" = x,
                "y" = y)

bugs_out <- bugs(bugsdata, inits = NULL, model.file = "bugs.txt",
                 parameters = c("mu", "b", "tau_r"),
                 n.chains = 1, n.burnin = 10000, n.iter = 20000)

#Something went wrong? use: debug = TRUE
bugs_out <- bugs(bugsdata, inits = NULL, model.file = "bugs.txt",
                 parameters = c("mu", "b", "tau_r"),
                 n.chains = 1, n.burnin = 10000, n.iter = 20000,
                 debug = TRUE)

#Convergence: 
#To use CODA for analyzing the MCMC convergence, we have to enable the c
#odaPkg option, which allows us to convert the output for CODA using read.bugs.

bugs_out <- bugs(bugsdata, inits = NULL, model.file = "bugs.txt",
                 parameters = c("mu", "b", "tau_r"),
                 n.chains = 1, n.burnin = 10000, n.iter = 20000,
                 codaPkg = TRUE)

out_coda = read.bugs(bugs_out)
xyplot(out_coda)

#Convergence: gelman rubin statistic: 
bugs_out <- bugs(bugsdata, inits = NULL, model.file = "bugs.txt",
                 parameters = c("mu", "b", "tau_r"),
                 n.chains = 2, n.burnin = 10000, n.iter = 20000,
                 codaPkg = TRUE)
out_coda = read.bugs(bugs_out)
gelman.diag(out_coda) 

#Inference: 
densityplot(out_coda)
bugs_out
