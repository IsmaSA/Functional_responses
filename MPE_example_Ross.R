rm(list=ls())

library(frair)
library(readxl)
library(plyr)
library(MuMIn)
library(emmeans)
library(car)
library(bbmle)
library(emdbook)
library(dataframes2xls)
library(AICcmodavg) 
library(ggplot2)
library(reshape2)
library(deSolve) 
library(FME)

params<-read.csv("params.csv")
params<-read.csv2("Ross_claws.csv")
# you should alter this file with the parameters for each claw treatment
# you could just run to the code three times, with all pairwise combinations of the three groups
# then you have to change less below (but just repeat the scripts over and over)

params$Std..Error <-as.numeric(params$Std..Error)
params$CI<-(params$Std..Error*1.96)

n<-params
n$Estimate <- as.numeric(n$Estimate)
n$CI<- as.numeric(n$CI)

Pars16 <- c(nilea = n[3,1], 
            nileh = n[4,1],
            omoa = n[5,1],
            omoh = n[6,1]) 

Pars16


parRanges16_IC <- data.frame(min = c(n[3,1]-n[3,7], n[4,1]-n[4,7], 
                                     n[5,1]-n[5,7], n[6,1]-n[6,7]), 
                             max = c(n[3,1]+n[3,7], n[4,1]+n[4,7], 
                                     n[5,1]+n[5,7], n[6,1]+n[6,7]))
rownames(parRanges16_IC)<- c("nilea","nileh","omoa","omoh")
parRanges16_IC

parCovar16 =diag(x = 1, nrow=4, ncol=4)

# The below codes just fit the model and determine CIs

MPE <- function (Pars,Nnile,Nomo,prey_density){
  model <- function(Time, State, Pars){ 
    with (as.list(c(State,Pars)), { 
      dNd <- -nilea*Nd*Nnile / (1 + nilea*nileh*Nd) - omoa*Nd*Nomo / (1 + omoa*omoh*Nd)  
      return(list(c(dNd))) 
    }) 
  }
  
  Time <- seq(0, 1, by = 0.01) ## e.g. 6 hours. must fit to your units of the handling times and attack rates (e.g. hours, seconds) 
  #Time<-24
  State <- c(Nd = prey_density) ## just some number, must be replaced with your experimental starting densities 
  
  results <- as.data.frame(ode(func = model, y = State, parms = Pars,method = rkMethod("rk45dp7"), times = Time))
  results$eaten<-(prey_density-results$Nd)
  return(results)
}

summary.sensRange<-function(object, ...) {
  
  npar <- attr(object, "npar")
  sens <- as.matrix(object[, -(1:npar)])
  x    <- attr(object, "x")
  names(x) <- NULL
  nx   <- attr(object, "nx")
  varnames  <- attr(object, "var")
  
  if (ncol(sens)>1)
    SumSens <- data.frame(
      x    = x,
      Mean = apply(sens, 2, FUN = mean),
      Sd   = apply(sens, 2, FUN = sd),
      Min  = apply(sens, 2, FUN = min),
      Max  = apply(sens, 2, FUN = max),
      q025  = apply(sens, 2, FUN=function(x)quantile(x, probs = 0.025)),
      q05  = apply(sens, 2, FUN=function(x)quantile(x, probs = 0.05)),
      q25  = apply(sens, 2, FUN=function(x)quantile(x, probs = 0.25)),
      q50  = apply(sens, 2, FUN=function(x)quantile(x, probs = 0.5)),
      q75  = apply(sens, 2, FUN=function(x)quantile(x, probs = 0.75)),
      q95  = apply(sens, 2, FUN=function(x)quantile(x, probs = 0.95)),
      q975  = apply(sens, 2, FUN=function(x)quantile(x, probs = 0.975))
    ) else
      SumSens <- data.frame(
        x    = x,
        Mean = mean(sens, na.rm = TRUE),
        Sd   = NA,
        Min  = min(sens, na.rm = TRUE),
        Max  = max(sens, na.rm = TRUE),
        q025  = quantile(sens, probs = 0.025, na.rm = TRUE),
        q05  = quantile(sens, probs = 0.05, na.rm = TRUE),
        q25  = quantile(sens, probs = 0.25, na.rm = TRUE),
        q50  = quantile(sens, probs = 0.5, na.rm = TRUE),
        q75  = quantile(sens, probs = 0.75, na.rm = TRUE),
        q95  = quantile(sens, probs = 0.95, na.rm = TRUE),
        q975  = quantile(sens, probs = 0.975, na.rm = TRUE)
      )
  
  rownames(SumSens) <- colnames(sens)
  attr(SumSens, "var") <- attr(object, "var")
  attr(SumSens, "nx")  <- attr(object, "nx")
  class(SumSens)<-c("summary.sensRange", "data.frame")
  
  return(SumSens)
}

conf_int<-function(Pars,parRanges,parCovar,Nnile,Nomo,prey_density){
  Pars2=Pars
  parRanges2=parRanges
  parCovar2=parCovar
  ##assume a normal distribution but I do not have the full variance coriance matrix for the parameters
  #SensR <- sensRange(func = MPE, parms = Pars, NC=Nb_C,NL=Nb_L,NS=Nb_S,prey_density= Daphnia, dist = "norm",
  #                   sensvar = "IS", parMean = Pars2, parCovar =parCovar2,  
  #                 parRange = parRanges2, num = 100)
  
  ##latin hypercube sampling algorithm
  SensR <- sensRange(func = MPE, parms = Pars, Nnile=Nnile, Nomo=Nomo, prey_density=prey_density, 
                     dist = "latin",
                     sensvar = "eaten",   
                     parRange = parRanges2, num = 100)
  Sens <- summary.sensRange(SensR, na.rm= TRUE)
  colnames(Sens)[1]<- c("Time")
  colnames(Sens)[2]<- c("eaten_predicted")
  res<-Sens[101,]
  
  res$Prey_density <-prey_density
  res$nile <- Nnile
  res$omo <- Nomo
  
  return(res)
}

# Then fit the model below to get predictions for your prey densities with each inter/intraspecific predator combination
# You adjust the numbers on the right for "no. pred. treatment1, no pred. treatment2, prey density"
# e.g., "1,1,12" is 1 x predator type 1, 0 x predator type 2, and 12 prey", so it calculates the predicted
# number of prey eaten for those treatments; "2,0,24" is "2 x predator type 1, 0 x predator type 2, and 24 prey"

# 0 vs 0

k1<-conf_int(Pars16,parRanges16_IC,parCovar16, 2,0,12)
k2<-conf_int(Pars16,parRanges16_IC,parCovar16, 2,0,24)
k3<-conf_int(Pars16,parRanges16_IC,parCovar16, 2,0,36)

# 1 vs 1
k4<-conf_int(Pars16,parRanges16_IC,parCovar16, 0,2,12)
k5<-conf_int(Pars16,parRanges16_IC,parCovar16, 0,2,24)
k6<-conf_int(Pars16,parRanges16_IC,parCovar16, 0,2,36)

# 0 vs 1
k7<-conf_int(Pars16,parRanges16_IC,parCovar16, 1,1,12)
k8<-conf_int(Pars16,parRanges16_IC,parCovar16, 1,1,24)
k9<-conf_int(Pars16,parRanges16_IC,parCovar16, 1,1,36)


# 0 vs 2
k10<-conf_int(Pars16,parRanges16_IC,parCovar16, 1,1,12)
k11<-conf_int(Pars16,parRanges16_IC,parCovar16, 1,1,24)
k12<-conf_int(Pars16,parRanges16_IC,parCovar16, 1,1,36)


# 1 vs 2
k13<-conf_int(Pars16,parRanges16_IC,parCovar16, 1,1,12)
k14<-conf_int(Pars16,parRanges16_IC,parCovar16, 1,1,24)
k15<-conf_int(Pars16,parRanges16_IC,parCovar16, 1,1,36)


# 2 vs 2
k16<-conf_int(Pars16,parRanges16_IC,parCovar16, 0,2,12)
k17<-conf_int(Pars16,parRanges16_IC,parCovar16, 0,2,24)
k18<-conf_int(Pars16,parRanges16_IC,parCovar16, 0,2,36)


results<-rbind(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,k18)
write_xlsx(results, "Non-trophic_result2.xlsx")



