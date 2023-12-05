# Name: Aryan Verma    UUN: s2512060

library(rjags) # Load jags and coda libraries
library(coda)
library(ggplot2) # ggplot2 for graphs and plots

#getwd()
#setwd("/Users/aryanverma/Desktop/ExtendedSP/Assignment 5/")

## Function for developing the effect dose graph
#  Input:
#  df <- Data from vin.txt file
#  samples <- gsamples generated from model
#  model <- the model we want to generate plot for
#  model_file <- integer value indicating which model- 1 or 2 
## Output: The Effect-Dose plot with model and 95% Credible Interval
effectDoseGraph <- function(df, samples, model, model_file){
  
  ## Calculating the posterior expected effect from model
  temp = do.call(rbind.data.frame, samples) ## Convert to dataframe
  lower<-rep(NA, 45) ##  lower bound of 95% CI
  mean<-rep(NA, 45) ## mean expected effect from model
  upper<-rep(NA, 45) ## upper bound of 95% CI
  for(j in 1:5){ ## Loop through experimental runs
    for(i in 1:9){ ## Loop through concentration values
      M = gsub(" ", "", paste("M[",j,"]")) ## M index
      if(model_file==1){
        mu = gsub(" ", "", paste("mu[",i,"]")) ## Model1 - use mu_i
      } else{
        mu = gsub(" ", "", paste("mu[",i,",",j,"]")) ## Model2 - use mu_i_j
      }
      predicted_effect = temp[M]*temp[mu] ## seperate column Mj * mu_i
      mean[(j-1)*9+i]<-mean(predicted_effect[[1]]) ## mean expected effect
      lower[(j-1)*9+i]<-quantile(predicted_effect[[1]], 0.025)[["2.5%"]] ##lower bound
      upper[(j-1)*9+i]<-quantile(predicted_effect[[1]], 0.975)[["97.5%"]] ##upper bound
    }
  }
  ## Create a data frame for plotting the effect against dose
  plot_data <- data.frame(
    Dose = rep(1:9, times = 5),  ## Assuming 5 experiments for each dose level
    Effect = as.vector(df$effect), ## The original effect from data
    Exp_Run = rep(1:5, each = 9), ## Experimental run
    expected_model_effect <- mean, ## Mean Expected effect from model
    min<-lower, ## lower bound for 95% CI
    max<-upper ## Upper bound for 95% CI
  )
  
  # Plotting the data
  p<-ggplot(plot_data, aes(x = Dose, y = Effect, color = factor(Exp_Run))) + ## data points
    geom_point() +
    geom_line(aes(x=Dose, y=expected_model_effect)) + ## Model as lines
    geom_ribbon(aes(x=Dose, ymin = min, ymax = max), alpha = 0.2)+  ## 95% credible intervals
    theme_minimal() +
    labs(title = paste("Effect vs. Dose: Model",model_file), ## Title of graph
         x = "Dose Index", ## x-axis of graph
         y = "Effect")  ## y-axis of graph
  plot(p) ## Plot the above graph
}

data_path <- "vin.txt" ## dataset file path
model_file1 <- "M1s2512060.jags" ## Model1 file path
model_file2 <- "M2s2512060.jags" ## Model2 file path

df = read.table(data_path, header=TRUE) ## Reading the data

data_matrix <- 
  with(df, ##dataframe to a matrix
          matrix(effect, 
                nrow = length(unique(conc)), ## Rows for Concentration
                ncol = length(unique(exper)) ## Cols for experiment
                ))

data <- list( ## Converting matrix to a list
  n_concentrations = 9, ## Concentration in rows
  n_experiments = 5, ## Experimental runs in columns
  Y = data_matrix ## Values as effect observed
)
model1 <- jags.model(model_file1, data = data) #Model Compiled

update(model1, 1000) #Burning in
samples <- coda.samples(model1, ## generating samples
                        c("x", "M","tau", "tau_o"), # Parameters
                        n.iter = 27000) # iterations = 27000

## Calculate effective sample size from generated samples
effective_sample_sizes <- effectiveSize(samples)
highest_ess <- names(which.max(effective_sample_sizes)) ##Highest ESS
lowest_ess <- names(which.min(effective_sample_sizes)) ##Lowest ESS

## Generate Plots for highest and lowest ESS
plot(samples[,highest_ess],main=paste("Highest ESS:",highest_ess))
plot(samples[,lowest_ess],main=paste("Lowest ESS:",lowest_ess))


## Creating Effect vs Dose graph for Model1
samples1 <- coda.samples(model1, ## generating samples with required variables
                        c("M","mu"), # Parameters
                        n.iter = 27000) # iterations = 27000
effectDoseGraph(df, samples1, model1, 1) ## Generate graph

## Creating effect vs Dose graph for Modified Model as Model 2
model2 <- jags.model(model_file2, data = data) #Model Compiled

update(model2, 1000) #Burning in
samples2 <- coda.samples(model2, ## generating samples
                        c("M", "mu"), # Parameters
                        n.iter = 27000) # iterations = 27000
effectDoseGraph(df, samples2, model2, 2)

# Comparing the models using DIC
model1 <- jags.model(model_file1, data = data, n.chains = 3) ## compile model1 with 3 chain
model2 <- jags.model(model_file2, data = data, n.chains = 3) ## compile model2 with 3 chain

# Burn-in and sample for model 1
update(model1, 1000)
samples1 <- coda.samples(model1, c("x", "M", "tau", "tau_o"), n.iter = 27000)

# Burn-in and sample for model 2
update(model2, 1000)
samples2 <- coda.samples(model2, c("x", "M", "tau", "tau_o"), n.iter = 27000)

dic1 <- dic.samples(model1, n.iter = 20000) ## Model1 DIC
dic2 <- dic.samples(model2, n.iter = 20000) ## Model2 DIC

#print(dic1) ##Printing DIC values
#print(dic2)
## As the penalized deviance is lower for the model 2, it is preferable
