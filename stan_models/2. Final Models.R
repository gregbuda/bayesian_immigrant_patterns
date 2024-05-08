
##############################
#0. Packages
##############################

#install.packages("rstan")
#install.packages("rstanarm")
#install.packages("tidyverse")
#install.packages("rstan")

library(rstan)
library(ggplot2)
library(bayesplot)
library(tidyverse)
library(readr)
library(loo)
theme_set(theme_bw())
options(mc.cores = parallel::detectCores())
setwd("C:/Users/gregb/Desktop/TFM")



#############################
#1. Prepare data
#############################

df <- read_csv("./intermediate_files/feature_matrix_final.csv")

#Target variable
df['pop_share'] <- df['population'] / df['census_total_pop']

#Unique CUSECS
j <- length(unique(df$cusec) )

#Scope of modeling is 5 countries: Italy, Morroco, China, Ecuador, Perú
df <- subset(df, cob %in% c('Italia','Marruecos','China','Ecuador','Perú'))

#Readjust cob_idx
df$cob_idx <- ifelse(df$cob == 'China',1,
                     ifelse(df$cob == 'Ecuador',2,
                            ifelse(df$cob == 'Italia',3,
                                   ifelse(df$cob == 'Marruecos',4,
                                          ifelse(df$cob == 'Perú',5,0)))))
i <- length(unique(df$cob) )

#BCN Overall baseline by nationality
     # Immigrant populations
df_total_im <- df %>%
  group_by(cob) %>%
  summarize(total_im = sum(population))
     # Total Population
df_total_pop <- df %>%
  group_by(cusec) %>%
  summarize(total_pop = max(census_total_pop)) %>%
  summarize(total_pop = sum(total_pop))
    #Overall density by nationality
result <- cross_join(df_total_im, df_total_pop)
result <- result %>%
  mutate(overall_density = total_im / total_pop)
    #Add to table
df <- left_join(df, result %>% select(cob, overall_density), by = "cob")


# Adjacency matrix
adj <- unlist(as.list(read.table("./intermediate_files/icar_adj.txt", header = FALSE)$V1))
# Neighborhood counts
num <- unlist(as.list(read.table("./intermediate_files/icar_num.txt", header = FALSE)$V1))
#Create Adjacency matrix
mungeCARdata4stan = function(adjBUGS,numBUGS) {
  N = length(numBUGS);
  nn = numBUGS;
  N_edges = length(adjBUGS) / 2;
  node1 = vector(mode="numeric", length=N_edges);
  node2 = vector(mode="numeric", length=N_edges);
  iAdj = 0;
  iEdge = 0;
  for (i in 1:N) {
    for (j in 1:nn[i]) {
      iAdj = iAdj + 1;
      if (i < adjBUGS[iAdj]) {
        iEdge = iEdge + 1;
        node1[iEdge] = i;
        node2[iEdge] = adjBUGS[iAdj];
      }
    }
  }
  return (list("N"=N,"N_edges"=N_edges,"node1"=node1,"node2"=node2));
}
nbs = mungeCARdata4stan(adj, num);
node1 = nbs$node1;
node2 = nbs$node2;
N_edges = nbs$N_edges;





##################################
#2. Intercept and Random Effects model
##################################

# for (country in c('Italia','Marruecos','China','Ecuador','Perú')) {
# 
#   X <- subset(df, cob==country)
#   X <- X[, !names(X) %in% c("cusec","cusec_idx", "cob","cob_idx","population","pop_share","overall_density")]
#   #X <- scale(X)
#   y <- subset(df, cob==country)$population
# 
#   #Model
#   data_list <- list(
#     #Dimensions
#     N = dim(X)[1],
#     j = j,
#     p=ncol(X),
#     #Data
#     X = X,
#     y = y,
#     #CAR
#     node1=node1,
#     node2=node2,
#     N_edges=N_edges,
#     #Indices
#     cusecs_indx = as.integer(subset(df, cob==country)$cusec_idx)
#   )
# 
#   model <- stan("poisson_intercept.stan",  data = data_list, chains=4, thin=2,
#                 iter=4000,  seed=1)
#   save(model, file = paste0("./model_intercept/model_intercept_", country, ".Rdata"))
# }




#############################
#3. Feature Selection Model
#############################

#Scaling the features is necessary here

# for (country in c('Italia','Marruecos','China','Ecuador','Perú')) {
#   
#   X <- subset(df, cob==country)
#   X <- X[, !names(X) %in% c("cusec","cusec_idx", "cob","cob_idx","population","pop_share","overall_density")]
#   X <- scale(X)
#   y <- subset(df, cob==country)$population
# 
#   Modeling
#   tau_values <- c(0.001, 0.1, 10)
#   for (tau in tau_values) {
# 
#         #Iterate over 3 possible values of tau
#     data_list <- list(
#       #Dimensions
#       N = dim(X)[1],
#       j = j,
#       p=ncol(X),
#       #Data
#       X = X,
#       y = y,
#       #CAR
#       node1=node1,
#       node2=node2,
#       N_edges=N_edges,
#       #Indices
#       cusecs_indx = as.integer(subset(df, cob==country)$cusec_idx),
#       #Penalty
#       tau=tau
#     )
# 
#     model <-  stan("poisson_horseshoe.stan", data = data_list, chains=4, thin=2,
#                    iter=4000,  seed=1, init = 0)
#     save(model, file = paste0("./model_horseshoe/model_horseshoe_", country,"_tau",tau, ".Rdata"))
# 
#   }}


#Initiate a binary matrix where we'll indicate in each setup which are the significant variables
binary_matrix <- matrix(NA, nrow = 3 * 5, ncol = 25)
counter <- 1

for (country in c('Italia','Marruecos','China','Ecuador','Peru')) {
  tau_values <- c(0.001, 0.1, 10)
  for (tau in tau_values) {

    load(paste0("./stan_final_models/model_horseshoe/model_horseshoe_", country,"_tau",tau, ".Rdata"))
    
    summary_stats <- summary(model)
    posterior_samples <- rstan::extract(model)
    num_betas <- 1 + dim(posterior_samples$betas)[2]
    
    binary_matrix[counter,1] <- country
    binary_matrix[counter,2] <- tau
    
    #Check if intercept is significant
    quantiles <- quantile(posterior_samples$beta0, c(0.025, 0.975))
    binary_matrix[counter,3] <- !(quantiles[1] <= 0 && quantiles[2] >= 0)
    
    #Check if other betas are significant
    for (j in 1:22) {
      quantiles <- quantile(posterior_samples$betas[,j], c(0.025, 0.975))
      binary_matrix[counter,j+3] <- !(quantiles[1] <= 0 && quantiles[2] >= 0)
    }
    
    counter <- counter + 1
  }}

binary_matrix<-as.data.frame(binary_matrix)
colnames(binary_matrix) <- c(c('country','tau','beta0'),colnames(df[, !names(df) %in% c("cusec","cusec_idx", "cob","cob_idx","population","pop_share","overall_density")]))

#How many significant variables are thewre per setup?
rowSums(as.matrix(binary_matrix[, 3:25]) == "TRUE")





#############################
#4. Model with CAR
#############################
#Example: https://mc-stan.org/users/documentation/case-studies/icar_stan.html

# for (country in c('Italia','Marruecos','China','Ecuador','Perú')) {
#   X <- subset(df, cob==country)
#   X <- X[, !names(X) %in% c("cusec","cusec_idx", "cob","cob_idx","population","pop_share","overall_density")]
#   X <- scale(X)
#   y <- subset(df, cob==country)$population
# 
#   tau_values <- c(0.001, 0.1, 10)
#   for (tau in tau_values) {
# 
#     #Relevant setup
#     setup_vector<-binary_matrix[binary_matrix$country == country & binary_matrix$tau == tau,]
#     #beta0 is significant?
#     beta0_sig <- as.logical(setup_vector[3])
#     #betas; which position indices are significant
#     sig_beta_positions <- which(as.logical(setup_vector[4:25]), arr.ind = TRUE)
# 
#     j <- length(unique(df$cusec) )
# 
#     #Iterate over 3 possible values of tau
#     data_list <- list(
#       #Dimensions
#       N = dim(X)[1],
#       j = j,
#       p=ncol(X),
#       #Data
#       X = X,
#       y = y,
#       #CAR
#       node1=node1,
#       node2=node2,
#       N_edges=N_edges,
#       #Indices
#       cusecs_indx = as.integer(subset(df, cob==country)$cusec_idx),
#       betas_to_estimate=as.integer(length(sig_beta_positions)),
#       parameter_indx = as.integer(sig_beta_positions),
#       sig_beta0=as.numeric(beta0_sig)
#     )
# 
#     model <-  stan("poisson_car.stan", data = data_list, chains=4, thin=2,
#                    iter=4000,  seed=1, init = 0)
#     save(model, file = paste0("./model_car/model_car_", country,"_tau",tau, ".Rdata"))
# 
#   }
# }



#############################
#5.Analysis
#############################

#5.1 Checking convergence

# for (country in c('Italia','Marruecos','China','Ecuador','Perú')) {
#   #Load all models
#   load(paste0("./stan_final_models/model_intercept/model_intercept_", country, ".Rdata"))
#   summary_output <- summary(model)
#   rhat_values <- summary_output$summary[, "Rhat"]
#   neff_values <- summary_output$summary[, "n_eff"]
#   print(country)
#   print(max(rhat_values))
#   print(min(neff_values))
#   load(paste0("./stan_final_models/model_car/model_car_", country,"_tau0.001.Rdata"))
#   summary_output <- summary(model)
#   rhat_values <- summary_output$summary[, "Rhat"]
#   neff_values <- summary_output$summary[, "n_eff"]
#   print(country)
#   print(max(rhat_values))
#   print(min(neff_values))
#   load(paste0("./stan_final_models/model_car/model_car_", country,"_tau0.1.Rdata"))
#   summary_output <- summary(model)
#   rhat_values <- summary_output$summary[, "Rhat"]
#   neff_values <- summary_output$summary[, "n_eff"]
#   print(country)
#   print(max(rhat_values))
#   print(min(neff_values))
#   load(paste0("./stan_final_models/model_car/model_car_", country,"_tau10.Rdata"))
#   summary_output <- summary(model)
#   rhat_values <- summary_output$summary[, "Rhat"]
#   neff_values <- summary_output$summary[, "n_eff"]
#   print(country)
#   print(max(rhat_values))
#   print(min(neff_values))
#   }
# 
#    #Save traceplot
# load(paste0("./stan_final_models/model_car/model_car_italia_tau0.1.Rdata"))
# traceplot(model)
# print(model)



#5.2 Compare LOOCV values

#for (country in c('Italia','Marruecos','China','Ecuador','Perú')) {
  # #Load all models
  # load(paste0("./stan_final_models/model_intercept/model_intercept_", country, ".Rdata"))
  # model0<-model
  # loo0 <- loo(model0, pars = "log_lik_fixed")
  # load(paste0("./stan_final_models/model_car/model_car_", country,"_tau0.001.Rdata"))
  # model1<-model
  # loo1 <- loo(model1, pars = "log_lik_fixed")
  # load(paste0("./stan_final_models/model_car/model_car_", country,"_tau0.1.Rdata"))
  # model2<-model
  # loo2 <- loo(model2, pars = "log_lik_fixed")
  # load(paste0("./stan_final_models/model_car/model_car_", country,"_tau10.Rdata"))
  # model3<-model
  # loo3 <- loo(model3, pars = "log_lik_fixed")
  # comp <- loo::loo_compare(loo0, loo1, loo2, loo3)
  # print(country)
  # print(loo0)
  # print(loo1)
  # print(loo2)
  # print(loo3)
  # print(comp, digits = 2)
  # print('\n')
  # }


#5.3 Coefficients for best models

# Function to get the estimate and CI
get_estimate <- function(fit) {
  # Extract the summary statistics
  summary_stats <- summary(fit)
  
  # Extract coefficients, mean, and confidence intervals
  coefficients <- rownames(summary_stats$summary)
  means <- summary_stats$summary[, "mean"]
  lower_ci <- summary_stats$summary[, "2.5%"]
  upper_ci <- summary_stats$summary[, "97.5%"]
  
  # Combine all information into a data frame
  coefficients_data <- data.frame(
    Coefficient = coefficients,
    Mean = means,
    Lower_CI = lower_ci,
    Upper_CI = upper_ci
  )
  
  return(coefficients_data)
}

#Function to get the fixed effect names
get_fe_names <- function(country, tau) {
  setup_vector<-binary_matrix[binary_matrix$country == country & binary_matrix$tau == tau,]
  first_row <- setup_vector[1, ]  
  selected_columns <- names(first_row)[first_row == TRUE]
  return(selected_columns)
}

#Italy
load(paste0("./stan_final_models/model_car/model_car_italia_tau0.1.Rdata"))
result <- get_estimate(model)
col_names<-get_fe_names('Italia','0.1')
result<- result[1:length(col_names),]
result$Coefficient<-col_names
result

#Marruecos
load(paste0("./stan_final_models/model_car/model_car_marruecos_tau0.1.Rdata"))
result <- get_estimate(model)
col_names<-get_fe_names('Marruecos','0.1')
result<- result[1:length(col_names),]
result$Coefficient<-col_names
result

#China
load(paste0("./stan_final_models/model_car/model_car_china_tau0.001.Rdata"))
result <- get_estimate(model)
col_names<-get_fe_names('China','0.001')
result<- result[1:length(col_names),]
result$Coefficient<-col_names
result

#Ecuador
load(paste0("./stan_final_models/model_car/model_car_ecuador_tau0.1.Rdata"))
result <- get_estimate(model)
col_names<-get_fe_names('Ecuador','0.1')
result<- result[1:length(col_names),]
result$Coefficient<-col_names
result

#Perú
load(paste0("./stan_final_models/model_car/model_car_peru_tau0.1.Rdata"))
result <- get_estimate(model)
col_names<-get_fe_names('Peru','0.1')
result<- result[1:length(col_names),]
result$Coefficient<-col_names
result


#5.4 Predictives
get_predictives <- function(model){
  #Fixed model
  generated_quantities <- as.array(model, "y_pred_fixed_effects")
  all_spatial_units <- dim(generated_quantities)[3]
  y_pred_all_units <- lapply(1:all_spatial_units, function(i) c(generated_quantities[, , i]))
  fixed <- sapply(y_pred_all_units, mean)
  #Fixed + ICAR model
  generated_quantities <- as.array(model, "y_pred_fixed_and_icar")
  all_spatial_units <- dim(generated_quantities)[3]
  y_pred_all_units <- lapply(1:all_spatial_units, function(i) c(generated_quantities[, , i]))
  fixed_icar <- sapply(y_pred_all_units, mean)
  #Full model
  generated_quantities <- as.array(model, "y_pred")
  all_spatial_units <- dim(generated_quantities)[3]
  y_pred_all_units <- lapply(1:all_spatial_units, function(i) c(generated_quantities[, , i]))
  full <- sapply(y_pred_all_units, mean)
  preds <- data.frame(cusec= unique(df$cusec), pred_fixed = fixed, pred_fixed_icar = fixed_icar, pred_full = full)
  return(preds)
}


#Italy
load(paste0("./stan_final_models/model_car/model_car_italia_tau0.1.Rdata"))
results <- get_predictives(model)
results$real <- subset(df, cob=='Italia')$population
results$census_total_pop <- subset(df, cob=='Italia')$census_total_pop
write.csv(results, file = "./predictives/preds_italy.csv", row.names = FALSE)


#Marruecos
load(paste0("./stan_final_models/model_car/model_car_marruecos_tau0.1.Rdata"))
results <- get_predictives(model)
results$real <- subset(df, cob=='Marruecos')$population
results$census_total_pop <- subset(df, cob=='Marruecos')$census_total_pop
write.csv(results, file = "./predictives/preds_morroco.csv", row.names = FALSE)

#China
load(paste0("./stan_final_models/model_car/model_car_china_tau0.001.Rdata"))
results <- get_predictives(model)
results$real <- subset(df, cob=='China')$population
results$census_total_pop <- subset(df, cob=='China')$census_total_pop
write.csv(results, file = "./predictives/preds_china.csv", row.names = FALSE)

#Ecuador
load(paste0("./stan_final_models/model_car/model_car_ecuador_tau0.1.Rdata"))
results <- get_predictives(model)
results$real <- subset(df, cob=='Ecuador')$population
results$census_total_pop <- subset(df, cob=='Ecuador')$census_total_pop
write.csv(results, file = "./predictives/preds_ecuador.csv", row.names = FALSE)

#Perú
load(paste0("./stan_final_models/model_car/model_car_peru_tau0.1.Rdata"))
results <- get_predictives(model)
results$real <- subset(df, cob=='Perú')$population
results$census_total_pop <- subset(df, cob=='Perú')$census_total_pop
write.csv(results, file = "./predictives/preds_peru.csv", row.names = FALSE)


ppc_dens_overlay(y, yrep_poisson[1:50, ])



#5.5 Posterior Predictive Check

#Perú
load(paste0("./stan_final_models/model_car/model_car_peru_tau0.1.Rdata"))
y <- subset(df, cob=='Perú')$population
xlim <- c(0, 200)

generated_quantities <- as.matrix(model, "y_pred")
plot <-ppc_dens_overlay(y, generated_quantities,  xlim=c(0, 200))
plot + xlim(xlim) 

generated_quantities <- as.matrix(model, "y_pred_fixed_effects")
plot <-ppc_dens_overlay(y, generated_quantities,  xlim=c(0, 200))
plot + xlim(xlim) 

dim(generated_quantities)



#############################
#6. Madrid
#############################

df <- read_csv("./intermediate_files/feature_matrix_madrid.csv")

#Target variable
df['pop_share'] <- df['population'] / df['census_total_pop']

#Unique CUSECS
j <- length(unique(df$cusec) )

#Scope of modeling is 5 countries: Italy, Morroco, China, Ecuador, Perú
df <- subset(df, cob %in% c('Italia','Marruecos','China','Ecuador','Perú'))

#Readjust cob_idx
df$cob_idx <- ifelse(df$cob == 'China',1,
                     ifelse(df$cob == 'Ecuador',2,
                            ifelse(df$cob == 'Italia',3,
                                   ifelse(df$cob == 'Marruecos',4,
                                          ifelse(df$cob == 'Perú',5,0)))))
i <- length(unique(df$cob) )

# Adjacency matrix
adj <- unlist(as.list(read.table("./intermediate_files/icar_adj_madrid.txt", header = FALSE)$V1))
# Neighborhood counts
num <- unlist(as.list(read.table("./intermediate_files/icar_num_madrid.txt", header = FALSE)$V1))

nbs = mungeCARdata4stan(adj, num);
node1 = nbs$node1;
node2 = nbs$node2;
N_edges = nbs$N_edges;

#Overall baseline by nationality
# Immigrant populations
df_total_im <- df %>%
  group_by(cob) %>%
  summarize(total_im = sum(population))
# Total Population
df_total_pop <- df %>%
  group_by(cusec) %>%
  summarize(total_pop = max(census_total_pop)) %>%
  summarize(total_pop = sum(total_pop))
#Overall density by nationality
result <- cross_join(df_total_im, df_total_pop)
result <- result %>%
  mutate(overall_density = total_im / total_pop)


####6.1. Horseshoe prior model
# for (country in c('Italia','Ecuador')) {
# 
#   X <- subset(df, cob==country)
#   X <- X[, !names(X) %in% c("cusec","cusec_idx", "cob","cob_idx","population","pop_share")]
#   X <- scale(X)
#   y <- subset(df, cob==country)$population
# 
#   tau <-0.1
#   
#   data_list <- list(
#     #Dimensions
#     N = dim(X)[1],
#     j = j,
#     p=ncol(X),
#     #Data
#     X = X,
#     y = y,
#     #CAR
#     node1=node1,
#     node2=node2,
#     N_edges=N_edges,
#     #Indices
#     cusecs_indx = as.integer(subset(df, cob==country)$cusec_idx),
#     #Penalty
#     tau=tau
#   )
# 
#   model <-  stan("./stan_final_models/poisson_horseshoe.stan", data = data_list, chains=4, thin=2,
#                  iter=4000,  seed=1, init = 0)
#   save(model, file = paste0("./stan_final_models/madrid/model_horseshoe_", country,"_tau0.1_madrid", ".Rdata"))
#   }


#Initiate a binary matrix where we'll indicate in each setup which are the significant variables
binary_matrix <- matrix(NA, nrow = 2, ncol = 25)
counter <- 1

for (country in c('Italia','Ecuador')) {
  tau_values <- c(0.1)
  for (tau in tau_values) {
    
    load(paste0("./stan_final_models/madrid/model_horseshoe_", country,"_tau0.1_madrid.Rdata"))
    
    summary_stats <- summary(model)
    posterior_samples <- rstan::extract(model)
    num_betas <- 1 + dim(posterior_samples$betas)[2]
    
    binary_matrix[counter,1] <- country
    binary_matrix[counter,2] <- tau
    
    #Check if intercept is significant
    quantiles <- quantile(posterior_samples$beta0, c(0.025, 0.975))
    binary_matrix[counter,3] <- !(quantiles[1] <= 0 && quantiles[2] >= 0)
    
    #Check if other betas are significant
    for (j in 1:22) {
      quantiles <- quantile(posterior_samples$betas[,j], c(0.025, 0.975))
      binary_matrix[counter,j+3] <- !(quantiles[1] <= 0 && quantiles[2] >= 0)
    }
    
    counter <- counter + 1
  }}

binary_matrix<-as.data.frame(binary_matrix)
colnames(binary_matrix) <- c(c('country','tau','beta0'),colnames(df[, !names(df) %in% c("cusec","cusec_idx", "cob","cob_idx","population","pop_share")]))

#How many significant variables are thewre per setup?
rowSums(as.matrix(binary_matrix[, 3:25]) == "TRUE")


####6.2 BYM model

# for (country in c('Italia','Ecuador')) {
#   X <- subset(df, cob==country)
#   X <- X[, !names(X) %in% c("cusec","cusec_idx", "cob","cob_idx","population","pop_share","overall_density")]
#   X <- scale(X)
#   y <- subset(df, cob==country)$population
# 
#   tau_values <- c(0.1)
#   for (tau in tau_values) {
# 
#     #Relevant setup
#     setup_vector<-binary_matrix[binary_matrix$country == country & binary_matrix$tau == tau,]
#     #beta0 is significant?
#     beta0_sig <- as.logical(setup_vector[3])
#     #betas; which position indices are significant
#     sig_beta_positions <- which(as.logical(setup_vector[4:25]), arr.ind = TRUE)
# 
#     j <- length(unique(df$cusec) )
# 
#     #Iterate over 3 possible values of tau
#     data_list <- list(
#       #Dimensions
#       N = dim(X)[1],
#       j = j,
#       p=ncol(X),
#       #Data
#       X = X,
#       y = y,
#       #CAR
#       node1=node1,
#       node2=node2,
#       N_edges=N_edges,
#       #Indices
#       cusecs_indx = as.integer(subset(df, cob==country)$cusec_idx),
#       betas_to_estimate=as.integer(length(sig_beta_positions)),
#       parameter_indx = as.integer(sig_beta_positions),
#       sig_beta0=as.numeric(beta0_sig)
#     )
# 
#     model <-  stan("./stan_final_models/poisson_car.stan", data = data_list, chains=4, thin=2,
#                    iter=4000,  seed=1, init = 0)
#     save(model, file = paste0("./stan_final_models/madrid/model_car_", country,"_tau",tau, "_madrid.Rdata"))
# 
#   }
# }

###6.3 Effects
#Italy
load(paste0("./stan_final_models/madrid/model_car_italia_tau0.1_madrid.Rdata"))
result <- get_estimate(model)
col_names<-get_fe_names('Italia','0.1')
result<- result[1:length(col_names),]
result$Coefficient<-col_names
result

#Ecuador
load(paste0("./stan_final_models/madrid/model_car_ecuador_tau0.1_madrid.Rdata"))
result <- get_estimate(model)
col_names<-get_fe_names('Ecuador','0.1')
result<- result[1:length(col_names),]
result$Coefficient<-col_names
result
  

###6.4 Predictives

#Ecuador
load(paste0("./stan_final_models/madrid/model_car_ecuador_tau0.1_madrid.Rdata"))
results <- get_predictives(model)
results$real <- subset(df, cob=='Ecuador')$population
results$census_total_pop <- subset(df, cob=='Ecuador')$census_total_pop
write.csv(results, file = "./predictives/preds_ecuador_madrid.csv", row.names = FALSE)
