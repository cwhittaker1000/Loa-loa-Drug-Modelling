# Import dataframe for use in likelihood function
individual_ivermectin_data$Timesteps <- individual_ivermectin_data$Time * 4

# Model Related Functions- Gives You "IVM_efficacy_model" (base structure) & "run_model" (runs the model)
source("sansL3_Ivermectin_Modelling.R")

# Prior Function
prior_function <- function(parameters, fitting_yn) {
  
  res <- 0
  if (fitting_yn["foi"]) {
    foi <- parameters[1]
    res <- res + dunif(foi, min = 0.000001, max = 1, log = T)
  }
  if (fitting_yn["epsilon"]) {
    epsilon <- parameters[2]
    res <- res + dnorm(epsilon, mean = 11.689, sd = 1, log = T)
  }
  if (fitting_yn["gamma"]) {
    gamma <- parameters[3]
    res <- res + dnorm(gamma, mean = 0.004014, sd = 0.0003, log = T)
  }
  if (fitting_yn["muMF_0"]) {
    muMF_0 <- parameters[4]
    res <- res + dnorm(muMF_0, mean = 0.006157, sd = 0.0005, log = T)
  }
  if (fitting_yn["mu_W0"]) {
    mu_W0 <- parameters[5]
    res <- res + dunif(mu_W0, min = 0.00001, max = 0.10000, log = T)
  }
  if (fitting_yn["beta_0"]) {
    beta_0 <- parameters[6]
    res <- res + dnorm(beta_0, mean = 0.2, sd = 0.05, log = T)
  }
  if (fitting_yn["mu_W1"]) {
    mu_W1 <- parameters[7]
    res <- res + dnorm(mu_W1, mean = 0.2, sd = 0.05, log = T)
  }
  if (fitting_yn["muMF_1_max"]) {
    muMF_1_max <- parameters[8]
    res <- res + dunif(muMF_1_max, min = 0.001, max = 19, log = T)
  }
  if (fitting_yn["rho"]) {
    rho <- parameters[9]
    res <- res + dunif(rho, min = 0.01, max = 10, log = T)
  }
  if (fitting_yn["beta_1_max"]) {
    beta_1_max <- parameters[10]
    res <- res + dunif(beta_1_max, min = 0.001, max = 19, log = T)
  }
  if (fitting_yn["lambda"]) {
    lambda <- parameters[11]
    res <- res + dunif(lambda, min = 0, max = 10, log = T)
  }
  
  # Return sum of prior probability densities
  prior_values <- unname(res)
  return(prior_values)
  
}


# Likelihood Function
log_lik_function <- function(parameters, treatment_time, recorded_data) {
  
  # Extracting Relevant Variables from Data
  follow_up_times <- recorded_data$Time
  #print(follow_up_times)
  
  recorded_mf_counts <- recorded_data$Mean_MF_Levels
  #print(recorded_mf_counts)
  
  number_individuals <- recorded_data$Number_Individuals
  #print(number_individuals)
  
  standard_error <- recorded_data$StdErr
  #print(standard_error)

  # Runs the Model Using Supplied Parameters and Extracts Entire MF Output
  out <- run_model(parameters, treatment_time, TRUE) 
  mf_counts <- out[, "MF"]
  
  # Extracts Model MF Output for Timepoints We Have Data For
  model_predicted_mf_counts <- mf_counts[(follow_up_times/0.1 + (treatment_time/0.1))]

  #print(follow_up_times/0.1 + (treatment_time/0.1))
  #print(model_predicted_mf_counts)

  # Calculates the Likelihood of the Data Given the Model MF Predictions (Average So Assume Normal Dist)
  individual_log_liks <- dnorm(recorded_mf_counts, model_predicted_mf_counts, sd = standard_error, log = T)
  #print(individual_log_liks)
  summed_log_lik <- sum(individual_log_liks)
  
  # Returns the Complete Loglikelihood
  return(summed_log_lik)
  
}


# Posterior Function
posterior_function <- function(parameters, fitting_yn, treatment_time, recorded_data) {
  
  # Calculate the Posterior Density
  posterior_output <- log_lik_function(parameters, treatment_time, recorded_data) +
                      prior_function(parameters, fitting_yn)
  # Return Posterior Density 
  return(posterior_output)

}

# Covariance Matrix Adaptation Function
joint_proposal_SD_adapter_R <- function(accepted_variable, current_iteration, iteration_adapting_began,
                                        current_scaling_factor, mu_previous, current_parameter_values,
                                        current_covariance_matrix) {
  
  # Defining various quantities
  iterations_since_adapting_began <- current_iteration - iteration_adapting_began
  cooldown <- (iterations_since_adapting_began + 1)^-0.99
  
  # Calculating the various components
  new_correlation_matrix <- ((1 - cooldown) * current_covariance_matrix) + (cooldown * ((t(current_parameter_values - mu_previous)) %*%  (current_parameter_values - mu_previous)))
  new_mu <- ((1 - cooldown) * mu_previous) + (cooldown * current_parameter_values)
  log_new_scaling_factor <- log(current_scaling_factor) + cooldown * (accepted_variable - 0.25)
  new_scaling_factor <- exp(log_new_scaling_factor)
  new_covariance_matrix <- new_scaling_factor * new_correlation_matrix
  
  # Outputting the relevant results
  output_list <- list()
  output_list[["New_Covariance_Matrix"]] = new_covariance_matrix
  output_list[["New_Mu"]] = new_mu
  output_list[["New_Scaling_Factor"]] = new_scaling_factor
  return(output_list)
  
} 
  
# Proposal Function
proposal_function <- function(parameters, fitting_yn, sigma) {
  
  eigenvalues <- eigen(sigma)
  minimum <- min(eigenvalues$values)
  
  if (minimum <= 0.0) {
    for (i in 1:sum(fitting_yn)) {
      for (j in 1:sum(fitting_yn)) {
        if (i == j) {
          sigma[i, j] = sigma[i, j] + 0.01;
        }
      }
    }
  }

  mu <- c()
  index <- 1
  for (i in 1:sum(fitting_yn)) {
      mu[index] <- parameters[i]
      index <- index + 1
  }
  ncols <- length(mu)
  output <- mu + rnorm(ncols) %*% chol(sigma)
  for (i in 1:length(mu)) {
    if (output[i] <= 0) {
      output[i] <- 0.00001
    }
  }
  return(output)
}

# Run MCMC Function 
MCMC_adapt <- function(start_adaptation, end_adaptation, fitting_yn, parameters, 
                       number_of_iterations, initial_sds, treatment_time, recorded_data) {
  
  # Storage for Output
  MCMC_output <- matrix(nrow = number_of_iterations + 1, ncol = sum(fitting_yn))
  Acceptances <- vector(length = number_of_iterations + 1)
  Acceptance_Ratio <- vector(length = number_of_iterations + 1)
  Accepted <- 0
  
  # Assigning Initial Values to First Row of MCMC_Output and Naming The Columns
  index <- 1
  name_vector <- c()
  for (i in 1:length(fitting_yn)) {
    if (fitting_yn[i] == T) {
      name_vector[index] <- names(fitting_yn[i])
      MCMC_output[1, index] <- parameters[i]
      index <- index + 1
    }
  }
  colnames(MCMC_output) <- name_vector
  
  # Components for Adaptive MCMC   
  current_scaling_factor <- 1
  current_mu <- matrix(nrow = 1, ncol = sum(fitting_yn))
  current_mu[1, ] <- MCMC_output[1, ]
  
  # Creating the Covariance Matrix
  current_covariance_matrix <- matrix(data = 0, nrow = sum(fitting_yn), ncol = sum(fitting_yn))
  tracker <- 1
  for (i in 1:length(fitting_yn)) {
    for (j in 1:length(fitting_yn)) {
      if (fitting_yn[i] == T) {
        if(i == j) {
          current_covariance_matrix[tracker, tracker] = initial_sds[i]
          tracker <- tracker + 1
        } 
      }
    }
  }

  # Calculating posterior density for initial values
  MCMC_parameters <- vector(mode = "numeric", length = length(fitting_yn))
  tracking_index <- 1
  for(i in 1:length(fitting_yn)) {
    if (fitting_yn[i] == T) {
      MCMC_parameters[i] <- MCMC_output[1, tracking_index] 
      tracking_index <- tracking_index + 1
    }
    else {
      MCMC_parameters[i] <- parameters[i]
    }
  }
  names(MCMC_parameters) <- names(parameters)
  current_posterior <- posterior_function(test_parameters, fitting_yn, treatment_time, recorded_data)
  if(is.nan(current_posterior)) {
    return("ERROR- INITIAL PARAMETERS ARE RETURNING NAN")
  }

  # Running the Actual MCMC
  for (i in 1:number_of_iterations){
    #print(i)
    
    # Proposing new values and calculating relevant posterior density
    #print(MCMC_output[i, ])
    proposed_parameter_values <- proposal_function(MCMC_output[i, ], fitting_yn, current_scaling_factor * current_covariance_matrix)
    #print(proposed_parameter_values)

    proposed_MCMC_parameters <- vector(mode = "numeric", length = length(fitting_yn))
    tracking_index_proposed <- 1
    for(j in 1:length(fitting_yn)) {
      if (fitting_yn[j] == T) {
        proposed_MCMC_parameters[j] <- proposed_parameter_values[tracking_index_proposed] 
        tracking_index_proposed <- tracking_index_proposed + 1
      }
      else {
        proposed_MCMC_parameters[j] <- parameters[j]
      }
    }
    names(proposed_MCMC_parameters) <- names(parameters)
    #print(proposed_MCMC_parameters)
    proposed_posterior <- posterior_function(proposed_MCMC_parameters, fitting_yn, treatment_time, recorded_data)
    #print(proposed_posterior)
    if(is.nan(proposed_posterior)) {
      proposed_posterior <- -Inf
    }
    #print(proposed_posterior)
    
    # Calculate likelihood ratio and evaluate acceptance/rejection
    likelihood_ratio <- exp(proposed_posterior - current_posterior);
    #print(likelihood_ratio)
    if(runif(1) < likelihood_ratio) {
      #print(proposed_parameter_values)
      MCMC_output[i + 1, ] <- proposed_parameter_values
      #print(MCMC_output[1:(i+1), ])
      Acceptances[i] <- 1
      Accepted <- 1
      current_posterior <- proposed_posterior
    } 
    else {
      MCMC_output[i + 1, ] <- MCMC_output[i, ]
      Acceptances[i] <- 0
      Accepted <- 0
      current_posterior <- current_posterior
    }
    Acceptance_Ratio[i] <- sum(Acceptances)/i
    
    # Adaptation of the Covariance Matrix Associated With Proposal Function
    if (i > start_adaptation & i < end_adaptation) {
      
      latest_parameter_values <- MCMC_output[i + 1, ]
      adapter_output <- joint_proposal_SD_adapter_R(Accepted, i, start_adaptation, current_scaling_factor, 
                                                    current_mu, latest_parameter_values, current_covariance_matrix)
      current_mu <- adapter_output[["New_Mu"]];
      #print(current_mu)
      current_covariance_matrix <- adapter_output[["New_Covariance_Matrix"]];
      #print(current_covariance_matrix)
      current_scaling_factor <- adapter_output[["New_Scaling_Factor"]];
      #print(current_scaling_factor)
      
    }
    
    # Outputting Results as the MCMC Runs
    if(i %% 100 == 0 | i == 1 | i == 2 | i == 3 | i == 4 | i == 5 | i == 6 | i == 7 | i == 8 | i == 9 | i == 10) {
      print(c("The iteration number is", i))
      print(c("The acceptance ratio is", Acceptance_Ratio[i]))
      print(c("Current mu is:", current_mu)) 
      print(c("The current covariance matrix \n", current_covariance_matrix))
    }
  }
  
  # Create List of Outputs and Return It
  list <- list(); 
  list[["MCMC_Output"]] <- MCMC_output; 
  list[["Acceptances"]] <- Acceptances
  return(list)
  
}

##############################################################################################################
### 
### Tester Code for Various functions
###
##############################################################################################################

fitting_yn <- c("foi" = T, "epsilon" = T, "gamma" = T, "muMF_0" = T, "mu_W0" = T,
                "beta_0" = F, "mu_W1" = F, "muMF_1_max" = T, "rho" = T, "beta_1_max" = T,
                "lambda" = T)
test_parameters <- c(foi = 0.0035, epsilon = 11.689, gamma = 0.004014,
                     muMF_0 = 0.006157, mu_W0 = 0.0002739, beta_0 = 0, mu_W1 = 0,
                     muMF_1_max = 2, rho = 1,
                     beta_1_max = 5, lambda = 0.1)

MCMC_adapt <- function(start_adaptation, end_adaptation, fitting_yn, parameters, 
                       number_of_iterations, initial_sds, treatment_time, recorded_data)

initial_sds <- c(0.0001, 0.1, 0.00005, 0.0001, 0.000004, 0, 0, 0.03, 0.03, 0.05, 0.003)
  
bloop <- MCMC_adapt(start_adaptation = 10, end_adaptation = 1000, fitting_yn = fitting_yn, 
                    parameters = test_parameters, number_of_iterations = 1000, initial_sds = initial_sds, 
                    treatment_time = 10, recorded_data = individual_ivermectin_data)

plot(bloop$MCMC_Output[, 1], type = "l")


# Tester Code for Prior Function
prior_function(test_parameters, fitting_yn) 

# Tester Code for Likelihood Function
log_lik_function(test_parameters, 10, individual_ivermectin_data)

# Tester Code for Posterior Function
posterior_function(test_parameters, fitting_yn, treatment_time = 10, individual_ivermectin_data)

# Tester Code for Proposal Function
parameters <- c(100, 200, 500)
fitting_yn_proposal_test <- c(T, F, T)
sigma <- matrix(c(10, 5, 5, 10), nrow = 2, ncol = 2, byrow = 2)
proposal_function(parameters, fitting_yn_proposal_test, sigma)

storage <- matrix(nrow = 10000, ncol = 2)
for (i in 1:10000) {
  storage[i, ] <- proposal_function(parameters, fitting_yn_proposal_test, sigma)
}
hist(storage[, 2])
plot(storage[, 1], storage[, 2])