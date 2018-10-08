# Standard Deviation Adapter for MCMCs Where You're Jointly Proposing Your Parameters


# IMPORTANT NOTE: MU_PREVIOUS, CURRENT_PARAMETER_VALUES AND CURRENT COVARIANCE MATRICES MUST ALL BE IN THE
# FORM OF R MATRICES FOR THIS TO WORK. The required dimensions of these matrices are as follows (rows, columns):

# Mu_Previous = (1, number_of_parameters) e.g. for 3 parameters it could be [3, 7, 9] (1 row, 3 columns)
# Current_Parameter_values = (1, number_of_parameters) e.g. for 3 parameters it could be [3, 7, 9] (1 row, 3 columns)
# Current_Covariance_Matrix = (number_of_parameters, number_of_parameters) e.g. for 3 parameters it'd be a 3 x 3 matrix.

# Accepted Variable = either 0 or 1. 0 if you rejected the most recently proposed parameter.
#                                    1 if you accepted the most recently proposed parameter.

# Current Iteration = the current iteration number the MCMC is on.

# Iteration Adapting Began = the iteration number that you started adapting on. I typically start
#                            start adapting from the beginning although you can start later on if required.
#                            (see Mu Previous for instances where you might need to)

# Current Scaling Factor = the value of "Scaling Factor" produced by the last run of this function.
#                          Need to provide an initial scaling factor in this case- this is just 1.

# Mu Previous = the value of "New Mu" produced by the last run of this function.
#               Need to provide initial values; typically I just assign whatever my initial values are
#               as the initial values. But this comes with the caveat that these initial values have to
#               be vaguely in the right ballpark. If they're unlikely to, consider running the MCMC for a
#               few thousand iterations, then setting the parameter values at that point as the initial
#               values for mu, and THEN start adapting.

# Current Parameter Values = just the current parameter values in your MCMC chain.

# Current Covariance Matrix = the output "New Covariance Matrix" from the previous run of this function.
#                             For the initial covariance matrix, just do a diagonal one (0s on non-diagonal elements)
#                             with a vaguely sensible set of variances along the diagonal.


joint_proposal_SD_adapter_R <- function(accepted_variable, current_iteration, iteration_adapting_began,
                                        current_scaling_factor, mu_previous, current_parameter_values,
                                        current_covariance_matrix) {

  # Calculating the number of iterations since you started adapting, and from that the cooldown factor.
  # Basically this algorithm adapts the covariance matrix specifying the size and shape of your jumps.
  # Over time however, it changes the covariance matrix less and less, so that eventually you get an
  # unchanging covariance matrix. This is what the cooldown factor is.
  iterations_since_adapting_began <- current_iteration - iteration_adapting_began
  cooldown <- (iterations_since_adapting_began + 1)^-0.6

  # Calculating the various components
  new_correlation_matrix = ((1 - cooldown) * current_covariance_matrix) + (cooldown * ((t(current_parameter_values - mu_previous)) %*%  (current_parameter_values - mu_previous)))
  new_mu = ((1 - cooldown) * mu_previous) + (cooldown * current_parameter_values)
  log_new_scaling_factor = log(current_scaling_factor) + cooldown * (accepted_variable - 0.25)
  new_scaling_factor = exp(log_new_scaling_factor)
  new_covariance_matrix = new_scaling_factor * new_correlation_matrix

  # Outputting the relevant results
  output_list <- list()
  output_list[[1]] = new_covariance_matrix
  output_list[[2]] = new_mu
  output_list[[3]] = new_scaling_factor

  return(output_list)


}


# # Tester Code
# 
# # Assigning initial values for mu and current parameters- currently vectors, need them to be matrices
# mu_prev_vec <- c(7, 3, 6)
# current_vec <- c(10, 5, 10)
# 
# # Converting the above to matrices
# mu_prev <- t(as.matrix(mu_prev_vec))
# current <- t(as.matrix(current_vec))
# 
# # Creating and filling a covariance matrix called bloop
# bloop <- matrix(nrow = 3, ncol =3)
# for (i in 1:3) {
#   for (j in 1:3) {
#     if (i == j) {
#       bloop[i, j] = 0.9
#     }
#     else {
#       bloop[i, j] = 0.2
#     }
#   }
# }
# 
# # Running the SD adapter using the values we have above- check the output!
# kloop <- joint_proposal_SD_adapter_R(1, 15, 10, log(2), mu_prev, current, bloop)
# 
# 
# 
# is.matrix(mu_prev)
# 
# 
# bloop
# sloop$New_Covariance_Matrix
# sloop$New_Mu









