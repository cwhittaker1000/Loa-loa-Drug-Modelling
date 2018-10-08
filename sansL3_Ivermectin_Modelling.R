# Loading required packages to run and fit the model
require(deSolve)
library(tictoc)
Pion_Overall <- read.csv("~/13th April 2018/Loa Loa/EpiLOA Data Analysis/Parameter Estimation and Model Building/IVM & ALB Efficacy Estimation/Pion_2018/Pion_Overall.csv")
individual_ivermectin_data <- Pion_Overall

# CONSIDER REMOVING THE LARVAL COMPARTMENTS AS IF WE JUST ASSUME PEOPLE'RE AT EQUILIBRIUM 
# THEN WE CAN PROBABLY JUST IGNORE IT. BUT THAT ASSUMES THAT L3 LARVAE DEVELOPING INTO ADULT WORMS
# DON'T EXPERIENCE ANY EFFECTS FROM TREATMENT. WHICH IS MAYBE FINE IF WE'RE ASSUMING THAT THE DRUG
# ONLY AFFECTS FERTILITY, NOT MORTALITY RATE. 

## Initial condition for:
# 1 F compartment
# 1 N compartment
# 1 MF compartment

# Model of Infection & MF Production
IVM_efficacy_model <- function(times, y, parameters, treatment_time, microfilaricidal){
  
  with(as.list(c(y, parameters)), {
    
    yvec <- y
    dyvec <- rep(0, 3)
    
    for (i in 1:3) { 
      
      if (times < treatment_time) {
        
        if (i == 1) { # i.e. fertile worm compartment
          dyvec[i] = (0.5 * foi) - (beta_0 * yvec[i]) - (mu_W0 * yvec[i])
        } 
        else if (i == 2) { # i.e. infertile worm compartment
          dyvec[i] = beta_0 * yvec[i - 1] - (mu_W0 * yvec[i])
        }
        else if (i == 3) { # i.e. microfilarie compartment
          dyvec[i] = (epsilon * yvec[i - 2]) - (muMF_0 * yvec[i])
        }
      }
      
      else if (times >= treatment_time) {
        
        t <- times - treatment_time 
        
        if (microfilaricidal) {
          muMF_1 <- muMF_1_max * exp(t * -rho)
        }
        else {
          muMF_1 <- 0
        }
        
        beta_1 <- beta_1_max * exp(-lambda * t)
        
        if (i == 1) { # i.e. fertile worm compartment
          dyvec[i] = (0.5 * foi) - (beta_0 + beta_1) * yvec[i] - (mu_W0 + mu_W1) * yvec[i]
        } 
        else if (i == 2) { # i.e. infertile worm compartment
          dyvec[i] = (beta_0 + beta_1) * yvec[i - 1] - (mu_W0 + mu_W1) * yvec[i]
        } 
        else if (i == 3) { # i.e. microfilarie compartment
          dyvec[i] = (epsilon * yvec[i - 2]) - (muMF_0 + muMF_1) * yvec[i] 
        }
      }
    }
    
    return(list(rbind(dyvec)))
  })
}

# Function for Running the Model & Outputting the Results
run_model <- function(parameters, treatment_time, microfilaricidal) {
  
  # Initial Values
  Fertile_Worms <- unname(parameters["foi"] / (2 * (parameters["beta_0"] + parameters["mu_W0"])))
  MF <- unname((parameters["epsilon"] * Fertile_Worms) / parameters["muMF_0"])
  
  initial_conditions <- c(Fertile_Worms, 0, MF)
  
  # Running the Model For Defined Time Period & Stepsize
  stepsize <- 0.1
  times <- seq(0, 370, stepsize)
  out <- rk4(y = initial_conditions, times = times, func = IVM_efficacy_model, parms = parameters, treatment_time = treatment_time, microfilaricidal = microfilaricidal)
  
  # Process & Store the Data (note 1st column is time so all columns shifted 1 compared to the model)
  Fertile_Worms <- out[, 2]
  Infertile_Worms <- out[, 3]
  MF <- out[, 4]
  t <- out[, 1]
  
  cbind(t, Fertile_Worms, Infertile_Worms, MF)
}

# Model Test Running
# Note: Weird behaviour with certain parameter values. 
# Omega can't really go less than -2 and shouldn't be above 0.
# Rho can't be negative, but other than that, we're alright (1000 means basically increased death rate doesn't decline)
# Lambda shouldn't be negative, as it makes the effect size increase over time.
# Beta_1_max shouldn't go above 24 for a lambda of 1 up to 1000.

### NOTES: TRANSFORM THE PARAMETER VALUES TO GET OUT OF THE DANGERIOUS TERRITORY WHERE IT RETURNS WEIRD VALUES <0 
### NOTES: ASK MARTIN IF THERE'S A WAY OF CURTAILING IT AT ZERO OR SOMETHING?
### NOTES: AVOIDED THROUGH DIFFERENT PARAMETERISATION
### NOTES: CAN I SOLVE FOR INITIAL CONDITIONS IN SUCH A WAY THAT IT MEANS I DON'T HAVE TO RUN TO EQUILIBRIUM AS LONG?
### NOTES: CHECK UNITS - I.E. MAKE SURE EVERYTHING'S IN DAYS, OR IN WEEKS, AS APPROPRIATE!!!

test_parameters <- c(foi = 0.0035, epsilon = 11.689, gamma = 0.004014,
                     muMF_0 = 0.006157, mu_W0 = 0.0002739, beta_0 = 0, mu_W1 = 0,
                     muMF_1_max = 2, rho = 0.1,
                     beta_1_max = 0.2, lambda = 0.1)

tic()
out <- run_model(MCMC_parameters, 10, TRUE)
toc()
plot(out[, "t"], out[, "MF"], type = "l", col = "red",  xlim = c(0, 400), ylim = c(0, max(out[, "MF"]) + 2000), lwd = 3)
points(10 + (individual_ivermectin_data$Time), individual_ivermectin_data$Mean_MF_Levels, pch = 20)



out[50:103, ]
individual_ivermectin_data$Time + (10/0.1)
out[individual_ivermectin_data$Timesteps + (10/0.1), ]


plot(out[, "t"], out[, "Fertile_Worms"], type = "l")

# Profile of microfilaricidal effect (two options for parameterisation, currently using the first one)
timepoints <- seq(0, 365, 1)
plot(timepoints, test_parameters["muMF_1_max"]* exp(timepoints * test_parameters["rho"]), type = "l", xlim = c(0, 10), ylab = "")

# Profile of macrofilaricidal effect
plot(timepoints, test_parameters["beta_1_max"] * exp(-test_parameters["lambda"] * timepoints), type = "l")


plot(out[, "t"], out[, "Fertile_Worms"], type = "l")
plot(out[, "t"], out[, "Infertile_Worms"], type = "l")

# Plotting Beta1 Profile
timepoints <- seq(0, 100, 1)
bloop <- test_parameters["beta_1_max"] * exp(-test_parameters["lambda"] * timepoints)
plot(timepoints, bloop, type = "l")

# Plotting muMF_1 Profile
timepoints <- seq(0, 100, 0.01)
sloop <- (timepoints + test_parameters["rho"]) ^ test_parameters["omega"]
plot(timepoints, sloop, type = "l", xlim = c(0, 10))

# Alternative parameterisation for muMF_1
# Issue - it doesn't always decline as fast at the inverse function above, BUT
# with the right prior (uniform, constraining the values it can take), we can make it so that
# the MCMC will favour values that produce really fast dropping curves. It's also a lot more stable than the 
# parameterisation that Maria-Gloria used and hence I'm going to go with this instead. 
plot(timepoints, test_parameters["muMF_1_max"]*(timepoints + 1)^test_parameters["rho"], type = "l", xlim = c(0, 10))

# Initial Conditions Equations
Fertile_Worms <- test_parameters["foi"] / (2 * (test_parameters["beta_0"] + test_parameters["mu_W0"]))
MF <- (test_parameters["epsilon"] * Fertile_Worms) / test_parameters["muMF_0"]

