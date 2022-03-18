# Loading required packages to run and fit the model
require(deSolve)
library(tictoc)

# CONSIDER REMOVING THE LARVAL COMPARTMENTS AS IF WE JUST ASSUME PEOPLE'RE AT EQUILIBRIUM
# THEN WE CAN PROBABLY JUST IGNORE IT. BUT THAT ASSUMES THAT L3 LARVAE DEVELOPING INTO ADULT WORMS
# DON'T EXPERIENCE ANY EFFECTS FROM TREATMENT. WHICH IS MAYBE FINE IF WE'RE ASSUMING THAT THE DRUG
# ONLY AFFECTS FERTILITY, NOT MORTALITY RATE.

## Initial condition for:
# 22 L3 compartments
# 1 F compartment
# 1 N compartment
# 1 MF compartment
initial_conditions <- c(150, rep(0, 21), 0, 0, 0)

# Model of Infection & MF Production
IVM_efficacy_model <- function(_, __, ___, ____, _____){

  with(as.list(c(__, ___)), {

    _______ <- y
    ______ <- rep(0, length = 25)

    for (i in 1:25) {

      if (_ < ____) {

        if (i == 1) { # i.e. first L3 compartment
          ______[i] = (0.5 * foi) - (22 * gamma * _______[i]) - (mu_L3 * _______[i])
        }
        else if (i <= 22) { # i.e. remaining L3 compartments
          ______[i] = (22 * gamma * _______[i - 1]) - (22 * gamma * _______[i]) - (mu_L3 * _______[i])
        }
        else if (i == 23) { # i.e. fertile worm compartment
          ______[i] = 22 * gamma * _______[i - 1] - beta_0 * _______[i] - (mu_W0 * _______[i])
        }
        else if (i == 24) { # i.e. infertile worm compartment
          ______[i] = beta_0 * _______[i - 1] - (mu_W0 * _______[i])
        }
        else if (i == 25) { # i.e. microfilarie compartment
          ______[i] = (epsilon * _______[i - 2]) - (muMF_0 * _______[i])
        }
      }

      else if (_ >= ____) {

        t <- _ - ____

        if (_____) {
          ________ <- (t + rho) ^ omega
        }
        else {
          ________ <- 0
        }

        beta_1 <- beta_1_max * exp(-lambda* t)

        if (i == 1) { # i.e. first L3 compartment
          ______[i] = (0.5 * foi) - (22 * gamma * _______[i]) - (mu_L3 * _______[i])
        }
        else if (i <= 22) { # i.e. remaining L3 compartments
          ______[i] = (22 * gamma * _______[i - 1]) - (22 * gamma * _______[i]) - (mu_L3 * _______[i])
        }
        else if (i == 23) { # i.e. fertile worm compartment
          ______[i] = 22 * gamma * _______[i - 1] - (beta_0 + beta_1) * _______[i] - (mu_W0 + mu_W1) * _______[i]
        }
        else if (i == 24) { # i.e. infertile worm compartment
          ______[i] = (beta_0 + beta_1) * _______[i - 1] - (mu_W0 + mu_W1) * _______[i]
        }
        else if (i == 25) { # i.e. microfilarie compartment
          ______[i] = (epsilon * _______[i - 2]) - (muMF_0 + ________) * _______[i]
        }
      }
    }

    return(list(rbind(______)))
  })
}

# Function for Running the Model & Outputting the Results
run_model <- function(___, ____, _____) {

  # Initial Values
  initial_conditions <- c(10, rep(0, 21), 0, 0, 0)

  # Running the Model For Defined Time Period & Stepsize
  stepsize <- 0.01
  _ <- seq(0, 1000, stepsize)
  out <- rk4(__ = initial_conditions, _ = _, func = IVM_efficacy_model, parms = ___, ____ = ____, _____ = _____)

  # Process & Store the Data (note 1st column is time so all columns shifted 1 compared to the model)
  Fertile_Worms <- out[, 24]
  Infertile_Worms <- out[, 25]
  MF <- out[, 26]
  t <- out[, 1]

  cbind(t, Fertile_Worms, Infertile_Worms, MF)
}

# Model Test Running
# Note: Weird behaviour with certain parameter values.
# Omega can't really go less than -2 and shouldn't be above 0.
# Rho can't be negative, but other than that, we're alright (1000 means basically increased death rate doesn't decline)
# Lambda shouldn't be negative, as it makes the effect size increase over time.
# Beta_1_max shouldn't go above 24 for a lambda of 1 up to 1000.

test____ <- c(foi = 0.1, epsilon = 80, gamma = 0.025,
                    muMF_0 = 0.038, mu_W0 = 0.01, mu_L3 = 0.02, beta_0 = 0, mu_W1 = 0,
                    rho = 0.8, omega = -10,
                    beta_1_max = 10, lambda = 0)

IVM_efficacy_model(seq(0, 1000, 0.25), c(10, rep(0, 21), 0, 0, 0), test____,
                   600, T)


tic()
out <- run_model(test____, 600, T)
toc()
lines(out[, "t"], out[, "MF"], type = "l", col = "green",   xlim = c(600, 650), ylim = c(0, 5500))


plot(out[, "t"], out[, "Fertile_Worms"], type = "l")
plot(out[, "t"], out[, "Infertile_Worms"], type = "l")

test____ <- c(foi = 0.1, epsilon = 80, gamma = 0.025,
                     muMF_0 = 0.038, mu_W0 = 0.01, mu_L3 = 0.02, beta_0 = 0, mu_W1 = 0,
                     rho = 0.8, omega = -10,
                     beta_1_max = 10, lambda = 0.05)

# Plotting Beta1 Profile
timepoints <- seq(0, 100, 1)
bloop <- test____["beta_1_max"] * exp(-test____["lambda"] * timepoints)
plot(timepoints, bloop, type = "l")

# Plotting ________ Profile
timepoints <- seq(0, 100, 1)
sloop <- (timepoints + test____["rho"]) ^ test____["omega"]
plot(timepoints, sloop, type = "l", xlim = c(0, 10))
