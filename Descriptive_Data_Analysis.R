# Preliminary Data Analysis and Plotting 

# Seb's Data 

## Overall Study 
Pion_Overall <- read.csv("~/13th April 2018/Loa Loa/EpiLOA Data Analysis/Parameter Estimation and Model Building/IVM & ALB Efficacy Estimation/Pion_2018/Pion_Overall.csv")
names(Pion_Overall) <- c("Time", "Mean_MF_Loads", "Number_Individuals")
weights <- Pion_Overall$Number_Individuals/sum(Pion_Overall$Number_Individuals)
for (i in 1:length(Pion_Overall$Time)) {
  if (weights[i] < 0.05) {
    weights[i] <- 0.05
  }
}
plot(Pion_Overall$Time, Pion_Overall$Mean_MF_Loads, cex = 2, ylim = c(0, 15000), pch = 20)
plot(log(Pion_Overall$Time + 1), Pion_Overall$Mean_MF_Loads, cex = weights*15, ylim = c(0, 15000), pch = 20)

## 200mg
Pion_200 <- read.csv("~/13th April 2018/Loa Loa/EpiLOA Data Analysis/Parameter Estimation and Model Building/IVM & ALB Efficacy Estimation/Pion_2018/Pion_200mg.csv")
names(Pion_200) <- c("Time", "Mean_MF_Loads", "Number_Individuals")
weights <- Pion_200$Number_Individuals/sum(Pion_200$Number_Individuals)
for (i in 1:length(Pion_200$Time)) {
  if (weights[i] < 0.05) {
    weights[i] <- 0.05
  }
}
plot(Pion_200$Time, Pion_200$Mean_MF_Loads, cex = 2, ylim = c(0, 15000), pch = 20)
plot(log(Pion_200$Time), Pion_200$Mean_MF_Loads, cex = 2, ylim = c(0, 15000), pch = 20)

## 150mg
Pion_150 <- read.csv("~/13th April 2018/Loa Loa/EpiLOA Data Analysis/Parameter Estimation and Model Building/IVM & ALB Efficacy Estimation/Pion_2018/Pion_150mg.csv")
names(Pion_150) <- c("Time", "Mean_MF_Loads", "Number_Individuals")
weights <- Pion_150$Number_Individuals/sum(Pion_150$Number_Individuals)
for (i in 1:length(Pion_150$Time)) {
  if (weights[i] < 0.05) {
    weights[i] <- 0.05
  }
}
plot(Pion_150$Time, Pion_150$Mean_MF_Loads, cex = 2, ylim = c(0, 15000), pch = 20)
plot(log(Pion_150$Time), Pion_150$Mean_MF_Loads, cex = 2, ylim = c(0, 15000), pch = 20)

# Plotting the Raw Counts for Each Treatment Categorisation
plot(Pion_Overall$Time, Pion_Overall$Mean_MF_Loads, type = "l", ylim = c(0, max(Pion_Overall$Mean_MF_Loads)), lwd = 2)
plot(log(Pion_Overall$Time + 1), Pion_Overall$Mean_MF_Loads, type = "l", ylim = c(0, 15000), lwd = 2)

Non_NA_200 <- !is.na(Pion_200$Mean_MF_Loads)
plot(Pion_200$Time[Non_NA_200], Pion_200$Mean_MF_Loads[Non_NA_200], type = "l", ylim = c(0, max(Pion_200$Mean_MF_Loads[Non_NA_200])), lwd = 2)
plot(log(Pion_200$Time + 1)[Non_NA_200], Pion_200$Mean_MF_Loads[Non_NA_200], type = "l", ylim = c(0, max(Pion_200$Mean_MF_Loads[Non_NA_200])), lwd = 2)

Non_NA_150 <- !is.na(Pion_150$Mean_MF_Loads)
plot(Pion_150$Time[Non_NA_150], Pion_150$Mean_MF_Loads[Non_NA_150], type = "l", ylim = c(0, max(Pion_150$Mean_MF_Loads[Non_NA_150])), lwd = 2)
plot(log(Pion_150$Time + 1)[Non_NA_150], Pion_150$Mean_MF_Loads[Non_NA_150], type = "l", ylim = c(0, max(Pion_150$Mean_MF_Loads[Non_NA_150])), lwd = 2)

plot(Pion_Overall$Time, Pion_Overall$Mean_MF_Loads, type = "l", ylim = c(0, 15000), lwd = 2)
lines(Pion_200$Time[Non_NA_200], Pion_200$Mean_MF_Loads[Non_NA_200], type = "l", ylim = c(0, max(Pion_200$Mean_MF_Loads[Non_NA_200])), lwd = 2, col = "red")
lines(Pion_150$Time[Non_NA_150], Pion_150$Mean_MF_Loads[Non_NA_150], type = "l", ylim = c(0, max(Pion_150$Mean_MF_Loads[Non_NA_150])), lwd = 2, col = "blue")

plot(Pion_Overall$Time, Pion_Overall$Mean_MF_Loads/Pion_Overall$Mean_MF_Loads[1], type = "l", ylim = c(0, 1), lwd = 2)
lines(Pion_200$Time[Non_NA_200], Pion_200$Mean_MF_Loads[Non_NA_200]/Pion_200$Mean_MF_Loads[Non_NA_200][1], type = "l", ylim = c(0, max(Pion_200$Mean_MF_Loads[Non_NA_200])), lwd = 2, col = "red")
lines(Pion_150$Time[Non_NA_150], Pion_150$Mean_MF_Loads[Non_NA_150]/Pion_150$Mean_MF_Loads[Non_NA_150][1], type = "l", ylim = c(0, max(Pion_150$Mean_MF_Loads[Non_NA_150])), lwd = 2, col = "blue")

plot(log(Pion_Overall$Time + 1), Pion_Overall$Mean_MF_Loads, type = "l", ylim = c(0, 17500), lwd = 2)
lines(log(Pion_200$Time + 1)[Non_NA_200], Pion_200$Mean_MF_Loads[Non_NA_200], type = "l", ylim = c(0, max(Pion_200$Mean_MF_Loads[Non_NA_200])), lwd = 2, col = "red")
lines(log(Pion_150$Time + 1)[Non_NA_150], Pion_150$Mean_MF_Loads[Non_NA_150], type = "l", ylim = c(0, max(Pion_150$Mean_MF_Loads[Non_NA_150])), lwd = 2, col = "blue")

plot(log(Pion_Overall$Time + 1), Pion_Overall$Mean_MF_Loads/Pion_Overall$Mean_MF_Loads[1], type = "l", ylim = c(0, 1), lwd = 2)
lines(log(Pion_200$Time + 1)[Non_NA_200], Pion_200$Mean_MF_Loads[Non_NA_200]/Pion_200$Mean_MF_Loads[Non_NA_200][1], type = "l", ylim = c(0, max(Pion_200$Mean_MF_Loads[Non_NA_200])), lwd = 2, col = "red")
lines(log(Pion_150$Time + 1)[Non_NA_150], Pion_150$Mean_MF_Loads[Non_NA_150]/Pion_150$Mean_MF_Loads[Non_NA_150][1], type = "l", ylim = c(0, max(Pion_150$Mean_MF_Loads[Non_NA_150])), lwd = 2, col = "blue")


# Kamgno 2016 Data

## 6 doses
Kamgno_6_doses <- read.csv("C:/Users/Charlie Whittaker/Documents/13th April 2018/Loa Loa/EpiLOA Data Analysis/Parameter Estimation and Model Building/IVM & ALB Efficacy Estimation/Kamgno_2016/6_Alb_Kamgno_2016.csv")
names(Kamgno_6_doses) <- c("Time", "Mean_MF_Loads", "Number_Individuals")
weights <- Kamgno_6_doses$Number_Individuals/sum(Kamgno_6_doses$Number_Individuals)
for (i in 1:length(Kamgno_6_doses$Time)) {
  if (weights[i] < 0.05) {
    weights[i] <- 0.05
  }
}
plot(Kamgno_6_doses$Time, Kamgno_6_doses$Mean_MF_Loads, cex = 2, ylim = c(0, max(Kamgno_6_doses$Mean_MF_Loads)), pch = 20)
plot(log(Kamgno_6_doses$Time), Kamgno_6_doses$Mean_MF_Loads, cex = 2, pch = 20)


## 2 doses
Kamgno_2_doses <- read.csv("C:/Users/Charlie Whittaker/Documents/13th April 2018/Loa Loa/EpiLOA Data Analysis/Parameter Estimation and Model Building/IVM & ALB Efficacy Estimation/Kamgno_2016/2_Alb_Kamgno_2016.csv")
names(Kamgno_2_doses) <- c("Time", "Mean_MF_Loads", "Number_Individuals")
weights <- Kamgno_2_doses$Number_Individuals/sum(Kamgno_2_doses$Number_Individuals)
for (i in 1:length(Kamgno_2_doses$Time)) {
  if (weights[i] < 0.05) {
    weights[i] <- 0.05
  }
}
plot(Kamgno_2_doses$Time, Kamgno_2_doses$Mean_MF_Loads, cex = 2, ylim = c(0, max(Kamgno_2_doses$Mean_MF_Loads)), pch = 20)
plot(log(Kamgno_2_doses$Time), Kamgno_2_doses$Mean_MF_Loads, cex = 2, pch = 20, ylim = c(0, max(Kamgno_2_doses$Mean_MF_Loads)))


## Control
Kamgno_Control <- read.csv("C:/Users/Charlie Whittaker/Documents/13th April 2018/Loa Loa/EpiLOA Data Analysis/Parameter Estimation and Model Building/IVM & ALB Efficacy Estimation/Kamgno_2016/Control_Alb_Kamgno_2016.csv")
names(Kamgno_Control) <- c("Time", "Mean_MF_Loads", "Number_Individuals")
weights <- Kamgno_Control$Number_Individuals/sum(Kamgno_Control$Number_Individuals)
for (i in 1:length(Kamgno_Control$Time)) {
  if (weights[i] < 0.05) {
    weights[i] <- 0.05
  }
}
plot(Kamgno_Control$Time, Kamgno_Control$Mean_MF_Loads, cex = 2, ylim = c(0, max(Kamgno_Control$Mean_MF_Loads)), pch = 20)
plot(log(Kamgno_Control$Time), Kamgno_Control$Mean_MF_Loads, cex = 2, pch = 20, ylim = c(0, max(Kamgno_Control$Mean_MF_Loads)))

# Plotting the Raw Counts for Each Treatment Arm
plot(Kamgno_Control$Time, Kamgno_Control$Mean_MF_Loads, 
     ylim = c(0, max(Kamgno_Control$Mean_MF_Loads)), type = "l", lwd = 2)
lines(Kamgno_2_doses$Time, Kamgno_2_doses$Mean_MF_Loads, 
      ylim = c(0, max(Kamgno_2_doses$Mean_MF_Loads)), type = "l", lwd = 2, col = "red")
lines(Kamgno_6_doses$Time, Kamgno_6_doses$Mean_MF_Loads,
      ylim = c(0, max(Kamgno_6_doses$Mean_MF_Loads)), type = "l", lwd = 2, col = "purple")

# Plotting the Normalised Counts for Each Treatment Arm
plot(Kamgno_Control$Time, Kamgno_Control$Mean_MF_Loads/Kamgno_Control$Mean_MF_Loads[1], 
     ylim = c(0, 1.5), type = "l", lwd = 2)
lines(Kamgno_2_doses$Time, Kamgno_2_doses$Mean_MF_Loads/Kamgno_2_doses$Mean_MF_Loads[1], 
      ylim = c(0, 1.5), type = "l", lwd = 2, col = "red")
lines(Kamgno_6_doses$Time, Kamgno_6_doses$Mean_MF_Loads/Kamgno_6_doses$Mean_MF_Loads[1],
      ylim = c(0, 1.5), type = "l", lwd = 2, col = "purple")
