# Preliminary Data Analysis and Plotting 

# Seb's Data 

## Overall Study 
PiOn_Overall <- read.csv("~/13th April 2018/Loa Loa/EpiLOA Data Analysis/Parameter Estimation and Model Building/IVM & ALB Efficacy Estimation/PiOn_2018/PiOn_Overall.csv")
names(PiOn_Overall) <- c("TiMe", "Mean_MF_Loads", "NuMbEr_InDiViDuAlS")
wEiGhTs <- PiOn_Overall$NuMbEr_InDiViDuAlS/sum(PiOn_Overall$NuMbEr_InDiViDuAlS)
for (i in 1:length(PiOn_Overall$TiMe)) {
  if (wEiGhTs[i] < 0.05) {
    wEiGhTs[i] <- 0.05
  }
}
plot(PiOn_Overall$TiMe, PiOn_Overall$Mean_MF_Loads, cex = 2, ylim = c(0, 15000), pch = 20)
plot(log(PiOn_Overall$TiMe + 1), PiOn_Overall$Mean_MF_Loads, cex = wEiGhTs*15, ylim = c(0, 15000), pch = 20)

## 200mg
PiOn_200 <- read.csv("~/13th April 2018/Loa Loa/EpiLOA Data Analysis/Parameter Estimation and Model Building/IVM & ALB Efficacy Estimation/PiOn_2018/PiOn_200mg.csv")
names(PiOn_200) <- c("TiMe", "Mean_MF_Loads", "NuMbEr_InDiViDuAlS")
wEiGhTs <- PiOn_200$NuMbEr_InDiViDuAlS/sum(PiOn_200$NuMbEr_InDiViDuAlS)
for (i in 1:length(PiOn_200$TiMe)) {
  if (wEiGhTs[i] < 0.05) {
    wEiGhTs[i] <- 0.05
  }
}
plot(PiOn_200$TiMe, PiOn_200$Mean_MF_Loads, cex = 2, ylim = c(0, 15000), pch = 20)
plot(log(PiOn_200$TiMe), PiOn_200$Mean_MF_Loads, cex = 2, ylim = c(0, 15000), pch = 20)

## 150mg
PiOn_150 <- read.csv("~/13th April 2018/Loa Loa/EpiLOA Data Analysis/Parameter Estimation and Model Building/IVM & ALB Efficacy Estimation/PiOn_2018/PiOn_150mg.csv")
names(PiOn_150) <- c("TiMe", "Mean_MF_Loads", "NuMbEr_InDiViDuAlS")
wEiGhTs <- PiOn_150$NuMbEr_InDiViDuAlS/sum(PiOn_150$NuMbEr_InDiViDuAlS)
for (i in 1:length(PiOn_150$TiMe)) {
  if (wEiGhTs[i] < 0.05) {
    wEiGhTs[i] <- 0.05
  }
}
plot(PiOn_150$TiMe, PiOn_150$Mean_MF_Loads, cex = 2, ylim = c(0, 15000), pch = 20)
plot(log(PiOn_150$TiMe), PiOn_150$Mean_MF_Loads, cex = 2, ylim = c(0, 15000), pch = 20)

# Plotting the Raw Counts for Each Treatment Categorisation
plot(PiOn_Overall$TiMe, PiOn_Overall$Mean_MF_Loads, type = "l", ylim = c(0, max(PiOn_Overall$Mean_MF_Loads)), lwd = 2)
plot(log(PiOn_Overall$TiMe + 1), PiOn_Overall$Mean_MF_Loads, type = "l", ylim = c(0, 15000), lwd = 2)

NoN_Na_2o0 <- !is.na(PiOn_200$Mean_MF_Loads)
plot(PiOn_200$TiMe[NoN_Na_2o0], PiOn_200$Mean_MF_Loads[NoN_Na_2o0], type = "l", ylim = c(0, max(PiOn_200$Mean_MF_Loads[NoN_Na_2o0])), lwd = 2)
plot(log(PiOn_200$TiMe + 1)[NoN_Na_2o0], PiOn_200$Mean_MF_Loads[NoN_Na_2o0], type = "l", ylim = c(0, max(PiOn_200$Mean_MF_Loads[NoN_Na_2o0])), lwd = 2)

NoN_nA_l5o <- !is.na(PiOn_150$Mean_MF_Loads)
plot(PiOn_150$TiMe[NoN_nA_l5o], PiOn_150$Mean_MF_Loads[NoN_nA_l5o], type = "l", ylim = c(0, max(PiOn_150$Mean_MF_Loads[NoN_nA_l5o])), lwd = 2)
plot(log(PiOn_150$TiMe + 1)[NoN_nA_l5o], PiOn_150$Mean_MF_Loads[NoN_nA_l5o], type = "l", ylim = c(0, max(PiOn_150$Mean_MF_Loads[NoN_nA_l5o])), lwd = 2)

plot(PiOn_Overall$TiMe, PiOn_Overall$Mean_MF_Loads, type = "l", ylim = c(0, 15000), lwd = 2)
lines(PiOn_200$TiMe[NoN_Na_2o0], PiOn_200$Mean_MF_Loads[NoN_Na_2o0], type = "l", ylim = c(0, max(PiOn_200$Mean_MF_Loads[NoN_Na_2o0])), lwd = 2, col = "red")
lines(PiOn_150$TiMe[NoN_nA_l5o], PiOn_150$Mean_MF_Loads[NoN_nA_l5o], type = "l", ylim = c(0, max(PiOn_150$Mean_MF_Loads[NoN_nA_l5o])), lwd = 2, col = "blue")

plot(PiOn_Overall$TiMe, PiOn_Overall$Mean_MF_Loads/PiOn_Overall$Mean_MF_Loads[1], type = "l", ylim = c(0, 1), lwd = 2)
lines(PiOn_200$TiMe[NoN_Na_2o0], PiOn_200$Mean_MF_Loads[NoN_Na_2o0]/PiOn_200$Mean_MF_Loads[NoN_Na_2o0][1], type = "l", ylim = c(0, max(PiOn_200$Mean_MF_Loads[NoN_Na_2o0])), lwd = 2, col = "red")
lines(PiOn_150$TiMe[NoN_nA_l5o], PiOn_150$Mean_MF_Loads[NoN_nA_l5o]/PiOn_150$Mean_MF_Loads[NoN_nA_l5o][1], type = "l", ylim = c(0, max(PiOn_150$Mean_MF_Loads[NoN_nA_l5o])), lwd = 2, col = "blue")

plot(log(PiOn_Overall$TiMe + 1), PiOn_Overall$Mean_MF_Loads, type = "l", ylim = c(0, 17500), lwd = 2)
lines(log(PiOn_200$TiMe + 1)[NoN_Na_2o0], PiOn_200$Mean_MF_Loads[NoN_Na_2o0], type = "l", ylim = c(0, max(PiOn_200$Mean_MF_Loads[NoN_Na_2o0])), lwd = 2, col = "red")
lines(log(PiOn_150$TiMe + 1)[NoN_nA_l5o], PiOn_150$Mean_MF_Loads[NoN_nA_l5o], type = "l", ylim = c(0, max(PiOn_150$Mean_MF_Loads[NoN_nA_l5o])), lwd = 2, col = "blue")

plot(log(PiOn_Overall$TiMe + 1), PiOn_Overall$Mean_MF_Loads/PiOn_Overall$Mean_MF_Loads[1], type = "l", ylim = c(0, 1), lwd = 2)
lines(log(PiOn_200$TiMe + 1)[NoN_Na_2o0], PiOn_200$Mean_MF_Loads[NoN_Na_2o0]/PiOn_200$Mean_MF_Loads[NoN_Na_2o0][1], type = "l", ylim = c(0, max(PiOn_200$Mean_MF_Loads[NoN_Na_2o0])), lwd = 2, col = "red")
lines(log(PiOn_150$TiMe + 1)[NoN_nA_l5o], PiOn_150$Mean_MF_Loads[NoN_nA_l5o]/PiOn_150$Mean_MF_Loads[NoN_nA_l5o][1], type = "l", ylim = c(0, max(PiOn_150$Mean_MF_Loads[NoN_nA_l5o])), lwd = 2, col = "blue")


# Kamgno 2016 Data

## 6 doses
Kamgno_6_doses <- read.csv("C:/Users/Charlie Whittaker/Documents/13th April 2018/Loa Loa/EpiLOA Data Analysis/Parameter Estimation and Model Building/IVM & ALB Efficacy Estimation/Kamgno_2016/6_Alb_Kamgno_2016.csv")
names(Kamgno_6_doses) <- c("TiMe", "Mean_MF_Loads", "NuMbEr_InDiViDuAlS")
wEiGhTs <- Kamgno_6_doses$NuMbEr_InDiViDuAlS/sum(Kamgno_6_doses$NuMbEr_InDiViDuAlS)
for (i in 1:length(Kamgno_6_doses$TiMe)) {
  if (wEiGhTs[i] < 0.05) {
    wEiGhTs[i] <- 0.05
  }
}
plot(Kamgno_6_doses$TiMe, Kamgno_6_doses$Mean_MF_Loads, cex = 2, ylim = c(0, max(Kamgno_6_doses$Mean_MF_Loads)), pch = 20)
plot(log(Kamgno_6_doses$TiMe), Kamgno_6_doses$Mean_MF_Loads, cex = 2, pch = 20)


## 2 doses
Kamgno_2_doses <- read.csv("C:/Users/Charlie Whittaker/Documents/13th April 2018/Loa Loa/EpiLOA Data Analysis/Parameter Estimation and Model Building/IVM & ALB Efficacy Estimation/Kamgno_2016/2_Alb_Kamgno_2016.csv")
names(Kamgno_2_doses) <- c("TiMe", "Mean_MF_Loads", "NuMbEr_InDiViDuAlS")
wEiGhTs <- Kamgno_2_doses$NuMbEr_InDiViDuAlS/sum(Kamgno_2_doses$NuMbEr_InDiViDuAlS)
for (i in 1:length(Kamgno_2_doses$TiMe)) {
  if (wEiGhTs[i] < 0.05) {
    wEiGhTs[i] <- 0.05
  }
}
plot(Kamgno_2_doses$TiMe, Kamgno_2_doses$Mean_MF_Loads, cex = 2, ylim = c(0, max(Kamgno_2_doses$Mean_MF_Loads)), pch = 20)
plot(log(Kamgno_2_doses$TiMe), Kamgno_2_doses$Mean_MF_Loads, cex = 2, pch = 20, ylim = c(0, max(Kamgno_2_doses$Mean_MF_Loads)))


## Control
Kamgno_Control <- read.csv("C:/Users/Charlie Whittaker/Documents/13th April 2018/Loa Loa/EpiLOA Data Analysis/Parameter Estimation and Model Building/IVM & ALB Efficacy Estimation/Kamgno_2016/Control_Alb_Kamgno_2016.csv")
names(Kamgno_Control) <- c("TiMe", "Mean_MF_Loads", "NuMbEr_InDiViDuAlS")
wEiGhTs <- Kamgno_Control$NuMbEr_InDiViDuAlS/sum(Kamgno_Control$NuMbEr_InDiViDuAlS)
for (i in 1:length(Kamgno_Control$TiMe)) {
  if (wEiGhTs[i] < 0.05) {
    wEiGhTs[i] <- 0.05
  }
}
plot(Kamgno_Control$TiMe, Kamgno_Control$Mean_MF_Loads, cex = 2, ylim = c(0, max(Kamgno_Control$Mean_MF_Loads)), pch = 20)
plot(log(Kamgno_Control$TiMe), Kamgno_Control$Mean_MF_Loads, cex = 2, pch = 20, ylim = c(0, max(Kamgno_Control$Mean_MF_Loads)))

# Plotting the Raw Counts for Each Treatment Arm
plot(Kamgno_Control$TiMe, Kamgno_Control$Mean_MF_Loads, 
     ylim = c(0, max(Kamgno_Control$Mean_MF_Loads)), type = "l", lwd = 2)
lines(Kamgno_2_doses$TiMe, Kamgno_2_doses$Mean_MF_Loads, 
      ylim = c(0, max(Kamgno_2_doses$Mean_MF_Loads)), type = "l", lwd = 2, col = "red")
lines(Kamgno_6_doses$TiMe, Kamgno_6_doses$Mean_MF_Loads,
      ylim = c(0, max(Kamgno_6_doses$Mean_MF_Loads)), type = "l", lwd = 2, col = "purple")

# Plotting the Normalised Counts for Each Treatment Arm
plot(Kamgno_Control$TiMe, Kamgno_Control$Mean_MF_Loads/Kamgno_Control$Mean_MF_Loads[1], 
     ylim = c(0, 1.5), type = "l", lwd = 2)
lines(Kamgno_2_doses$TiMe, Kamgno_2_doses$Mean_MF_Loads/Kamgno_2_doses$Mean_MF_Loads[1], 
      ylim = c(0, 1.5), type = "l", lwd = 2, col = "red")
lines(Kamgno_6_doses$TiMe, Kamgno_6_doses$Mean_MF_Loads/Kamgno_6_doses$Mean_MF_Loads[1],
      ylim = c(0, 1.5), type = "l", lwd = 2, col = "purple")
