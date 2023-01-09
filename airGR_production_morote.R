# ////////////////////////////////////////////////////////////////////////////////////////////////////////////
# INSTITUTO TECNOLOGICO DE COSTA RICA
# Escuela de Ingenieria en Construccion
# Construction Engineering School
# https://www.tec.ac.cr
# Eng. MSc. Maikel Mendez Morales
# Email: maikel.mendez@gmail.com; mamendez@itcr.ac.cr
# https://orcid.org/0000-0003-1919-141X
# https://www.scopus.com/authid/detail.uri?authorId=51665581300
# https://scholar.google.com/citations?user=JnmSVFYAAAAJ&hl=en
# https://www.youtube.com/c/maikelmendez
# https://twitter.com/MaikelMendezM
# https://github.com/maikelonu
# Skype: maikel.mendez
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////

#-------------------------------------------------------------------------------------------------------------------
# INFO: This script is intended for the RAW WGS84-netCDF masking, resampling, dating, reprojecting, conversion and
# raster-stacking of (0.44� x 0.44�) CORDEX data from RCA4 (HadGEM2, EARTH, MIROC, CanESM2 and MPI) to CRTM05 GeotiFF
# at (25 x 25 km) for all climatic regions in Costa Rica. It also generates the ensemble-mean of the entire group.
#-------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------
# MANUSCRIPT TITLE:
# "Hydrological Response of Tropical Catchments to Climate Change as Modeled by the GR2M Model: A Case Study in Costa Rica"
#
#-------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------
# INPUT FILES:
# ASC raster.format:
# Binary R format:
#-------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------
# MANUSCRIPT FIGURES:
# None
#-------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------
# OUTPUT FILES:
# 
#-------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------
# Workspace is cleared
#rm(list = ls())

# CRAN libraries are loaded
require(airGR)
require(airGRteaching)
require(DescTools)
require(dplyr)
require(ggplot2)
require(lubridate)
require(matrixStats)
require(pastecs)
require(plyr)
require(reshape)
require(reshape2)
require(tidyr)
require(viridis)
require(weathermetrics)

# Working directory is defined
setwd("~/Documents/Github")

# Scientific notation is disabled
options(scipen=999)

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: Watershed morote
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////

#--------------------------------------------------------------------------------------------------------------
# SUBBLOCK: Observed Calibration and Validation
#--------------------------------------------------------------------------------------------------------------

# Watershed compiled files are loaded
BasinObs <- read.table("airGR_morote.txt",header=T,sep="\t",quote="")

# ET0 is converted from mm/day to mm/month
BasinObs$E <- BasinObs$E*30

# Historical DATE vectors are created according to calibration period
dates_obs <- seq(as.Date("1961-01-01"), as.Date("1990-12-01"), by = "1 month")

# Dates are transformed to POSIX class
dates_obs <- as.POSIXlt(dates_obs, format="%Y-%m-%d")

# Watershed compiled data.frame is generated
BasinObs$DatesR <- dates_obs

# data.frame is requested
View(BasinObs)

# Flowpeaks are compensated due to overestimation
factorC <- 0.70 # 0.0 to 1.00 proportional factor is created

# A temporal vector is created
basin_temp_01 <- BasinObs$Qmm*factorC

# Maximum flowpeak threshold is defined (mm/month)
basin_temp_02 <-ifelse(BasinObs$Qmm >= 300, basin_temp_01, BasinObs$Qmm)

# Compensated flows are replaced
BasinObs$Qmm <- basin_temp_02

# GR2M InputsModel object is created
InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR2M,
                                 DatesR = BasinObs$DatesR,
                                 Precip = BasinObs$P,
                                 PotEvap = BasinObs$E)

#--------------------------------------------------------------------------------------------------------------
# SUBBLOCK: Historical Modelling
#--------------------------------------------------------------------------------------------------------------

# CALIBRATION ++++++++++

# GR2M run period is selected
Ind_Run <- seq(which(format(BasinObs$DatesR, format = "%Y-%m")=="1970-05"),
               which(format(BasinObs$DatesR, format = "%Y-%m")=="1975-05"))

# GR2M RunOptions object is created
RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR2M,
                               InputsModel = InputsModel,
                               IndPeriod_Run = Ind_Run)

# Optimum parameters are defined
Param <- c(X1 =  772.7843, X2 = 0.9675)

# GR2M OutputsModel object is created
OutputsModel <- RunModel_GR2M(InputsModel = InputsModel,
                              RunOptions = RunOptions,
                              Param = Param)

# A GR2M summary plot is requested
plot(OutputsModel, Qobs = BasinObs$Qmm[Ind_Run],
            cex.axis =0.9, cex.lab = 0.9, cex.leg = 0.9,
            BasinArea=929.40, which = c("all"))

# Nash-Sutcliffe Efficiency is defined
InputsCrit  <- CreateInputsCrit(FUN_CRIT = ErrorCrit_NSE,
                                InputsModel = InputsModel,
                                RunOptions = RunOptions,
                                Obs = BasinObs$Qmm[Ind_Run])

# Nash-Sutcliffe Efficiency is requested
OutputsCrit <- ErrorCrit_NSE(InputsCrit = InputsCrit,
                             OutputsModel = OutputsModel)


# KGE Efficiency is defined
InputsCrit02  <- CreateInputsCrit(FUN_CRIT = ErrorCrit_KGE,
                                  InputsModel = InputsModel,
                                  RunOptions = RunOptions,
                                  Obs = BasinObs$Qmm[Ind_Run])

# KGE Efficiency is requested
OutputsCrit02 <- ErrorCrit_KGE(InputsCrit = InputsCrit02,
                               OutputsModel = OutputsModel)

# VALIDATION ++++++++++

# GR2M run period is selected
Ind_Run02 <- seq(which(format(BasinObs$DatesR, format = "%Y-%m")=="1975-06"),
                 which(format(BasinObs$DatesR, format = "%Y-%m")=="1981-05"))

# GR2M RunOptions02 object is created
RunOptions02 <- CreateRunOptions(FUN_MOD = RunModel_GR2M,
                                 InputsModel = InputsModel,
                                 IndPeriod_Run = Ind_Run02)

# GR2M OutputsModel02 object is created
OutputsModel02 <- RunModel_GR2M(InputsModel = InputsModel,
                                RunOptions = RunOptions02,
                                Param = Param)

# A GR2M summary plot is requested
plot(OutputsModel02, Qobs = BasinObs$Qmm[Ind_Run02],
     cex.axis =0.9, cex.lab = 0.9, cex.leg = 0.9,
     BasinArea=929.40, which = c("all"))

# Nash-Sutcliffe Efficiency is defined
InputsCritVal  <- CreateInputsCrit(FUN_CRIT = ErrorCrit_NSE,
                                   InputsModel = InputsModel,
                                   RunOptions = RunOptions02,
                                   Obs = BasinObs$Qmm[Ind_Run02])

# Nash-Sutcliffe Efficiency is requested
OutputsCritVal <- ErrorCrit_NSE(InputsCrit = InputsCritVal,
                                OutputsModel = OutputsModel02)


# KGE Efficiency is defined
InputsCritVal02  <- CreateInputsCrit(FUN_CRIT = ErrorCrit_KGE,
                                     InputsModel = InputsModel,
                                     RunOptions = RunOptions02,
                                     Obs = BasinObs$Qmm[Ind_Run02])

# KGE Efficiency is requested
OutputsCritVal02 <- ErrorCrit_KGE(InputsCrit = InputsCritVal02,
                                  OutputsModel = OutputsModel02)

#--------------------------------------------------------------------------------------------------------------
# SUBBLOCK: Parameter Optimization
#--------------------------------------------------------------------------------------------------------------

# GR2M CalibOptions object is created
CalibOptions <- CreateCalibOptions(FUN_MOD = RunModel_GR2M,
                                   FUN_CALIB = Calibration_Michel)

# GR2M OutputsCalib object is created
OutputsCalib <- Calibration(InputsModel = InputsModel,
                            RunOptions = RunOptions,
                            InputsCrit = InputsCrit, # Observations are included here
                            CalibOptions = CalibOptions,
                            FUN_MOD = RunModel_GR2M,
                            FUN_CALIB = Calibration_Michel)

# OutputsCalib outputs are requested
OutputsCalib$ParamFinalR
OutputsCalib$CritFinal

# A final GR2M simulation is executed
Param <- OutputsCalib$ParamFinalR
OutputsModel <- RunModel(InputsModel = InputsModel, RunOptions = RunOptions,
                         Param = Param, FUN = RunModel_GR2M)

# A GR2M summary plot is requested
plot(OutputsModel, Qobs = BasinObs$Qmm[Ind_Run])

# A simple mass balance is requested
sum(OutputsModel$Qsim)
sum(BasinObs$Qmm[Ind_Run])
sum(OutputsModel$Qsim)/sum(BasinObs$Qmm[Ind_Run])

#--------------------------------------------------------------------------------------------------------------
# SUBBLOCK: 1961-1990 Historical Modelling
#--------------------------------------------------------------------------------------------------------------

# Watershed compiled files are loaded
BasinObs_1961_1990 <- read.table("airGR_morote_1961_1990.txt",header=T,sep="\t",quote="")

# ET0 is converted from mm/day to mm/month
BasinObs_1961_1990$E <- BasinObs_1961_1990$E*30

# Historical date vectors are created
dates_obs_1961_1990 <- seq(as.Date("1961-01-01"), as.Date("1990-12-01"), by = "1 month")

# Dates are transformed to POSIX class
dates_obs_1961_1990 <- as.POSIXlt(dates_obs_1961_1990, format="%Y-%m-%d")

# Watershed compiled data.frame is generated
BasinObs_1961_1990$DatesR <- dates_obs_1961_1990

# GR2M InputsModel object is created
InputsModel_1961_1990 <- CreateInputsModel(FUN_MOD = RunModel_GR2M,
                                           DatesR = BasinObs_1961_1990$DatesR,
                                           Precip = BasinObs_1961_1990$P,
                                           PotEvap = BasinObs_1961_1990$E)

# GR2M run period is selected
Ind_Run_1961_1990 <- seq(which(format(BasinObs_1961_1990$DatesR, format = "%Y-%m")=="1961-01"),
                         which(format(BasinObs_1961_1990$DatesR, format = "%Y-%m")=="1990-12"))

# GR2M RunOptions object is created
RunOptions_1961_1990 <- CreateRunOptions(FUN_MOD = RunModel_GR2M,
                                         InputsModel = InputsModel_1961_1990,
                                         IndPeriod_Run = Ind_Run_1961_1990)

# GR2M OutputsModel object is created
OutputsModel_1961_1990 <- RunModel_GR2M(InputsModel = InputsModel_1961_1990,
                                        RunOptions = RunOptions_1961_1990,
                                        Param = Param)

# A GR2M summary plot is requested
plot(OutputsModel_1961_1990, Qobs = BasinObs_1961_1990$Qmm[Ind_Run_1961_1990],
            cex.axis =0.9, cex.lab = 0.9, cex.leg = 0.9,
            BasinArea=929.40, which = c("all"))

# Mass balance of the period is calculated
mass_1961_1990 <- sum(OutputsModel_1961_1990$Qsim)
AE_1961_1990 <- sum(OutputsModel_1961_1990$AE)

# Mass balance of the period is requested
mass_1961_1990

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: Future Simulation RCP85
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////

#--------------------------------------------------------------------------------------------------------------
# SUBBLOCK: RCP85 2011-2040 Future Simulation RCP85
#--------------------------------------------------------------------------------------------------------------

# Watershed compiled files are loaded
BasinObs_2011_2040 <- read.table("airGR_morote_2011_2040_RCP85.txt",header=T,sep="\t",quote="")

# ET0 is converted from mm/day to mm/month
BasinObs_2011_2040$E <- BasinObs_2011_2040$E*30

# Historical date vectors are created
dates_obs_2011_2040 <- seq(as.Date("2011-01-01"), as.Date("2040-12-01"), by = "1 month")

# Dates are transformed to POSIX class
dates_obs_2011_2040 <- as.POSIXlt(dates_obs_2011_2040, format="%Y-%m-%d")

# Watershed compiled data.frame is generated
BasinObs_2011_2040$DatesR <- dates_obs_2011_2040

# GR2M InputsModel object is created
InputsModel_2011_2040 <- CreateInputsModel(FUN_MOD = RunModel_GR2M,
                                           DatesR = BasinObs_2011_2040$DatesR,
                                           Precip = BasinObs_2011_2040$P,
                                           PotEvap = BasinObs_2011_2040$E)

# GR2M run period is selected
Ind_Run_2011_2040 <- seq(which(format(BasinObs_2011_2040$DatesR, format = "%Y-%m")=="2011-01"),
                         which(format(BasinObs_2011_2040$DatesR, format = "%Y-%m")=="2040-12"))

# GR2M RunOptions object is created
RunOptions_2011_2040 <- CreateRunOptions(FUN_MOD = RunModel_GR2M,
                                         InputsModel = InputsModel_2011_2040,
                                         IndPeriod_Run = Ind_Run_2011_2040)

# GR2M OutputsModel object is created
OutputsModel_2011_2040 <- RunModel_GR2M(InputsModel = InputsModel_2011_2040,
                                        RunOptions = RunOptions_2011_2040,
                                        Param = Param)

# A GR2M summary plot is requested
plot(OutputsModel_2011_2040, Qobs = BasinObs_2011_2040$Qmm[Ind_Run_2011_2040],
            cex.axis =0.9, cex.lab = 0.9, cex.leg = 0.9,
            BasinArea=929.40, which = c("all"))

# Mass balance of the period is calculated
mass_2011_2040 <- sum(OutputsModel_2011_2040$Qsim)
AE_2011_2040 <- sum(OutputsModel_2011_2040$AE)

# Mass balance of the period is requested
mass_2011_2040

#--------------------------------------------------------------------------------------------------------------
# SUBBLOCK: RCP85 2041_2070 Future Simulation RCP85
#--------------------------------------------------------------------------------------------------------------

# Watershed compiled files are loaded
BasinObs_2041_2070 <- read.table("airGR_morote_2041_2070_RCP85.txt",header=T,sep="\t",quote="")

# ET0 is converted from mm/day to mm/month
BasinObs_2041_2070$E <- BasinObs_2041_2070$E*30

# Historical date vectors are created
dates_obs_2041_2070 <- seq(as.Date("2041-01-01"), as.Date("2070-12-01"), by = "1 month")

# Dates are transformed to POSIX class
dates_obs_2041_2070 <- as.POSIXlt(dates_obs_2041_2070, format="%Y-%m-%d")

# Watershed compiled data.frame is generated
BasinObs_2041_2070$DatesR <- dates_obs_2041_2070

# GR2M InputsModel object is created
InputsModel_2041_2070 <- CreateInputsModel(FUN_MOD = RunModel_GR2M,
                                           DatesR = BasinObs_2041_2070$DatesR,
                                           Precip = BasinObs_2041_2070$P,
                                           PotEvap = BasinObs_2041_2070$E)

# GR2M run period is selected
Ind_Run_2041_2070 <- seq(which(format(BasinObs_2041_2070$DatesR, format = "%Y-%m")=="2041-01"),
                         which(format(BasinObs_2041_2070$DatesR, format = "%Y-%m")=="2070-12"))

# GR2M RunOptions object is created
RunOptions_2041_2070 <- CreateRunOptions(FUN_MOD = RunModel_GR2M,
                                         InputsModel = InputsModel_2041_2070,
                                         IndPeriod_Run = Ind_Run_2041_2070)

# GR2M OutputsModel object is created
OutputsModel_2041_2070 <- RunModel_GR2M(InputsModel = InputsModel_2041_2070,
                                        RunOptions = RunOptions_2041_2070,
                                        Param = Param)

# A GR2M summary plot is requested
plot(OutputsModel_2041_2070, Qobs = BasinObs_2041_2070$Qmm[Ind_Run_2041_2070],
            cex.axis =0.9, cex.lab = 0.9, cex.leg = 0.9,
            BasinArea=929.40, which = c("all"))

# Mass balance of the period is calculated
mass_2041_2070 <- sum(OutputsModel_2041_2070$Qsim)
AE_2041_2070 <- sum(OutputsModel_2041_2070$AE)

# Mass balance of the period is requested
mass_2041_2070

#--------------------------------------------------------------------------------------------------------------
# SUBBLOCK: RCP85 2071_2100 Future Simulation RCP85
#--------------------------------------------------------------------------------------------------------------

# Watershed compiled files are loaded
BasinObs_2071_2100 <- read.table("airGR_morote_2071_2100_RCP85.txt",header=T,sep="\t",quote="")

# ET0 is converted from mm/day to mm/month
BasinObs_2071_2100$E <- BasinObs_2071_2100$E*30

# Historical date vectors are created
dates_obs_2071_2100 <- seq(as.Date("2071-01-01"), as.Date("2100-12-01"), by = "1 month")

# Dates are transformed to POSIX class
dates_obs_2071_2100 <- as.POSIXlt(dates_obs_2071_2100, format="%Y-%m-%d")

# Watershed compiled data.frame is generated
BasinObs_2071_2100$DatesR <- dates_obs_2071_2100

# GR2M InputsModel object is created
InputsModel_2071_2100 <- CreateInputsModel(FUN_MOD = RunModel_GR2M,
                                           DatesR = BasinObs_2071_2100$DatesR,
                                           Precip = BasinObs_2071_2100$P,
                                           PotEvap = BasinObs_2071_2100$E)

# GR2M run period is selected
Ind_Run_2071_2100 <- seq(which(format(BasinObs_2071_2100$DatesR, format = "%Y-%m")=="2071-01"),
                         which(format(BasinObs_2071_2100$DatesR, format = "%Y-%m")=="2100-12"))

# GR2M RunOptions object is created
RunOptions_2071_2100 <- CreateRunOptions(FUN_MOD = RunModel_GR2M,
                                         InputsModel = InputsModel_2071_2100,
                                         IndPeriod_Run = Ind_Run_2071_2100)

# GR2M OutputsModel object is created
OutputsModel_2071_2100 <- RunModel_GR2M(InputsModel = InputsModel_2071_2100,
                                        RunOptions = RunOptions_2071_2100,
                                        Param = Param)

# A GR2M summary plot is requested
plot(OutputsModel_2071_2100, Qobs = BasinObs_2071_2100$Qmm[Ind_Run_2071_2100],
            cex.axis =0.9, cex.lab = 0.9, cex.leg = 0.9,
            BasinArea=929.40, which = c("all"))

# Mass balance of the period is calculated
mass_2071_2100 <- sum(OutputsModel_2071_2100$Qsim)
AE_2071_2100 <- sum(OutputsModel_2071_2100$AE)

# Mass balance of the period is requested
mass_2071_2100

# 2011_2040 PBIAS is calculated with respect to 1961-1990
PBIAS_2011_2040_RCP85 <- ((mass_2011_2040 - mass_1961_1990) / mass_1961_1990)*100
PBIAS_2011_2040_RCP85

# 2041_2070 PBIAS is calculated with respect to 1961-1990
PBIAS_2041_2070_RCP85 <- ((mass_2041_2070 - mass_1961_1990) / mass_1961_1990)*100
PBIAS_2041_2070_RCP85

# 2041_2070 PBIAS is calculated with respect to 1961-1990
PBIAS_2071_2100_RCP85 <- ((mass_2071_2100 - mass_1961_1990) / mass_1961_1990)*100
PBIAS_2071_2100_RCP85

# 2011_2040 AE_PBIAS is calculated with respect to 1961-1990
AE_PBIAS_2011_2040_RCP85 <- ((AE_2011_2040 - AE_1961_1990) / AE_1961_1990)*100
AE_PBIAS_2011_2040_RCP85

# 2041_2070 AE_PBIAS is calculated with respect to 1961-1990
AE_PBIAS_2041_2070_RCP85 <- ((AE_2041_2070 - AE_1961_1990) / AE_1961_1990)*100
AE_PBIAS_2041_2070_RCP85

# 2041_2070 AE_PBIAS is calculated with respect to 1961-1990
AE_PBIAS_2071_2100_RCP85 <- ((AE_2071_2100 - AE_1961_1990) / AE_1961_1990)*100
AE_PBIAS_2071_2100_RCP85

#--------------------------------------------------------------------------------------------------------------
# SUBBLOCK: Data Processing RCP85
#--------------------------------------------------------------------------------------------------------------

# An optimization data.frame is created
df.opti_rcp85 <- BasinObs

# lubridate Library functions are applied to df_WCD to create new columns contaning:
# YEAR, YEAR_CH, MONTH, MONTH_CH
df.opti_rcp85$YEAR <- year(df.opti_rcp85$DatesR) # Years component of a date-time
df.opti_rcp85$YEAR_CH <- as.character(year(df.opti_rcp85$DatesR)) # Years component of a date-time as character
df.opti_rcp85$MONTH <- month(df.opti_rcp85$DatesR, label = FALSE) # Months component of a date-time
df.opti_rcp85$MONTH_CH <- month(df.opti_rcp85$DatesR, label = TRUE) #

# Relevant variables are selected
df.opti_rcp85 <- df.opti_rcp85[, c(1,9)]

# NA values are deleted
df.opti_rcp85 <-  df.opti_rcp85[complete.cases(df.opti_rcp85), ]

# A dummy variable is created
df.opti_rcp85$value <- c("value")

# A new data.frame is created
df.opti_rcp85.plot <- as.data.frame(dcast(df.opti_rcp85, MONTH_CH  ~ value, mean, value.var = "Qmm"))

# A dummy variable is created
df.opti_rcp85.plot$class <- c("obs")

# **********************************************
# An modelling data.frame is created
df.modelling_rcp85 <- data.frame(c(OutputsModel$Qsim))
df.modelling_rcp85$Dates <- OutputsModel$DatesR

# data.frame names are replaced
names(df.modelling_rcp85) <- c("Qsim", "DatesR")

# lubridate Library functions are applied to df_WCD to create new columns contaning:
# YEAR, YEAR_CH, MONTH, MONTH_CH
df.modelling_rcp85$YEAR <- year(df.modelling_rcp85$DatesR) # Years component of a date-time
df.modelling_rcp85$YEAR_CH <- as.character(year(df.modelling_rcp85$DatesR)) # Years component of a date-time as character
df.modelling_rcp85$MONTH <- month(df.modelling_rcp85$DatesR, label = FALSE) # Months component of a date-time
df.modelling_rcp85$MONTH_CH <- month(df.modelling_rcp85$DatesR, label = TRUE) #

# Relevant variables are selected
df.modelling_rcp85 <- df.modelling_rcp85[, c(1,6)]

# A dummy variable is created
df.modelling_rcp85$value <- c("value")

# A new data.frame is created
df.modelling_rcp85 <- as.data.frame(dcast(df.modelling_rcp85, MONTH_CH  ~ value, mean, value.var = "Qsim"))

# A dummy variable is created
df.modelling_rcp85$class <- c("opti")

# A production data.frame is created
df.modelling_rcp85.cal <- df.modelling_rcp85

# **********************************************
# An simulation 1961_1990 data.frame is created
df.modelling_rcp85 <- data.frame(c(OutputsModel_1961_1990$Qsim))
df.modelling_rcp85$Dates <- OutputsModel_1961_1990$DatesR

# data.frame names are replaced
names(df.modelling_rcp85) <- c("Qsim", "DatesR")

# lubridate Library functions are applied to df_WCD to create new columns contaning:
# YEAR, YEAR_CH, MONTH, MONTH_CH
df.modelling_rcp85$YEAR <- year(df.modelling_rcp85$DatesR) # Years component of a date-time
df.modelling_rcp85$YEAR_CH <- as.character(year(df.modelling_rcp85$DatesR)) # Years component of a date-time as character
df.modelling_rcp85$MONTH <- month(df.modelling_rcp85$DatesR, label = FALSE) # Months component of a date-time
df.modelling_rcp85$MONTH_CH <- month(df.modelling_rcp85$DatesR, label = TRUE) #

# Relevant variables are selected
df.modelling_rcp85 <- df.modelling_rcp85[, c(1,6)]

# A dummy variable is created
df.modelling_rcp85$value <- c("value")

# A new data.frame is created
df.modelling_rcp85 <- as.data.frame(dcast(df.modelling_rcp85, MONTH_CH  ~ value, mean, value.var = "Qsim"))

# A dummy variable is created
df.modelling_rcp85$class <- c("1961_1990")

# A production data.frame is created
df.sim_1961_1990_rcp85 <- df.modelling_rcp85

# **********************************************
# An simulation 2011_2040 data.frame is created
df.modelling_rcp85 <- data.frame(c(OutputsModel_2011_2040$Qsim))
df.modelling_rcp85$Dates <- OutputsModel_2011_2040$DatesR

# data.frame names are replaced
names(df.modelling_rcp85) <- c("Qsim", "DatesR")

# lubridate Library functions are applied to df_WCD to create new columns contaning:
# YEAR, YEAR_CH, MONTH, MONTH_CH
df.modelling_rcp85$YEAR <- year(df.modelling_rcp85$DatesR) # Years component of a date-time
df.modelling_rcp85$YEAR_CH <- as.character(year(df.modelling_rcp85$DatesR)) # Years component of a date-time as character
df.modelling_rcp85$MONTH <- month(df.modelling_rcp85$DatesR, label = FALSE) # Months component of a date-time
df.modelling_rcp85$MONTH_CH <- month(df.modelling_rcp85$DatesR, label = TRUE) #

# Relevant variables are selected
df.modelling_rcp85 <- df.modelling_rcp85[, c(1,6)]

# A dummy variable is created
df.modelling_rcp85$value <- c("value")

# A new data.frame is created
df.modelling_rcp85 <- as.data.frame(dcast(df.modelling_rcp85, MONTH_CH  ~ value, mean, value.var = "Qsim"))

# A dummy variable is created
df.modelling_rcp85$class <- c("2011_2040")

# A production data.frame is created
df.sim_2011_2040_rcp85 <- df.modelling_rcp85

# **********************************************
# An simulation 2041_2070 data.frame is created
df.modelling_rcp85 <- data.frame(c(OutputsModel_2041_2070$Qsim))
df.modelling_rcp85$Dates <- OutputsModel_2041_2070$DatesR

# data.frame names are replaced
names(df.modelling_rcp85) <- c("Qsim", "DatesR")

# lubridate Library functions are applied to df_WCD to create new columns contaning:
# YEAR, YEAR_CH, MONTH, MONTH_CH
df.modelling_rcp85$YEAR <- year(df.modelling_rcp85$DatesR) # Years component of a date-time
df.modelling_rcp85$YEAR_CH <- as.character(year(df.modelling_rcp85$DatesR)) # Years component of a date-time as character
df.modelling_rcp85$MONTH <- month(df.modelling_rcp85$DatesR, label = FALSE) # Months component of a date-time
df.modelling_rcp85$MONTH_CH <- month(df.modelling_rcp85$DatesR, label = TRUE) #

# Relevant variables are selected
df.modelling_rcp85 <- df.modelling_rcp85[, c(1,6)]

# A dummy variable is created
df.modelling_rcp85$value <- c("value")

# A new data.frame is created
df.modelling_rcp85 <- as.data.frame(dcast(df.modelling_rcp85, MONTH_CH  ~ value, mean, value.var = "Qsim"))

# A dummy variable is created
df.modelling_rcp85$class <- c("2041_2070")

# A production data.frame is created
df.sim_2041_2070_rcp85 <- df.modelling_rcp85

# **********************************************
# An simulation 2071_2100 data.frame is created
df.modelling_rcp85 <- data.frame(c(OutputsModel_2071_2100$Qsim))
df.modelling_rcp85$Dates <- OutputsModel_2071_2100$DatesR

# data.frame names are replaced
names(df.modelling_rcp85) <- c("Qsim", "DatesR")

# lubridate Library functions are applied to df_WCD to create new columns contaning:
# YEAR, YEAR_CH, MONTH, MONTH_CH
df.modelling_rcp85$YEAR <- year(df.modelling_rcp85$DatesR) # Years component of a date-time
df.modelling_rcp85$YEAR_CH <- as.character(year(df.modelling_rcp85$DatesR)) # Years component of a date-time as character
df.modelling_rcp85$MONTH <- month(df.modelling_rcp85$DatesR, label = FALSE) # Months component of a date-time
df.modelling_rcp85$MONTH_CH <- month(df.modelling_rcp85$DatesR, label = TRUE) #

# Relevant variables are selected
df.modelling_rcp85 <- df.modelling_rcp85[, c(1,6)]

# A dummy variable is created
df.modelling_rcp85$value <- c("value")

# A new data.frame is created
df.modelling_rcp85 <- as.data.frame(dcast(df.modelling_rcp85, MONTH_CH  ~ value, mean, value.var = "Qsim"))

# A dummy variable is created
df.modelling_rcp85$class <- c("2071_2100")

# A production data.frame is created
df.sim_2071_2100_rcp85 <- df.modelling_rcp85

# **********************************************
# An export data.frame is created
df.export_rcp85 <- rbind(df.sim_1961_1990_rcp85, df.sim_2011_2040_rcp85, df.sim_2041_2070_rcp85, df.sim_2071_2100_rcp85)
#df.export_rcp85 <- rbind(df.opti_rcp85.plot, df.modelling_rcp85.cal, df.sim_1961_1990_rcp85, df.sim_2011_2040_rcp85, df.sim_2041_2070_rcp85, df.sim_2071_2100_rcp85)

# A repetition variable is created
df.export_rcp85$REP <- rep(df.sim_1961_1990_rcp85$value, 4)

# A ratio variable is created
df.export_rcp85$RATIO <- ((df.export_rcp85$value / df.export_rcp85$REP) -1)*100

# Irrelevant variables are deleted
df.export_rcp85 <- df.export_rcp85[, c(-4)]

# data.frame is converted to long format
df.export_rcp85 <- melt(df.export_rcp85)

#==================================
# Mass-Balance-Components RCP85
#==================================

# $PotEvap [numeric] series of input potential evapotranspiration [mm/month]
# $Precip  [numeric] series of input total precipitation [mm/month]
# $AE	     [numeric] series of actual evapotranspiration [mm/month]
# $Perc    [numeric] series of percolation (P2) [mm/month]
# $Qsim    [numeric] series of simulated discharge [mm/month]

# Historical
vector.hist <- c((sum(OutputsModel_1961_1990$Precip))/30,
                 (sum(OutputsModel_1961_1990$AE))/30,
                 (sum(OutputsModel_1961_1990$Perc))/30,
                 (sum(OutputsModel_1961_1990$Qsim))/30)

# RCP85 2011_2040
vector.rcp85.2011_2040 <- c((sum(OutputsModel_2011_2040$Precip))/30,
                            (sum(OutputsModel_2011_2040$AE))/30,
                            (sum(OutputsModel_2011_2040$Perc))/30,
                            (sum(OutputsModel_2011_2040$Qsim))/30)

# RCP85 2041_2070
vector.rcp85.2041_2070 <- c((sum(OutputsModel_2041_2070$Precip))/30,
                            (sum(OutputsModel_2041_2070$AE))/30,
                            (sum(OutputsModel_2041_2070$Perc))/30,
                            (sum(OutputsModel_2041_2070$Qsim))/30)

# RCP85 2071_2100
vector.rcp85.2071_2100 <- c((sum(OutputsModel_2071_2100$Precip))/30,
                            (sum(OutputsModel_2071_2100$AE))/30,
                            (sum(OutputsModel_2071_2100$Perc))/30,
                            (sum(OutputsModel_2071_2100$Qsim))/30)

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: Future Simulation RCP45
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////

#--------------------------------------------------------------------------------------------------------------
# SUBBLOCK: RCP45 2011-2040 Future Simulation RCP45
#--------------------------------------------------------------------------------------------------------------

# Watershed compiled files are loaded
BasinObs_2011_2040 <- read.table("airGR_morote_2011_2040_RCP45.txt",header=T,sep="\t",quote="")

# ET0 is converted from mm/day to mm/month
BasinObs_2011_2040$E <- BasinObs_2011_2040$E*30

# Historical date vectors are created
dates_obs_2011_2040 <- seq(as.Date("2011-01-01"), as.Date("2040-12-01"), by = "1 month")

# Dates are transformed to POSIX class
dates_obs_2011_2040 <- as.POSIXlt(dates_obs_2011_2040, format="%Y-%m-%d")

# Watershed compiled data.frame is generated
BasinObs_2011_2040$DatesR <- dates_obs_2011_2040

# GR2M InputsModel object is created
InputsModel_2011_2040 <- CreateInputsModel(FUN_MOD = RunModel_GR2M,
                                           DatesR = BasinObs_2011_2040$DatesR,
                                           Precip = BasinObs_2011_2040$P,
                                           PotEvap = BasinObs_2011_2040$E)

# GR2M run period is selected
Ind_Run_2011_2040 <- seq(which(format(BasinObs_2011_2040$DatesR, format = "%Y-%m")=="2011-01"),
                         which(format(BasinObs_2011_2040$DatesR, format = "%Y-%m")=="2040-12"))

# GR2M RunOptions object is created
RunOptions_2011_2040 <- CreateRunOptions(FUN_MOD = RunModel_GR2M,
                                         InputsModel = InputsModel_2011_2040,
                                         IndPeriod_Run = Ind_Run_2011_2040)

# GR2M OutputsModel object is created
OutputsModel_2011_2040 <- RunModel_GR2M(InputsModel = InputsModel_2011_2040,
                                        RunOptions = RunOptions_2011_2040,
                                        Param = Param)

# A GR2M summary plot is requested
plot(OutputsModel_2011_2040, Qobs = BasinObs_2011_2040$Qmm[Ind_Run_2011_2040],
     cex.axis =0.9, cex.lab = 0.9, cex.leg = 0.9,
     BasinArea=929.40, which = c("all"))

# Mass balance of the period is calculated
mass_2011_2040 <- sum(OutputsModel_2011_2040$Qsim)
AE_2011_2040 <- sum(OutputsModel_2011_2040$AE)

# Mass balance of the period is requested
mass_2011_2040

#--------------------------------------------------------------------------------------------------------------
# SUBBLOCK: RCP45 2041_2070 Future Simulation RCP45
#--------------------------------------------------------------------------------------------------------------

# Watershed compiled files are loaded
BasinObs_2041_2070 <- read.table("airGR_morote_2041_2070_RCP45.txt",header=T,sep="\t",quote="")

# ET0 is converted from mm/day to mm/month
BasinObs_2041_2070$E <- BasinObs_2041_2070$E*30

# Historical date vectors are created
dates_obs_2041_2070 <- seq(as.Date("2041-01-01"), as.Date("2070-12-01"), by = "1 month")

# Dates are transformed to POSIX class
dates_obs_2041_2070 <- as.POSIXlt(dates_obs_2041_2070, format="%Y-%m-%d")

# Watershed compiled data.frame is generated
BasinObs_2041_2070$DatesR <- dates_obs_2041_2070

# GR2M InputsModel object is created
InputsModel_2041_2070 <- CreateInputsModel(FUN_MOD = RunModel_GR2M,
                                           DatesR = BasinObs_2041_2070$DatesR,
                                           Precip = BasinObs_2041_2070$P,
                                           PotEvap = BasinObs_2041_2070$E)

# GR2M run period is selected
Ind_Run_2041_2070 <- seq(which(format(BasinObs_2041_2070$DatesR, format = "%Y-%m")=="2041-01"),
                         which(format(BasinObs_2041_2070$DatesR, format = "%Y-%m")=="2070-12"))

# GR2M RunOptions object is created
RunOptions_2041_2070 <- CreateRunOptions(FUN_MOD = RunModel_GR2M,
                                         InputsModel = InputsModel_2041_2070,
                                         IndPeriod_Run = Ind_Run_2041_2070)

# GR2M OutputsModel object is created
OutputsModel_2041_2070 <- RunModel_GR2M(InputsModel = InputsModel_2041_2070,
                                        RunOptions = RunOptions_2041_2070,
                                        Param = Param)

# A GR2M summary plot is requested
plot(OutputsModel_2041_2070, Qobs = BasinObs_2041_2070$Qmm[Ind_Run_2041_2070],
     cex.axis =0.9, cex.lab = 0.9, cex.leg = 0.9,
     BasinArea=929.40, which = c("all"))

# Mass balance of the period is calculated
mass_2041_2070 <- sum(OutputsModel_2041_2070$Qsim)
AE_2041_2070 <- sum(OutputsModel_2041_2070$AE)

# Mass balance of the period is requested
mass_2041_2070

#--------------------------------------------------------------------------------------------------------------
# SUBBLOCK: RCP45 2071_2100 Future Simulation RCP45
#--------------------------------------------------------------------------------------------------------------

# Watershed compiled files are loaded
BasinObs_2071_2100 <- read.table("airGR_morote_2071_2100_RCP45.txt",header=T,sep="\t",quote="")

# ET0 is converted from mm/day to mm/month
BasinObs_2071_2100$E <- BasinObs_2071_2100$E*30

# Historical date vectors are created
dates_obs_2071_2100 <- seq(as.Date("2071-01-01"), as.Date("2100-12-01"), by = "1 month")

# Dates are transformed to POSIX class
dates_obs_2071_2100 <- as.POSIXlt(dates_obs_2071_2100, format="%Y-%m-%d")

# Watershed compiled data.frame is generated
BasinObs_2071_2100$DatesR <- dates_obs_2071_2100

# GR2M InputsModel object is created
InputsModel_2071_2100 <- CreateInputsModel(FUN_MOD = RunModel_GR2M,
                                           DatesR = BasinObs_2071_2100$DatesR,
                                           Precip = BasinObs_2071_2100$P,
                                           PotEvap = BasinObs_2071_2100$E)

# GR2M run period is selected
Ind_Run_2071_2100 <- seq(which(format(BasinObs_2071_2100$DatesR, format = "%Y-%m")=="2071-01"),
                         which(format(BasinObs_2071_2100$DatesR, format = "%Y-%m")=="2100-12"))

# GR2M RunOptions object is created
RunOptions_2071_2100 <- CreateRunOptions(FUN_MOD = RunModel_GR2M,
                                         InputsModel = InputsModel_2071_2100,
                                         IndPeriod_Run = Ind_Run_2071_2100)

# GR2M OutputsModel object is created
OutputsModel_2071_2100 <- RunModel_GR2M(InputsModel = InputsModel_2071_2100,
                                        RunOptions = RunOptions_2071_2100,
                                        Param = Param)

# A GR2M summary plot is requested
plot(OutputsModel_2071_2100, Qobs = BasinObs_2071_2100$Qmm[Ind_Run_2071_2100],
     cex.axis =0.9, cex.lab = 0.9, cex.leg = 0.9,
     BasinArea=929.40, which = c("all"))

# Mass balance of the period is calculated
mass_2071_2100 <- sum(OutputsModel_2071_2100$Qsim)
AE_2071_2100 <- sum(OutputsModel_2071_2100$AE)

# Mass balance of the period is requested
mass_2071_2100

# 2011_2040 PBIAS is calculated with respect to 1961-1990
PBIAS_2011_2040_RCP45 <- ((mass_2011_2040 - mass_1961_1990) / mass_1961_1990)*100
PBIAS_2011_2040_RCP45

# 2041_2070 PBIAS is calculated with respect to 1961-1990
PBIAS_2041_2070_RCP45 <- ((mass_2041_2070 - mass_1961_1990) / mass_1961_1990)*100
PBIAS_2041_2070_RCP45

# 2041_2070 PBIAS is calculated with respect to 1961-1990
PBIAS_2071_2100_RCP45 <- ((mass_2071_2100 - mass_1961_1990) / mass_1961_1990)*100
PBIAS_2071_2100_RCP45

# 2011_2040 AE_PBIAS is calculated with respect to 1961-1990
AE_PBIAS_2011_2040_RCP45 <- ((AE_2011_2040 - AE_1961_1990) / AE_1961_1990)*100
AE_PBIAS_2011_2040_RCP45

# 2041_2070 AE_PBIAS is calculated with respect to 1961-1990
AE_PBIAS_2041_2070_RCP45 <- ((AE_2041_2070 - AE_1961_1990) / AE_1961_1990)*100
AE_PBIAS_2041_2070_RCP45

# 2041_2070 AE_PBIAS is calculated with respect to 1961-1990
AE_PBIAS_2071_2100_RCP45 <- ((AE_2071_2100 - AE_1961_1990) / AE_1961_1990)*100
AE_PBIAS_2071_2100_RCP45

#--------------------------------------------------------------------------------------------------------------
# SUBBLOCK: Data Processing RCP45
#--------------------------------------------------------------------------------------------------------------

# An optimization data.frame is created
df.opti_rcp45 <- BasinObs

# lubridate Library functions are applied to df_WCD to create new columns contaning:
# YEAR, YEAR_CH, MONTH, MONTH_CH
df.opti_rcp45$YEAR <- year(df.opti_rcp45$DatesR) # Years component of a date-time
df.opti_rcp45$YEAR_CH <- as.character(year(df.opti_rcp45$DatesR)) # Years component of a date-time as character
df.opti_rcp45$MONTH <- month(df.opti_rcp45$DatesR, label = FALSE) # Months component of a date-time
df.opti_rcp45$MONTH_CH <- month(df.opti_rcp45$DatesR, label = TRUE) #

# Relevant variables are selected
df.opti_rcp45 <- df.opti_rcp45[, c(1,9)]

# NA values are deleted
df.opti_rcp45 <-  df.opti_rcp45[complete.cases(df.opti_rcp45), ]

# A dummy variable is created
df.opti_rcp45$value <- c("value")

# A new data.frame is created
df.opti_rcp45.plot <- as.data.frame(dcast(df.opti_rcp45, MONTH_CH  ~ value, mean, value.var = "Qmm"))

# A dummy variable is created
df.opti_rcp45.plot$class <- c("obs")

# **********************************************
# An modelling data.frame is created
df.modelling_rcp45 <- data.frame(c(OutputsModel$Qsim))
df.modelling_rcp45$Dates <- OutputsModel$DatesR

# data.frame names are replaced
names(df.modelling_rcp45) <- c("Qsim", "DatesR")

# lubridate Library functions are applied to df_WCD to create new columns contaning:
# YEAR, YEAR_CH, MONTH, MONTH_CH
df.modelling_rcp45$YEAR <- year(df.modelling_rcp45$DatesR) # Years component of a date-time
df.modelling_rcp45$YEAR_CH <- as.character(year(df.modelling_rcp45$DatesR)) # Years component of a date-time as character
df.modelling_rcp45$MONTH <- month(df.modelling_rcp45$DatesR, label = FALSE) # Months component of a date-time
df.modelling_rcp45$MONTH_CH <- month(df.modelling_rcp45$DatesR, label = TRUE) #

# Relevant variables are selected
df.modelling_rcp45 <- df.modelling_rcp45[, c(1,6)]

# A dummy variable is created
df.modelling_rcp45$value <- c("value")

# A new data.frame is created
df.modelling_rcp45 <- as.data.frame(dcast(df.modelling_rcp45, MONTH_CH  ~ value, mean, value.var = "Qsim"))

# A dummy variable is created
df.modelling_rcp45$class <- c("opti")

# A production data.frame is created
df.modelling_rcp45.cal <- df.modelling_rcp45

# **********************************************
# An simulation 1961_1990 data.frame is created
df.modelling_rcp45 <- data.frame(c(OutputsModel_1961_1990$Qsim))
df.modelling_rcp45$Dates <- OutputsModel_1961_1990$DatesR

# data.frame names are replaced
names(df.modelling_rcp45) <- c("Qsim", "DatesR")

# lubridate Library functions are applied to df_WCD to create new columns contaning:
# YEAR, YEAR_CH, MONTH, MONTH_CH
df.modelling_rcp45$YEAR <- year(df.modelling_rcp45$DatesR) # Years component of a date-time
df.modelling_rcp45$YEAR_CH <- as.character(year(df.modelling_rcp45$DatesR)) # Years component of a date-time as character
df.modelling_rcp45$MONTH <- month(df.modelling_rcp45$DatesR, label = FALSE) # Months component of a date-time
df.modelling_rcp45$MONTH_CH <- month(df.modelling_rcp45$DatesR, label = TRUE) #

# Relevant variables are selected
df.modelling_rcp45 <- df.modelling_rcp45[, c(1,6)]

# A dummy variable is created
df.modelling_rcp45$value <- c("value")

# A new data.frame is created
df.modelling_rcp45 <- as.data.frame(dcast(df.modelling_rcp45, MONTH_CH  ~ value, mean, value.var = "Qsim"))

# A dummy variable is created
df.modelling_rcp45$class <- c("1961_1990")

# A production data.frame is created
df.sim_1961_1990_rcp45 <- df.modelling_rcp45

# **********************************************
# An simulation 2011_2040 data.frame is created
df.modelling_rcp45 <- data.frame(c(OutputsModel_2011_2040$Qsim))
df.modelling_rcp45$Dates <- OutputsModel_2011_2040$DatesR

# data.frame names are replaced
names(df.modelling_rcp45) <- c("Qsim", "DatesR")

# lubridate Library functions are applied to df_WCD to create new columns contaning:
# YEAR, YEAR_CH, MONTH, MONTH_CH
df.modelling_rcp45$YEAR <- year(df.modelling_rcp45$DatesR) # Years component of a date-time
df.modelling_rcp45$YEAR_CH <- as.character(year(df.modelling_rcp45$DatesR)) # Years component of a date-time as character
df.modelling_rcp45$MONTH <- month(df.modelling_rcp45$DatesR, label = FALSE) # Months component of a date-time
df.modelling_rcp45$MONTH_CH <- month(df.modelling_rcp45$DatesR, label = TRUE) #

# Relevant variables are selected
df.modelling_rcp45 <- df.modelling_rcp45[, c(1,6)]

# A dummy variable is created
df.modelling_rcp45$value <- c("value")

# A new data.frame is created
df.modelling_rcp45 <- as.data.frame(dcast(df.modelling_rcp45, MONTH_CH  ~ value, mean, value.var = "Qsim"))

# A dummy variable is created
df.modelling_rcp45$class <- c("2011_2040")

# A production data.frame is created
df.sim_2011_2040_rcp45 <- df.modelling_rcp45

# **********************************************
# An simulation 2041_2070 data.frame is created
df.modelling_rcp45 <- data.frame(c(OutputsModel_2041_2070$Qsim))
df.modelling_rcp45$Dates <- OutputsModel_2041_2070$DatesR

# data.frame names are replaced
names(df.modelling_rcp45) <- c("Qsim", "DatesR")

# lubridate Library functions are applied to df_WCD to create new columns contaning:
# YEAR, YEAR_CH, MONTH, MONTH_CH
df.modelling_rcp45$YEAR <- year(df.modelling_rcp45$DatesR) # Years component of a date-time
df.modelling_rcp45$YEAR_CH <- as.character(year(df.modelling_rcp45$DatesR)) # Years component of a date-time as character
df.modelling_rcp45$MONTH <- month(df.modelling_rcp45$DatesR, label = FALSE) # Months component of a date-time
df.modelling_rcp45$MONTH_CH <- month(df.modelling_rcp45$DatesR, label = TRUE) #

# Relevant variables are selected
df.modelling_rcp45 <- df.modelling_rcp45[, c(1,6)]

# A dummy variable is created
df.modelling_rcp45$value <- c("value")

# A new data.frame is created
df.modelling_rcp45 <- as.data.frame(dcast(df.modelling_rcp45, MONTH_CH  ~ value, mean, value.var = "Qsim"))

# A dummy variable is created
df.modelling_rcp45$class <- c("2041_2070")

# A production data.frame is created
df.sim_2041_2070_rcp45 <- df.modelling_rcp45

# **********************************************
# An simulation 2071_2100 data.frame is created
df.modelling_rcp45 <- data.frame(c(OutputsModel_2071_2100$Qsim))
df.modelling_rcp45$Dates <- OutputsModel_2071_2100$DatesR

# data.frame names are replaced
names(df.modelling_rcp45) <- c("Qsim", "DatesR")

# lubridate Library functions are applied to df_WCD to create new columns contaning:
# YEAR, YEAR_CH, MONTH, MONTH_CH
df.modelling_rcp45$YEAR <- year(df.modelling_rcp45$DatesR) # Years component of a date-time
df.modelling_rcp45$YEAR_CH <- as.character(year(df.modelling_rcp45$DatesR)) # Years component of a date-time as character
df.modelling_rcp45$MONTH <- month(df.modelling_rcp45$DatesR, label = FALSE) # Months component of a date-time
df.modelling_rcp45$MONTH_CH <- month(df.modelling_rcp45$DatesR, label = TRUE) #

# Relevant variables are selected
df.modelling_rcp45 <- df.modelling_rcp45[, c(1,6)]

# A dummy variable is created
df.modelling_rcp45$value <- c("value")

# A new data.frame is created
df.modelling_rcp45 <- as.data.frame(dcast(df.modelling_rcp45, MONTH_CH  ~ value, mean, value.var = "Qsim"))

# A dummy variable is created
df.modelling_rcp45$class <- c("2071_2100")

# A production data.frame is created
df.sim_2071_2100_rcp45 <- df.modelling_rcp45

# **********************************************
# An export data.frame is created
df.export_rcp45 <- rbind(df.sim_1961_1990_rcp45, df.sim_2011_2040_rcp45, df.sim_2041_2070_rcp45, df.sim_2071_2100_rcp45)
#df.export_rcp45 <- rbind(df.opti_rcp45.plot, df.modelling_rcp45.cal, df.sim_1961_1990_rcp45, df.sim_2011_2040_rcp45, df.sim_2041_2070_rcp45, df.sim_2071_2100_rcp45)

# A repetition variable is created
df.export_rcp45$REP <- rep(df.sim_1961_1990_rcp45$value, 4)

# A ratio variable is created
df.export_rcp45$RATIO <- ((df.export_rcp45$value / df.export_rcp45$REP) -1)*100

# Irrelevant variables are deleted
df.export_rcp45 <- df.export_rcp45[, c(-4)]

# data.frame is converted to long format
df.export_rcp45 <- melt(df.export_rcp45)

#==================================
# Mass-Balance-Components RCP45
#==================================

# $PotEvap [numeric] series of input potential evapotranspiration [mm/month]
# $Precip  [numeric] series of input total precipitation [mm/month]
# $AE	     [numeric] series of actual evapotranspiration [mm/month]
# $Perc    [numeric] series of percolation (P2) [mm/month]
# $Qsim    [numeric] series of simulated discharge [mm/month]

# Historical
vector.hist <- c((sum(OutputsModel_1961_1990$Precip))/30,
                 (sum(OutputsModel_1961_1990$AE))/30,
                 (sum(OutputsModel_1961_1990$Perc))/30,
                 (sum(OutputsModel_1961_1990$Qsim))/30)

# RCP45 2011_2040
vector.rcp45.2011_2040 <- c((sum(OutputsModel_2011_2040$Precip))/30,
                            (sum(OutputsModel_2011_2040$AE))/30,
                            (sum(OutputsModel_2011_2040$Perc))/30,
                            (sum(OutputsModel_2011_2040$Qsim))/30)

# RCP45 2041_2070
vector.rcp45.2041_2070 <- c((sum(OutputsModel_2041_2070$Precip))/30,
                            (sum(OutputsModel_2041_2070$AE))/30,
                            (sum(OutputsModel_2041_2070$Perc))/30,
                            (sum(OutputsModel_2041_2070$Qsim))/30)

# RCP45 2071_2100
vector.rcp45.2071_2100 <- c((sum(OutputsModel_2071_2100$Precip))/30,
                            (sum(OutputsModel_2071_2100$AE))/30,
                            (sum(OutputsModel_2071_2100$Perc))/30,
                            (sum(OutputsModel_2071_2100$Qsim))/30)

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////
# BLOCK: Future Simulation RCP26
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////

#--------------------------------------------------------------------------------------------------------------
# SUBBLOCK: RCP26 2011-2040 Future Simulation RCP26
#--------------------------------------------------------------------------------------------------------------

# Watershed compiled files are loaded
BasinObs_2011_2040 <- read.table("airGR_morote_2011_2040_RCP26.txt",header=T,sep="\t",quote="")

# ET0 is converted from mm/day to mm/month
BasinObs_2011_2040$E <- BasinObs_2011_2040$E*30

# Historical date vectors are created
dates_obs_2011_2040 <- seq(as.Date("2011-01-01"), as.Date("2040-12-01"), by = "1 month")

# Dates are transformed to POSIX class
dates_obs_2011_2040 <- as.POSIXlt(dates_obs_2011_2040, format="%Y-%m-%d")

# Watershed compiled data.frame is generated
BasinObs_2011_2040$DatesR <- dates_obs_2011_2040

# GR2M InputsModel object is created
InputsModel_2011_2040 <- CreateInputsModel(FUN_MOD = RunModel_GR2M,
                                           DatesR = BasinObs_2011_2040$DatesR,
                                           Precip = BasinObs_2011_2040$P,
                                           PotEvap = BasinObs_2011_2040$E)

# GR2M run period is selected
Ind_Run_2011_2040 <- seq(which(format(BasinObs_2011_2040$DatesR, format = "%Y-%m")=="2011-01"),
                         which(format(BasinObs_2011_2040$DatesR, format = "%Y-%m")=="2040-12"))

# GR2M RunOptions object is created
RunOptions_2011_2040 <- CreateRunOptions(FUN_MOD = RunModel_GR2M,
                                         InputsModel = InputsModel_2011_2040,
                                         IndPeriod_Run = Ind_Run_2011_2040)

# GR2M OutputsModel object is created
OutputsModel_2011_2040 <- RunModel_GR2M(InputsModel = InputsModel_2011_2040,
                                        RunOptions = RunOptions_2011_2040,
                                        Param = Param)

# A GR2M summary plot is requested
plot(OutputsModel_2011_2040, Qobs = BasinObs_2011_2040$Qmm[Ind_Run_2011_2040],
     cex.axis =0.9, cex.lab = 0.9, cex.leg = 0.9,
     BasinArea=929.40, which = c("all"))

# Mass balance of the period is calculated
mass_2011_2040 <- sum(OutputsModel_2011_2040$Qsim)
AE_2011_2040 <- sum(OutputsModel_2011_2040$AE)

# Mass balance of the period is requested
mass_2011_2040

#--------------------------------------------------------------------------------------------------------------
# SUBBLOCK: RCP26 2041_2070 Future Simulation RCP26
#--------------------------------------------------------------------------------------------------------------

# Watershed compiled files are loaded
BasinObs_2041_2070 <- read.table("airGR_morote_2041_2070_RCP26.txt",header=T,sep="\t",quote="")

# ET0 is converted from mm/day to mm/month
BasinObs_2041_2070$E <- BasinObs_2041_2070$E*30

# Historical date vectors are created
dates_obs_2041_2070 <- seq(as.Date("2041-01-01"), as.Date("2070-12-01"), by = "1 month")

# Dates are transformed to POSIX class
dates_obs_2041_2070 <- as.POSIXlt(dates_obs_2041_2070, format="%Y-%m-%d")

# Watershed compiled data.frame is generated
BasinObs_2041_2070$DatesR <- dates_obs_2041_2070

# GR2M InputsModel object is created
InputsModel_2041_2070 <- CreateInputsModel(FUN_MOD = RunModel_GR2M,
                                           DatesR = BasinObs_2041_2070$DatesR,
                                           Precip = BasinObs_2041_2070$P,
                                           PotEvap = BasinObs_2041_2070$E)

# GR2M run period is selected
Ind_Run_2041_2070 <- seq(which(format(BasinObs_2041_2070$DatesR, format = "%Y-%m")=="2041-01"),
                         which(format(BasinObs_2041_2070$DatesR, format = "%Y-%m")=="2070-12"))

# GR2M RunOptions object is created
RunOptions_2041_2070 <- CreateRunOptions(FUN_MOD = RunModel_GR2M,
                                         InputsModel = InputsModel_2041_2070,
                                         IndPeriod_Run = Ind_Run_2041_2070)

# GR2M OutputsModel object is created
OutputsModel_2041_2070 <- RunModel_GR2M(InputsModel = InputsModel_2041_2070,
                                        RunOptions = RunOptions_2041_2070,
                                        Param = Param)

# A GR2M summary plot is requested
plot(OutputsModel_2041_2070, Qobs = BasinObs_2041_2070$Qmm[Ind_Run_2041_2070],
     cex.axis =0.9, cex.lab = 0.9, cex.leg = 0.9,
     BasinArea=929.40, which = c("all"))

# Mass balance of the period is calculated
mass_2041_2070 <- sum(OutputsModel_2041_2070$Qsim)
AE_2041_2070 <- sum(OutputsModel_2041_2070$AE)

# Mass balance of the period is requested
mass_2041_2070

#--------------------------------------------------------------------------------------------------------------
# SUBBLOCK: RCP26 2071_2100 Future Simulation RCP26
#--------------------------------------------------------------------------------------------------------------

# Watershed compiled files are loaded
BasinObs_2071_2100 <- read.table("airGR_morote_2071_2100_RCP26.txt",header=T,sep="\t",quote="")

# ET0 is converted from mm/day to mm/month
BasinObs_2071_2100$E <- BasinObs_2071_2100$E*30

# Historical date vectors are created
dates_obs_2071_2100 <- seq(as.Date("2071-01-01"), as.Date("2100-12-01"), by = "1 month")

# Dates are transformed to POSIX class
dates_obs_2071_2100 <- as.POSIXlt(dates_obs_2071_2100, format="%Y-%m-%d")

# Watershed compiled data.frame is generated
BasinObs_2071_2100$DatesR <- dates_obs_2071_2100

# GR2M InputsModel object is created
InputsModel_2071_2100 <- CreateInputsModel(FUN_MOD = RunModel_GR2M,
                                           DatesR = BasinObs_2071_2100$DatesR,
                                           Precip = BasinObs_2071_2100$P,
                                           PotEvap = BasinObs_2071_2100$E)

# GR2M run period is selected
Ind_Run_2071_2100 <- seq(which(format(BasinObs_2071_2100$DatesR, format = "%Y-%m")=="2071-01"),
                         which(format(BasinObs_2071_2100$DatesR, format = "%Y-%m")=="2100-12"))

# GR2M RunOptions object is created
RunOptions_2071_2100 <- CreateRunOptions(FUN_MOD = RunModel_GR2M,
                                         InputsModel = InputsModel_2071_2100,
                                         IndPeriod_Run = Ind_Run_2071_2100)

# GR2M OutputsModel object is created
OutputsModel_2071_2100 <- RunModel_GR2M(InputsModel = InputsModel_2071_2100,
                                        RunOptions = RunOptions_2071_2100,
                                        Param = Param)

# A GR2M summary plot is requested
plot(OutputsModel_2071_2100, Qobs = BasinObs_2071_2100$Qmm[Ind_Run_2071_2100],
     cex.axis =0.9, cex.lab = 0.9, cex.leg = 0.9,
     BasinArea=929.40, which = c("all"))

# Mass balance of the period is calculated
mass_2071_2100 <- sum(OutputsModel_2071_2100$Qsim)
AE_2071_2100 <- sum(OutputsModel_2071_2100$AE)

# Mass balance of the period is requested
mass_2071_2100

# 2011_2040 PBIAS is calculated with respect to 1961-1990
PBIAS_2011_2040_RCP26 <- ((mass_2011_2040 - mass_1961_1990) / mass_1961_1990)*100
PBIAS_2011_2040_RCP26

# 2041_2070 PBIAS is calculated with respect to 1961-1990
PBIAS_2041_2070_RCP26 <- ((mass_2041_2070 - mass_1961_1990) / mass_1961_1990)*100
PBIAS_2041_2070_RCP26

# 2041_2070 PBIAS is calculated with respect to 1961-1990
PBIAS_2071_2100_RCP26 <- ((mass_2071_2100 - mass_1961_1990) / mass_1961_1990)*100
PBIAS_2071_2100_RCP26

# 2011_2040 AE_PBIAS is calculated with respect to 1961-1990
AE_PBIAS_2011_2040_RCP26 <- ((AE_2011_2040 - AE_1961_1990) / AE_1961_1990)*100
AE_PBIAS_2011_2040_RCP26

# 2041_2070 AE_PBIAS is calculated with respect to 1961-1990
AE_PBIAS_2041_2070_RCP26 <- ((AE_2041_2070 - AE_1961_1990) / AE_1961_1990)*100
AE_PBIAS_2041_2070_RCP26

# 2041_2070 AE_PBIAS is calculated with respect to 1961-1990
AE_PBIAS_2071_2100_RCP26 <- ((AE_2071_2100 - AE_1961_1990) / AE_1961_1990)*100
AE_PBIAS_2071_2100_RCP26

#--------------------------------------------------------------------------------------------------------------
# SUBBLOCK: Data Processing RCP26
#--------------------------------------------------------------------------------------------------------------

# An optimization data.frame is created
df.opti_rcp26 <- BasinObs

# lubridate Library functions are applied to df_WCD to create new columns contaning:
# YEAR, YEAR_CH, MONTH, MONTH_CH
df.opti_rcp26$YEAR <- year(df.opti_rcp26$DatesR) # Years component of a date-time
df.opti_rcp26$YEAR_CH <- as.character(year(df.opti_rcp26$DatesR)) # Years component of a date-time as character
df.opti_rcp26$MONTH <- month(df.opti_rcp26$DatesR, label = FALSE) # Months component of a date-time
df.opti_rcp26$MONTH_CH <- month(df.opti_rcp26$DatesR, label = TRUE) #

# Relevant variables are selected
df.opti_rcp26 <- df.opti_rcp26[, c(1,9)]

# NA values are deleted
df.opti_rcp26 <-  df.opti_rcp26[complete.cases(df.opti_rcp26), ]

# A dummy variable is created
df.opti_rcp26$value <- c("value")

# A new data.frame is created
df.opti_rcp26.plot <- as.data.frame(dcast(df.opti_rcp26, MONTH_CH  ~ value, mean, value.var = "Qmm"))

# A dummy variable is created
df.opti_rcp26.plot$class <- c("obs")

# **********************************************
# An modelling data.frame is created
df.modelling_rcp26 <- data.frame(c(OutputsModel$Qsim))
df.modelling_rcp26$Dates <- OutputsModel$DatesR

# data.frame names are replaced
names(df.modelling_rcp26) <- c("Qsim", "DatesR")

# lubridate Library functions are applied to df_WCD to create new columns contaning:
# YEAR, YEAR_CH, MONTH, MONTH_CH
df.modelling_rcp26$YEAR <- year(df.modelling_rcp26$DatesR) # Years component of a date-time
df.modelling_rcp26$YEAR_CH <- as.character(year(df.modelling_rcp26$DatesR)) # Years component of a date-time as character
df.modelling_rcp26$MONTH <- month(df.modelling_rcp26$DatesR, label = FALSE) # Months component of a date-time
df.modelling_rcp26$MONTH_CH <- month(df.modelling_rcp26$DatesR, label = TRUE) #

# Relevant variables are selected
df.modelling_rcp26 <- df.modelling_rcp26[, c(1,6)]

# A dummy variable is created
df.modelling_rcp26$value <- c("value")

# A new data.frame is created
df.modelling_rcp26 <- as.data.frame(dcast(df.modelling_rcp26, MONTH_CH  ~ value, mean, value.var = "Qsim"))

# A dummy variable is created
df.modelling_rcp26$class <- c("opti")

# A production data.frame is created
df.modelling_rcp26.cal <- df.modelling_rcp26

# **********************************************
# An simulation 1961_1990 data.frame is created
df.modelling_rcp26 <- data.frame(c(OutputsModel_1961_1990$Qsim))
df.modelling_rcp26$Dates <- OutputsModel_1961_1990$DatesR

# data.frame names are replaced
names(df.modelling_rcp26) <- c("Qsim", "DatesR")

# lubridate Library functions are applied to df_WCD to create new columns contaning:
# YEAR, YEAR_CH, MONTH, MONTH_CH
df.modelling_rcp26$YEAR <- year(df.modelling_rcp26$DatesR) # Years component of a date-time
df.modelling_rcp26$YEAR_CH <- as.character(year(df.modelling_rcp26$DatesR)) # Years component of a date-time as character
df.modelling_rcp26$MONTH <- month(df.modelling_rcp26$DatesR, label = FALSE) # Months component of a date-time
df.modelling_rcp26$MONTH_CH <- month(df.modelling_rcp26$DatesR, label = TRUE) #

# Relevant variables are selected
df.modelling_rcp26 <- df.modelling_rcp26[, c(1,6)]

# A dummy variable is created
df.modelling_rcp26$value <- c("value")

# A new data.frame is created
df.modelling_rcp26 <- as.data.frame(dcast(df.modelling_rcp26, MONTH_CH  ~ value, mean, value.var = "Qsim"))

# A dummy variable is created
df.modelling_rcp26$class <- c("1961_1990")

# A production data.frame is created
df.sim_1961_1990_rcp26 <- df.modelling_rcp26

# **********************************************
# An simulation 2011_2040 data.frame is created
df.modelling_rcp26 <- data.frame(c(OutputsModel_2011_2040$Qsim))
df.modelling_rcp26$Dates <- OutputsModel_2011_2040$DatesR

# data.frame names are replaced
names(df.modelling_rcp26) <- c("Qsim", "DatesR")

# lubridate Library functions are applied to df_WCD to create new columns contaning:
# YEAR, YEAR_CH, MONTH, MONTH_CH
df.modelling_rcp26$YEAR <- year(df.modelling_rcp26$DatesR) # Years component of a date-time
df.modelling_rcp26$YEAR_CH <- as.character(year(df.modelling_rcp26$DatesR)) # Years component of a date-time as character
df.modelling_rcp26$MONTH <- month(df.modelling_rcp26$DatesR, label = FALSE) # Months component of a date-time
df.modelling_rcp26$MONTH_CH <- month(df.modelling_rcp26$DatesR, label = TRUE) #

# Relevant variables are selected
df.modelling_rcp26 <- df.modelling_rcp26[, c(1,6)]

# A dummy variable is created
df.modelling_rcp26$value <- c("value")

# A new data.frame is created
df.modelling_rcp26 <- as.data.frame(dcast(df.modelling_rcp26, MONTH_CH  ~ value, mean, value.var = "Qsim"))

# A dummy variable is created
df.modelling_rcp26$class <- c("2011_2040")

# A production data.frame is created
df.sim_2011_2040_rcp26 <- df.modelling_rcp26

# **********************************************
# An simulation 2041_2070 data.frame is created
df.modelling_rcp26 <- data.frame(c(OutputsModel_2041_2070$Qsim))
df.modelling_rcp26$Dates <- OutputsModel_2041_2070$DatesR

# data.frame names are replaced
names(df.modelling_rcp26) <- c("Qsim", "DatesR")

# lubridate Library functions are applied to df_WCD to create new columns contaning:
# YEAR, YEAR_CH, MONTH, MONTH_CH
df.modelling_rcp26$YEAR <- year(df.modelling_rcp26$DatesR) # Years component of a date-time
df.modelling_rcp26$YEAR_CH <- as.character(year(df.modelling_rcp26$DatesR)) # Years component of a date-time as character
df.modelling_rcp26$MONTH <- month(df.modelling_rcp26$DatesR, label = FALSE) # Months component of a date-time
df.modelling_rcp26$MONTH_CH <- month(df.modelling_rcp26$DatesR, label = TRUE) #

# Relevant variables are selected
df.modelling_rcp26 <- df.modelling_rcp26[, c(1,6)]

# A dummy variable is created
df.modelling_rcp26$value <- c("value")

# A new data.frame is created
df.modelling_rcp26 <- as.data.frame(dcast(df.modelling_rcp26, MONTH_CH  ~ value, mean, value.var = "Qsim"))

# A dummy variable is created
df.modelling_rcp26$class <- c("2041_2070")

# A production data.frame is created
df.sim_2041_2070_rcp26 <- df.modelling_rcp26

# **********************************************
# An simulation 2071_2100 data.frame is created
df.modelling_rcp26 <- data.frame(c(OutputsModel_2071_2100$Qsim))
df.modelling_rcp26$Dates <- OutputsModel_2071_2100$DatesR

# data.frame names are replaced
names(df.modelling_rcp26) <- c("Qsim", "DatesR")

# lubridate Library functions are applied to df_WCD to create new columns contaning:
# YEAR, YEAR_CH, MONTH, MONTH_CH
df.modelling_rcp26$YEAR <- year(df.modelling_rcp26$DatesR) # Years component of a date-time
df.modelling_rcp26$YEAR_CH <- as.character(year(df.modelling_rcp26$DatesR)) # Years component of a date-time as character
df.modelling_rcp26$MONTH <- month(df.modelling_rcp26$DatesR, label = FALSE) # Months component of a date-time
df.modelling_rcp26$MONTH_CH <- month(df.modelling_rcp26$DatesR, label = TRUE) #

# Relevant variables are selected
df.modelling_rcp26 <- df.modelling_rcp26[, c(1,6)]

# A dummy variable is created
df.modelling_rcp26$value <- c("value")

# A new data.frame is created
df.modelling_rcp26 <- as.data.frame(dcast(df.modelling_rcp26, MONTH_CH  ~ value, mean, value.var = "Qsim"))

# A dummy variable is created
df.modelling_rcp26$class <- c("2071_2100")

# A production data.frame is created
df.sim_2071_2100_rcp26 <- df.modelling_rcp26

# **********************************************
# An export data.frame is created
df.export_rcp26 <- rbind(df.sim_1961_1990_rcp26, df.sim_2011_2040_rcp26, df.sim_2041_2070_rcp26, df.sim_2071_2100_rcp26)
#df.export_rcp26 <- rbind(df.opti_rcp26.plot, df.modelling_rcp26.cal, df.sim_1961_1990_rcp26, df.sim_2011_2040_rcp26, df.sim_2041_2070_rcp26, df.sim_2071_2100_rcp26)

# A repetition variable is created
df.export_rcp26$REP <- rep(df.sim_1961_1990_rcp26$value, 4)

# A ratio variable is created
df.export_rcp26$RATIO <- ((df.export_rcp26$value / df.export_rcp26$REP) -1)*100

# Irrelevant variables are deleted
df.export_rcp26 <- df.export_rcp26[, c(-4)]

# data.frame is converted to long format
df.export_rcp26 <- melt(df.export_rcp26)

#==================================
# Mass-Balance-Components RCP26
#==================================

# $PotEvap [numeric] series of input potential evapotranspiration [mm/month]
# $Precip  [numeric] series of input total precipitation [mm/month]
# $AE	     [numeric] series of actual evapotranspiration [mm/month]
# $Perc    [numeric] series of percolation (P2) [mm/month]
# $Qsim    [numeric] series of simulated discharge [mm/month]

# Historical
vector.hist <- c((sum(OutputsModel_1961_1990$Precip))/30,
                 (sum(OutputsModel_1961_1990$AE))/30,
                 (sum(OutputsModel_1961_1990$Perc))/30,
                 (sum(OutputsModel_1961_1990$Qsim))/30)

# RCP26 2011_2040
vector.rcp26.2011_2040 <- c((sum(OutputsModel_2011_2040$Precip))/30,
                            (sum(OutputsModel_2011_2040$AE))/30,
                            (sum(OutputsModel_2011_2040$Perc))/30,
                            (sum(OutputsModel_2011_2040$Qsim))/30)

# RCP26 2041_2070
vector.rcp26.2041_2070 <- c((sum(OutputsModel_2041_2070$Precip))/30,
                            (sum(OutputsModel_2041_2070$AE))/30,
                            (sum(OutputsModel_2041_2070$Perc))/30,
                            (sum(OutputsModel_2041_2070$Qsim))/30)

# RCP26 2071_2100
vector.rcp26.2071_2100 <- c((sum(OutputsModel_2071_2100$Precip))/30,
                            (sum(OutputsModel_2071_2100$AE))/30,
                            (sum(OutputsModel_2071_2100$Perc))/30,
                            (sum(OutputsModel_2071_2100$Qsim))/30)

#--------------------------------------------------------------------------------------------------------------
# SUBBLOCK: ggplot2
#--------------------------------------------------------------------------------------------------------------
# A color palette is created based on http://colorbrewer2.org/
scale04 <- c("#010101", "#4daf4a", "#377eb8", "#e41a1c")

# Absolute flow values (mm/month) are selected
df.export_hist_select <- df.export_rcp26[1:12, ]
df.export_rcp26_select <- df.export_rcp26[1:48, ]
df.export_rcp45_select <- df.export_rcp45[1:48, ]
df.export_rcp85_select <- df.export_rcp85[1:48, ]

# RCP labels are added
df.export_hist_select$RCP <- c("HIST")
df.export_rcp26_select$RCP <- c("RCP26")
df.export_rcp45_select$RCP <- c("RCP45")
df.export_rcp85_select$RCP <- c("RCP85")

# A compiled data.frame is created

df.export_total_select <- rbind(#df.export_hist_select,
                                df.export_rcp26_select,
                                df.export_rcp45_select,
                                df.export_rcp85_select)

# A comparison hydrograph is created
g01 <- ggplot() +
        geom_point(aes(x = MONTH_CH,y = value,shape = class,colour = class),data=df.export_total_select,size = 3.5) +
        facet_grid(facets = RCP ~ .) +
        geom_line(aes(x = MONTH_CH,y = value,colour = class,group = class),data=df.export_total_select,size = 1.2) +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
        scale_color_manual(values = scale04) +
        scale_shape_manual(values = c(16, 15, 18, 17))+
        #ggtitle("G01. morote. To be used") +
        xlab("Month of the year") +
        ylab("Q (mm/month)") +
        theme_bw() +
        theme(axis.text.x = element_text(vjust = 0.5,angle = 0.0), 
              text=element_text(size=20,  family="serif"), legend.position="bottom")

# A comparison hydrograph is requested
g01

# A comparison hydrograph is created
g02 <- ggplot() +
  geom_point(aes(x = MONTH_CH,y = value,shape = class,colour = class),data=df.export_rcp85,size = 3.5) +
  geom_line(aes(x = MONTH_CH,y = value,colour = class,linetype = variable, group = class),data=df.export_rcp85,size = 1.2 ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  facet_grid(facets = variable ~ ., scales = 'free_y') +
  scale_color_manual(values = scale04) +
  scale_shape_manual(values = c(16, 15, 18, 17))+
  ggtitle("G02. morote. To be used") +
  xlab("Month of the year") +
  ylab("Q (mm/month)") +
  theme_bw()

# A comparison hydrograph is requested
g02

# RCP mass balance vectors are created
Period_2011_2040 <- c(PBIAS_2011_2040_RCP26, PBIAS_2011_2040_RCP45, PBIAS_2011_2040_RCP85)
Period_2041_2070 <- c(PBIAS_2041_2070_RCP26, PBIAS_2041_2070_RCP45, PBIAS_2041_2070_RCP85)
Period_2071_1200 <- c(PBIAS_2071_2100_RCP26, PBIAS_2071_2100_RCP45, PBIAS_2071_2100_RCP85)

# A mass balance data.frame is created
df.mass.balance <- data.frame(Period_2011_2040, Period_2041_2070, Period_2071_1200)

# data.frame rownames are replaced
row.names(df.mass.balance) <- c("RCP26", "RCP45", "RCP85")

# A data.frame is requested
View(df.mass.balance)

# RCP mass balance vectors are created
AE_Period_2011_2040 <- c(AE_PBIAS_2011_2040_RCP26, AE_PBIAS_2011_2040_RCP45, AE_PBIAS_2011_2040_RCP85)
AE_Period_2041_2070 <- c(AE_PBIAS_2041_2070_RCP26, AE_PBIAS_2041_2070_RCP45, AE_PBIAS_2041_2070_RCP85)
AE_Period_2071_1200 <- c(AE_PBIAS_2071_2100_RCP26, AE_PBIAS_2071_2100_RCP45, AE_PBIAS_2071_2100_RCP85)

# A mass balance data.frame is created
df.mass.balance_AE <- data.frame(AE_Period_2011_2040, AE_Period_2041_2070, AE_Period_2071_1200)

# data.frame rownames are replaced
row.names(df.mass.balance_AE) <- c("RCP26", "RCP45", "RCP85")

# A data.frame is requested
View(df.mass.balance_AE)

# Relevant data.frames are exported to *.CSV
write.csv(df.export_total_select,"morote_flux_airGR.csv")
write.csv(df.mass.balance,"morote_mass_balance_airGR.csv")

# Calibration dispersed data.frame is created
df.cal.disper <- data.frame(OutputsModel$Qsim, BasinObs$Qmm[Ind_Run])

# data.frame variable names are replaced
names(df.cal.disper) <- c("Qsim", "Qobs")

# Stage variable is added
df.cal.disper$stage <- c("Calibration")

# Catchment variable is added
df.cal.disper$Catchment <- c("Morote")

# Validation dispersed data.frame is created
df.val.disper <- data.frame(OutputsModel02$Qsim, BasinObs$Qmm[Ind_Run02])

# data.frame variable names are replaced
names(df.val.disper) <- c("Qsim", "Qobs")

# Stage variable is added
df.val.disper$stage <- c("Validation")

# Catchment variable is added
df.val.disper$Catchment <- c("Morote")

# A compiled data.frame is created
df.calval.disper <- rbind(df.cal.disper, df.val.disper)

# Relevant data.frames are exported to *.CSV
write.csv(df.calval.disper,"Morote_calval_airGR.csv")

# RCP85 Compiled mass-balance data.frame
df.mass.rcp85 <- data.frame(vector.hist,
                            vector.rcp85.2011_2040,
                            vector.rcp85.2041_2070,
                            vector.rcp85.2071_2100)

# RCP45 Compiled mass-balance data.frame
df.mass.rcp45 <- data.frame(vector.hist,
                            vector.rcp45.2011_2040,
                            vector.rcp45.2041_2070,
                            vector.rcp45.2071_2100)

# RCP26 Compiled mass-balance data.frame
df.mass.rcp26 <- data.frame(vector.hist,
                            vector.rcp26.2011_2040,
                            vector.rcp26.2041_2070,
                            vector.rcp26.2071_2100)

#--------------------------------------------------------------------------------------------------------------
# SUBBLOCK: Calibration/Validation vs Observation Hydrographs
#--------------------------------------------------------------------------------------------------------------

# **********************************************
# An observed 1961_1990 data.frame is created
df.modelling_obs <- data.frame(c(BasinObs$Qmm))
df.modelling_obs$Dates <- BasinObs$DatesR

# data.frame names are replaced
names(df.modelling_obs) <- c("Qobs", "DatesR")
df.modelling_obs <- df.modelling_obs[c(Ind_Run, Ind_Run02), ]

# lubridate Library functions are applied to df_WCD to create new columns contaning:
# YEAR, YEAR_CH, MONTH, MONTH_CH
df.modelling_obs$YEAR <- year(df.modelling_obs$DatesR) # Years component of a date-time
df.modelling_obs$YEAR_CH <- as.character(year(df.modelling_obs$DatesR)) # Years component of a date-time as character
df.modelling_obs$MONTH <- month(df.modelling_obs$DatesR, label = FALSE) # Months component of a date-time
df.modelling_obs$MONTH_CH <- month(df.modelling_obs$DatesR, label = TRUE) #

# Relevant variables are selected
df.modelling_obs <- df.modelling_obs[, c(1,6)]

# A dummy variable is created
df.modelling_obs$value <- c("value")

# A new data.frame is created
df.modelling_obs <- as.data.frame(dcast(df.modelling_obs, MONTH_CH  ~ value, mean, value.var = "Qobs", na.rm = TRUE))

# A dummy variable is created
df.modelling_obs$class <- c("Obs")

# **********************************************
# An calibration 1961_1990 data.frame is created
df.modelling_cal <- data.frame(c(OutputsModel$Qsim))
df.modelling_cal$Dates <- OutputsModel$DatesR

# data.frame names are replaced
names(df.modelling_cal) <- c("Qobs", "DatesR")

# lubridate Library functions are applied to df_WCD to create new columns contaning:
# YEAR, YEAR_CH, MONTH, MONTH_CH
df.modelling_cal$YEAR <- year(df.modelling_cal$DatesR) # Years component of a date-time
df.modelling_cal$YEAR_CH <- as.character(year(df.modelling_cal$DatesR)) # Years component of a date-time as character
df.modelling_cal$MONTH <- month(df.modelling_cal$DatesR, label = FALSE) # Months component of a date-time
df.modelling_cal$MONTH_CH <- month(df.modelling_cal$DatesR, label = TRUE) #

# Relevant variables are selected
df.modelling_cal <- df.modelling_cal[, c(1,6)]

# A dummy variable is created
df.modelling_cal$value <- c("value")

# A new data.frame is created
df.modelling_cal <- as.data.frame(dcast(df.modelling_cal, MONTH_CH  ~ value, mean, value.var = "Qobs", na.rm = TRUE))

# A dummy variable is created
df.modelling_cal$class <- c("Cal")

# **********************************************
# An validation 1961_1990 data.frame is created
df.modelling_val <- data.frame(c(OutputsModel02$Qsim))
df.modelling_val$Dates <- OutputsModel02$DatesR

# data.frame names are replaced
names(df.modelling_val) <- c("Qobs", "DatesR")

# lubridate Library functions are applied to df_WCD to create new columns contaning:
# YEAR, YEAR_CH, MONTH, MONTH_CH
df.modelling_val$YEAR <- year(df.modelling_val$DatesR) # Years component of a date-time
df.modelling_val$YEAR_CH <- as.character(year(df.modelling_val$DatesR)) # Years component of a date-time as character
df.modelling_val$MONTH <- month(df.modelling_val$DatesR, label = FALSE) # Months component of a date-time
df.modelling_val$MONTH_CH <- month(df.modelling_val$DatesR, label = TRUE) #

# Relevant variables are selected
df.modelling_val <- df.modelling_val[, c(1,6)]

# A dummy variable is created
df.modelling_val$value <- c("value")

# A new data.frame is created
df.modelling_val <- as.data.frame(dcast(df.modelling_val, MONTH_CH  ~ value, mean, value.var = "Qobs", na.rm = TRUE))

# A dummy variable is created
df.modelling_val$class <- c("Val")

# A compiled data.frame is created
df.modelling_hydro <- rbind(df.modelling_obs, df.modelling_cal, df.modelling_val)

# A catchment variable is created
df.modelling_hydro$Catchment <- c("Morote")

# Relevant data.frames are exported to *.CSV
write.csv(df.modelling_hydro,"Morote_calval_airGR_hydro.csv")

# Minimum July-Mass Balance is calculated per RCP
df.export_JUL <- subset(df.export_total_select, MONTH_CH == "Jul")
df.export_JUL$ratio <- ((df.export_JUL [, 4] / df.export_JUL[1,4])-1)*100
View(df.export_JUL)

# Maximum Octuber-Mass Balance is calculated per RCP
df.export_OCT <- subset(df.export_total_select, MONTH_CH == "Oct")
df.export_OCT$ratio <- ((df.export_OCT [, 4] / df.export_OCT[1,4])-1)*100
View(df.export_OCT)

#----------------------------------------------------------------------------------------------------
# END OF CODE
#----------------------------------------------------------------------------------------------------
