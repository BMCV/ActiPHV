#PLEASE DEFINE INPUT PARAMETERS

#file path to ActiPHVV.R
scriptPath = "C:/Users/Nasrin/Documents/MoBiMaster/Masterarbeit_frontiers/Code/"

#if this script is in the same directory as ActiPHVV.R, uncomment the following line:
#scriptPath = getwd()

#input directory (tracking data of 1 image sequence)
#input files should be tab-separated txt-files
d = "C:/Users/Nasrin/Documents/MoBiMaster/Masterarbeit_frontiers/Data/InputData/Example2/"

#output directory. If it does not exist yet, it will be created.
o = "C:/Users/Nasrin/Documents/MoBiMaster/Masterarbeit_frontiers/Test_full_analysis/Example2_gauss/"
ifelse(!dir.exists(o), dir.create(o), FALSE)

#framerate in seconds per frame
f = 0.1

#pixel scale in micrometer per pixel
px = 0.128

#minimum phase length in seconds
l = 1

#START ANALYSIS
setwd(scriptPath)
source("ActiPHVV_Abgabe.R")
