##################################################
################ActiPHVV Analysis#################
##################################################
##requires R version 3.3.0
##by Nasrin Bopp
##################################################

#load packages
#The FKF package is required for ActiPHV
library(FKF) 

#additional packages for graphical output
library(ggplot2)
library(reshape2) 

#################################
####Check input parameters#######
#################################

#check, whether all arguments are provided, else provide default
if(exists("d") == FALSE){stop("Please supply the input directory.", call.=FALSE)}
if(exists("o") == FALSE){stop("Please supply the output directory.", call.=FALSE)}
if(exists("f") == FALSE){stop("Please supply the framerate in seconds per frame.", call.=FALSE)}
if(exists("px") == FALSE){stop("Please supply the pixel scale in micrometer per pixel.", call.=FALSE)}
if(exists("l") == FALSE){stop("Please supply the minimum phase length.", call.=FALSE)}

#rename input parameters for readability
rawDataFolder = d
outputFolder = o
secondsPerFrame = f
convertToMicrometer = px*(1/secondsPerFrame)
minimumPhaseLength = l

#check, if the last character of the passed file paths is a "/". If not, add it.
test.rawDataFolder =  unname(unlist(strsplit(rawDataFolder,"")))[length(unname(unlist(strsplit(rawDataFolder,""))))]
test.outputFolder =  unname(unlist(strsplit(outputFolder,"")))[length(unname(unlist(strsplit(outputFolder,""))))]

if(test.rawDataFolder != "/"){
  rawDataFolder = paste(rawDataFolder,"/", sep = "")
}

if(test.outputFolder != "/"){
  outputFolder = paste(outputFolder, "/", sep = "")
}

#get the filenames of the tracking data 
imgDataNames = unname(unlist(strsplit(list.files(rawDataFolder), ".txt")))

#################################
#######Kalman parameters#########
#################################

#kalman.fixedVariance: variance assumed for all measured values
#The value chosen is based on the manual filament tracking accuracy (2*sd)
kalman.fixedVariance = 2.6 

#Other parameters required for Kalman filtering
kalman.dt = matrix(c(0,0,0,0), ncol = 1)
kalman.ct = matrix(c(0,0), ncol = 1)
kalman.Tt = array(c(1,0,0,0,
                    0,1,0,0,
                    1,0,1,0,
                    0,1,0,1), dim = c(4,4,1)) # = A state transition matrix
kalman.Zt = array(c(1,0,0,1,0,0,0,0), dim = c(2,4,1)) # = C matrix in measurement equation
#function for HHt matrix to easily change the added errors
kalman.makeHHTarray = function(noiseX,noiseY,noiseVeloX,noiseVeloY){
  theArray = array(c(noiseX,0,0,0,
                     0,noiseY,0,0,
                     0,0,noiseVeloX,0,
                     0,0,0,noiseVeloY), dim = c(4,4,1))
  return(theArray)
}

#################################
####Start ActiPHV Analysis#######
#################################

#initialize variables to save the analysis results
allVeloAndError = NULL
phaseResults = NULL
finalCuts.all = NULL

#Perform ActiPHV analysis on each txt file (tracking data) provided in the input directory 
for (i in 1:length(imgDataNames)){
  #load raw data
  input = paste(rawDataFolder,imgDataNames[i], ".txt", sep = "")
  all = read.table(input, header = TRUE)

  #save data into variables (only for readability)
  timeSteps = all[,1]*f
  xCoord = all[,2]
  yCoord = all[,3]
  
  #calculate velocity if there are at least 15 tracked coordinates
  #if there are not enough tracked coordinates, the filament is skipped
  if(length(yCoord) >= 15){
    velocity = sqrt((diff(xCoord)^2+diff(yCoord)^2))
  } else {
      next  
    }
  
    #Perform Kalman filtering
    kalman.a0 = c(xCoord[1], yCoord[1], mean(velocity, na.rm = T), mean(velocity, na.rm = T))
    kalman.P0 = matrix(c(max(velocity),0,0,0,
                         0,max(velocity),0,0,
                         0,0,0,0,
                         0,0,0,0), nrow = 4, ncol = 4, byrow = TRUE)
    kalman.HHt = kalman.makeHHTarray(4,4,2,2)
    kalman.GGt = array(c(kalman.fixedVariance,0,0,kalman.fixedVariance), dim = c(2,2,1))
    kalman.yt = t(matrix(c(xCoord,yCoord), ncol = 2, nrow = length(xCoord)))
    
    #The Kalman-filtered velocity signal is saved in the variable velocity.kalman
    velocity.kalman = fkf(kalman.a0, kalman.P0, kalman.dt, kalman.ct,
                          kalman.Tt, kalman.Zt, kalman.HHt, kalman.GGt,
                          kalman.yt, check.input = TRUE)
    velocity.kalman = sapply(1:ncol(velocity.kalman$att), function(x) sqrt((velocity.kalman$att[3,x])^2 + 
                                                       (velocity.kalman$att[4,x])^2))
    
    #The Kalman filter needs a few iterations to achieve good estimates
    #Therefore the first five values estimated by the filter are removed
    velocity.kalman[1:5] <- NA
    
    ######################################
    ######### CALCULATE PHASES ###########
    ######################################
    
    #To prevent that too many maxima are detected, a Gauss filter is applied
    coeff.gauss = c(0.063327,0.093095,0.122589,0.144599,0.152781,0.144599,0.122589,0.093095,0.063327) #sigma = 3, kernel = 9
    velocity.kalman.gauss = filter(velocity.kalman, coeff.gauss, method = "convolution", sides = 2, circular = FALSE)

    #**************************************************
    #Determine maxima of the gaussian filtered signal**
    #**************************************************
    
    #compute derivatives for filtered signals
    velocity.1stDev = filter(velocity.kalman.gauss, c(-1,0,1), method = "convolution", sides = 2, circular = FALSE)
    velocity.2ndDev = filter(velocity.1stDev, c(-1,0,1), method = "convolution", sides = 2, circular = FALSE)
  
    #search for zero crossings = change of sign between two values 
    velocity.1stDev.zeros = sapply(1:(length(velocity.1stDev) - 1), function(x) (velocity.1stDev[x]*velocity.1stDev[x+1]) < 0)
    velocity.1stDev.zeros = which(velocity.1stDev.zeros == TRUE) 
  
    #find sign of the second derivative at the zero crossing positions of the first derivative
    velocity.2ndDev.signs = NULL
    for(d in 1:(length(velocity.2ndDev) - 1)){
      if((is.na(velocity.2ndDev[d]) || is.na(velocity.2ndDev[d+1])) == TRUE){
        velocity.2ndDev.signs[d] = NA
        next
      } else {
        if((velocity.2ndDev[d] * velocity.2ndDev[d+1] < 0) == TRUE){
          velocity.2ndDev.signs[d] = 0
        } else {
          if(velocity.2ndDev[d] < 0){
            velocity.2ndDev.signs[d] = -1
          } else {
            velocity.2ndDev.signs[d] = 1
          }
        }
      } 
    }
    
    velocity.maxima = velocity.1stDev.zeros[which(velocity.2ndDev.signs[velocity.1stDev.zeros] < 0)]

    #*********************************************
    #Determine the final phases
    #*********************************************
    
    #The velocity signal is split into fractions at the maxima positions given in the variable velocity.maxima
    #For the merging algorithm, the first and the last position of the signal is added
    velocity.maxima = c(1,velocity.maxima,length(velocity.kalman))
    
    #vector to save the maxima positions at which the final cuts will be
    finalCuts = vector() 
    
    #The following code compares every velocity fraction to the following one starting with 
    #the first fraction. If there is no significant difference between the fractions, the 
    #fractions are merged and then compared to the following fraction. If there is a significant
    #difference, the position that separates the significantly different fractions is saved 
    #in the vector finalCuts
    
    if(length(velocity.maxima) > 2){
      #Initialization
      starter = 1  #starter = first position of the first fraction 
      stopper = 2 #stopper = last position of the first fraction/first position of the second fraction
      ender = 3 #ender = last position of the second fraction
      
      for(t in 1:(length(velocity.maxima) - 1)){
        firstPhase = velocity.kalman[velocity.maxima[starter]:velocity.maxima[stopper]]
        secondPhase = velocity.kalman[velocity.maxima[stopper]:velocity.maxima[ender]]
        phaseDiff = t.test(firstPhase,secondPhase,na.action=na.omit)
        if(phaseDiff$p.value > 0.05){
          starter = starter #value does not change
          stopper = ender
          ender = ender + 1
        } else {
          finalCuts = c(finalCuts, stopper)
          starter = stopper
          stopper = ender
          ender = ender + 1
        }
        if(ender == (length(velocity.maxima) + 1)){break}
      }
      }
      
    # The positions of the final cuts in the input vector are: velocity.maxima[finalCuts]
    # These positions are used to get the correct cut time steps from the input data
    finalCuts.tempResult = timeSteps[velocity.maxima[finalCuts]]
    
    finalCuts.all = rbind(finalCuts.all,cbind(finalCuts.tempResult,rep(imgDataNames[i],length(finalCuts))))  
    
    ######################################
    #### CALCULATE PHASE STATISTICS ######
    ######################################
    
  #calculate velocity +- var for each phase and save it in an output table
  phaseMean = phaseSD = phaseVar = phaseLength = NULL
  phaseMean.kalman = phaseVar.kalman = phaseSD.kalman = phaseLength.kalman.withoutNA = NULL
  
  if(length(velocity.maxima) > 2){
    finalCuts = c(1,velocity.maxima[finalCuts],length(velocity.kalman))
    for(c in 1:(length(finalCuts) - 1)){
      phaseMean = c(phaseMean, mean(velocity[finalCuts[c]:finalCuts[c+1]], na.rm = TRUE))
      phaseVar = c(phaseVar, var(velocity[finalCuts[c]:finalCuts[c+1]], na.rm = TRUE))
      phaseSD = c(phaseSD,sd(velocity[finalCuts[c]:finalCuts[c+1]], na.rm = TRUE))
      phaseLength = c(phaseLength, length(velocity[finalCuts[c]:finalCuts[c+1]]))
      
      phaseMean.kalman = c(phaseMean.kalman, mean(velocity.kalman[finalCuts[c]:finalCuts[c+1]], na.rm = TRUE))
      phaseVar.kalman = c(phaseVar.kalman, var(velocity.kalman[finalCuts[c]:finalCuts[c+1]], na.rm = TRUE))
      phaseSD.kalman = c(phaseSD.kalman,sd(velocity.kalman[finalCuts[c]:finalCuts[c+1]], na.rm = TRUE))
      
      #after Kalman filtering the first five estimated values were set to NA,
      #because the filter needs a few iteration to give good results
      #here we calculate the number of values left per phase after filtering
      temp = velocity.kalman[finalCuts[c]:finalCuts[c+1]]
      phaseLength.kalman.withoutNA = c(phaseLength.kalman.withoutNA, length(temp[!is.na(temp)]))
      }
  } else {
    phaseMean = mean(velocity, na.rm = TRUE)
    phaseVar = var(velocity, na.rm = TRUE)
    phaseSD = sd(velocity, na.rm = TRUE)
    phaseLength = length(velocity)
    
    phaseMean.kalman = mean(velocity.kalman, na.rm = TRUE)
    phaseVar.kalman = var(velocity.kalman, na.rm = TRUE)
    phaseSD.kalman = sd(velocity.kalman, na.rm = TRUE)
    phaseLength.kalman.withoutNA = length(velocity.kalman[!is.na(velocity.kalman)])
  }
  
  phaseResults = rbind(phaseResults, cbind(rep(imgDataNames[i],length(phaseMean)), phaseMean, phaseVar,
                                           phaseSD, phaseLength, phaseMean.kalman, phaseVar.kalman, phaseSD.kalman, phaseLength.kalman.withoutNA))

  #save velo data of all tracks in one dataframe
  allVeloAndError_part = cbind(timeSteps[2:length(xCoord)],velocity, 
                               velocity.kalman[1:length(velocity)],velocity.kalman.gauss[1:length(velocity)],
                               rep(imgDataNames[i],length(velocity)))
  allVeloAndError = rbind(allVeloAndError, allVeloAndError_part)
  }

######################################
#format dataframes
######################################
#velo results of different filters
allVeloAndError = as.data.frame(allVeloAndError)
names(allVeloAndError) = c("time","unfilteredVelo","filteredVelo","gaussKalmanVelo","ID")
allVeloAndError$ID = as.factor(allVeloAndError$ID)
allVeloAndError$filteredVeloMicromPerSec = as.numeric(as.character(allVeloAndError$filteredVelo))*convertToMicrometer

#format the final cut data
finalCuts.all = as.data.frame(finalCuts.all)
names(finalCuts.all) = c("timeStep","ID")

#phase results
phaseResults = as.data.frame(phaseResults)
names(phaseResults) = c("ID","unfilteredVelo","var_unfiltered","sd_unfiltered","phaseLength",
                        "filteredVelo","var_filtered","sd_filtered","phaseLength_afterFiltering")
phaseResults$phaseLength = as.numeric(as.character(phaseResults$phaseLength)) * secondsPerFrame
phaseResults$phaseLength_afterFiltering = as.numeric(as.character(phaseResults$phaseLength_afterFiltering)) * secondsPerFrame
phaseResults$veloMicromPerSec = as.numeric(as.character(phaseResults$filteredVelo))*convertToMicrometer
phaseResults$var_veloMicromPerSec = as.numeric(as.character(phaseResults$var_filtered))*convertToMicrometer

#remove phases < minimumPhaseLength
phaseResults = phaseResults[phaseResults$phaseLength_afterFiltering >= minimumPhaseLength,]

#sort columns to better order 
#phaseResults = phaseResults[,c(1,5,9,2:4,6:10)]
phaseResults = data.frame(phaseResults$ID,phaseResults$phaseLength,phaseResults$phaseLength_afterFiltering,
           phaseResults$veloMicromPerSec,phaseResults$var_veloMicromPerSec,
           phaseResults$filteredVelo,phaseResults$var_filtered,phaseResults$sd_filtered,
           phaseResults$unfilteredVelo,phaseResults$var_unfiltered,phaseResults$sd_unfiltered)
names(phaseResults) = c("ID","phaseLength","phaseLength_afterFiltering",
                        "veloMicromPerSec","var_veloMicromPerSec",
                        "filteredVelo","var_filtered","sd_filtered",
                        "unfilteredVelo","var_unfiltered","sd_unfiltered")

##calculate maximum velocity fraction
x = as.numeric(as.character(phaseResults$filteredVelo))
totalSum = sum(x, na.rm = T)
sortedVekt = sort(x, decreasing = T)
result = 0
i = 1
while(result < 0.1*totalSum){
  result = result + sortedVekt[i]
  i = i+1
}
result = mean(sortedVekt[1:i])
result_sd = sd(sortedVekt[1:i])
result.micromPerSec = result*convertToMicrometer
result_sd.micromPerSec = result_sd*convertToMicrometer
noOfFilInMaxFraction = i

######################################
############MAKE PLOTS################
######################################

#only plot the filament velocities which were included in the phaseResults
allVeloAndError = allVeloAndError[allVeloAndError$ID %in% phaseResults$ID,]

#Merge velocity and finalCuts
allVeloAndError$joinVariable = paste(allVeloAndError$time,allVeloAndError$ID, sep = "_")
finalCuts.all$joinVariable = paste(finalCuts.all$timeStep,finalCuts.all$ID, sep = "_")
plotData = merge(allVeloAndError,finalCuts.all, all.x = TRUE)

#convert from factor to numeric 
plotData$time = as.numeric(as.character(plotData$time))
plotData$kalman = as.numeric(as.character(plotData$filteredVelo)) * convertToMicrometer
plotData$gaussKalmanVelo = as.numeric(as.character(plotData$gaussKalmanVelo)) * convertToMicrometer
plotData$unfilteredVelo = as.numeric(as.character(plotData$unfilteredVelo)) * convertToMicrometer
plotData$timeStep = as.numeric(as.character(plotData$timeStep))

#Create plot with the txt-file titles as headings
this.plot = ggplot(plotData, aes(x = time, y = kalman)) +
    geom_line() +
    facet_wrap(~ID,ncol = 7) +
    #labs(title = "Trajectory Cuts") +
    xlab("Time [s]") + ylab("Velocity [µm/s]") +
    ylim(0,15) + theme_bw() +
    geom_vline(aes(xintercept = timeStep), color = "red")

suppressWarnings(ggsave(paste(outputFolder,"PhaseIdentification.pdf", sep =""), width = 30, height = 20, units = "cm"))

#Create the same plot without headings
this.plot = ggplot(plotData, aes(x = time, y = kalman)) +
  geom_line() +
  facet_wrap(~ID,ncol = 7) +
  #labs(title = "Trajectory Cuts") +
  xlab("Time [s]") + ylab("Velocity [µm/s]") +
  ylim(0,15) + theme_bw() +
  geom_vline(aes(xintercept = timeStep), color = "red") +
  theme(strip.background = element_blank(), strip.text.x = element_blank()) #remove facet labels

suppressWarnings(ggsave(paste(outputFolder,"PhaseIdentification_noHeadings.pdf", sep =""), width = 30, height = 20, units = "cm"))


this.plot = ggplot(plotData, aes(x = time)) +
    geom_line(aes(y = kalman),colour = "orange") +
    geom_line(aes(y = gaussKalmanVelo), colour = "blue") + 
    facet_wrap(~ID,ncol = 7) +
    #labs(title = "Trajectory Cuts") +
    xlab("Time [s]") + ylab("Velocity [µm/s]") +
    ylim(0,15) + theme_bw() +
    geom_vline(aes(xintercept = timeStep), color = "red")
suppressWarnings(ggsave(paste(outputFolder,"PhaseIdentification_gaussFilteredSignal.pdf", sep =""), width = 30, height = 20, units = "cm"))

this.plot = ggplot(plotData, aes(x = time, y = unfilteredVelo)) +
  geom_line() +
  facet_wrap(~ID,ncol = 7) +
  #labs(title = "Trajectory Cuts") +
  xlab("Time [s]") + ylab("Velocity [µm/s]") +
  ylim(0,15) + theme_bw() +
  geom_vline(aes(xintercept = timeStep), color = "red")

suppressWarnings(ggsave(paste(outputFolder,"PhaseIdentification_rawSignal.pdf", sep =""), width = 30, height = 20, units = "cm"))

######################################
############SAVE RESULTS##############
######################################
#make filenames
file.phaseResults = paste(outputFolder,"PhaseResults.txt", sep ="")
file.timeDiff = paste(outputFolder,"TimeDiffBetweenExtrema.txt", sep ="")
file.maxVeloPerFilament = paste(outputFolder,"MaxVeloPerFilament.txt", sep ="")
file.dataForPlot = paste(outputFolder,"DataForPlot.txt", sep ="")
file.maxFraction = paste(outputFolder,"veloOfMaxFraction.txt", sep ="")

write.table(phaseResults, file.phaseResults , sep = "\t", row.names = FALSE)
write.table(plotData,file.dataForPlot,row.names = F, sep = "\t")

#txt file with maximum phase results
fileConn<-file(file.maxFraction)
writeLines(c("Max. velo in [px/frame]",result, "\n",
             "SD of max. velo in [px/frame]", result_sd,"\n",
             "Max. velo in [µm/s]",result.micromPerSec,"\n",
             "SD of max. velo in [µm/s]",result_sd.micromPerSec,"\n",
             "No. of filaments in the fastest fraction:", noOfFilInMaxFraction), 
           fileConn)
close(fileConn)

print("Analysis finished. Please find the results in the output directory.")

#reset input parameters
d = rawDataFolder
o = outputFolder
f = secondsPerFrame
px = convertToMicrometer*f
l = minimumPhaseLength


