# Example for Jags-Ymet-XmetMulti-Mrobust.R 
#------------------------------------------------------------------------------- 
# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
#rm(list=ls())  # Careful! This clears all of R's memory!
setwd("/Users/carloscastillo/Documents/RMIT/Semester 4/Bayesian Statistics/Project/bayesian-final-project")

##myData = read.csv("Assignment2PropertyPrices.csv")
yName = "count" ; xName = c("Total_Persons_Persons",	"Median_age_persons",	"Median_mortgage_repay_monthly", 	"Median_tot_hhd_inc_weekly",  	"Average_household_size","P_PGrad_Deg_Total","P_GradDip_and_GradCert_Total" , "P_BachDeg_Total","Certificates","Not Sated","Not Sated","Prv_Trs","Non_Motor_Trs")
#yName = "SalePrice" ; xName = c("Area",	"Bedrooms",	"Bathrooms", 	"CarParks")

fileNameRoot = "Task6"
numSavedSteps = 5000 ; thinSteps=2

graphFileType = "eps" 
#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("Part_B_Info.R")
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
#startTime = proc.time()
xPred = c(250 ,	3 ,	2 ,1 , 0 )
#xPred = c(600 ,	2 ,	2 ,	1  )
colnames(xPred) = c("beta[1]","beta[2]","beta[3]","beta[4]","beta[5]","beta[6]","beta[7]","beta[8]","beta[9]","beta[10]","beta[11]","beta[12]","beta[13]")
#colnames(xPred) = c("beta[1]","beta[2]","beta[3]","beta[4]")
mcmcCoda = genMCMC( data=myData , xName=xName , yName=yName , 
                    numSavedSteps=numSavedSteps , thinSteps=thinSteps , 
                    saveName=fileNameRoot , xPred = xPred )
#stopTime = proc.time()
#duration = stopTime - startTime
#show(duration)
#------------------------------------------------------------------------------- 
# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) # get all parameter names
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName , 
            saveName=fileNameRoot , saveType=graphFileType )
}
graphics.off()
#------------------------------------------------------------------------------- 
# Get summary statistics of chain:

summaryInfo = smryMCMC( mcmcCoda , 
                        saveName=fileNameRoot  )
show(summaryInfo)
# Display posterior information:
plotMCMC( mcmcCoda , data=myData , xName=xName , yName=yName , 
          pairsPlot=TRUE , showCurve=FALSE ,
          saveName=fileNameRoot , saveType=graphFileType )
#------------------------------------------------------------------------------- 
