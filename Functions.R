# Working Functions

# Functions for the Thornthwaite-Mather
#
# 3 Functions to calculate SWE and excess when soil is drying, 
#   wetting, and wetting above capacity
#
soildrying<-function(AWprev,dP,AWC){
  AW<-AWprev*exp(dP/AWC)
  excess<-0.0
  c(AW,excess)
}
soil_wetting_above_capacity<-function(AWprev,dP,AWC){
  AW<-AWC
  excess<-AWprev+dP-AWC
  c(AW,excess)
}
soilwetting<-function(AWprev,dP,AWC){
  AW<-AWprev+dP
  excess<-0.0
  c(AW,excess)
}

#
# We now paste together the remainder of the Lab 2-4 solutions into a function 
# TMWBModel() which will be called for the TopSlope, MidSlope, and BotSlope 
# HRU Objects that in turn returns an object with the same structure. Everything between 
# "START OF MODEL FUNCTION" and "END OF MODEL FUNCTION" is the 
# long function we are building
#
#  START OF MODEL FUNCTION
# Start of Model
TMWB_Model=function(fnc_TMWB,fnc_slope=0, 
                    fnc_aspect=0,func_DAWC=.3,
                    func_z=1000,fnc_fcres=.3,TempBias=-3) {
  # Energy Balance based Snow Accumulation 
  # and Melt model from the EcoHydRology package.
  attach(fnc_TMWB)
  SNO_Energy=SnowMelt(date, P, MaxTemp+TempBias, 
                      MinTemp+TempBias, myflowgage$declat, 
                      slope = fnc_slope, aspect = fnc_aspect, tempHt = 1, 
                      windHt = 2, groundAlbedo = 0.25,
                      SurfEmissiv = 0.95, windSp = 2, forest = 0, 
                      startingSnowDepth_m = 0, startingSnowDensity_kg_m3=450)
  # We will update the -3 in the above to be a lapse rate adjustment
  detach(fnc_TMWB)
  fnc_TMWB$SNO=SNO_Energy$SnowWaterEq_mm
  fnc_TMWB$SNOmlt=SNO_Energy$SnowMelt_mm
  
  
  # Snow accumulation and melt only depend on the surface attributes and weather, and as such, can run at the beginning, independent of the daily calculated ET, TMWB, and the linear reservoir Storage Discharge (Qmm). 
  # Similarly, Potential ET (PET) only depends on the surface attributes and weather, and as such, can run at the beginning.
  
  attach(fnc_TMWB)
  fnc_TMWB$Albedo=.23
  fnc_TMWB$Albedo[fnc_TMWB$SNO>0]=.95
  PET=PET_fromTemp(Jday=(1+as.POSIXlt(date)$yday),Tmax_C = MaxTemp,Tmin_C = MinTemp,lat_radians = myflowgage$declat*pi/180) * 1000
  fnc_TMWB$PET=PET
  detach(fnc_TMWB)
  rm(list="PET")
  
  # Those processes that are dependent on prior days conditions, we run as a loop through each of the days.
  
  fnc_TMWB$AWC=func_DAWC*func_z
  # Oh, this we want to vary some of these around our watershed!
  fnc_TMWB$dP = 0 # Initializing Net Precipitation
  fnc_TMWB$ET = 0 # Initializing ET
  fnc_TMWB$AW = 0 # Initializing AW
  fnc_TMWB$Excess = 0 # Initializing Excess
  
  
  # Loop to calculate AW and Excess
  attach(fnc_TMWB)
  for (t in 2:length(AW)){
    # This is where ET and Net Precipitation is now calculated
    # Update this to reflect the ET model described above
    ET[t] = AW[t-1]/AWC[t-1]*PET[t]
    dP[t] = SNO_Energy$Rain_mm[t] - ET[t] + 
      SNO_Energy$SnowMelt_mm[t] + HillslopeAboveExcess[t]
    
    # From here onward, everything is the same as Week2’s lab
    if (dP[t]<=0) {
      values<-soildrying(AW[t-1],dP[t],AWC[t])
    } else if((dP[t]>0) & (AW[t-1]+dP[t])<=AWC[t]) {
      values<-soilwetting(AW[t-1],dP[t],AWC[t])
    } else {
      values<-soil_wetting_above_capacity(AW[t-1],dP[t],AWC[t])
    }
    AW[t]<-values[1]
    Excess[t]<-values[2]
  }
  fnc_TMWB$AW=AW
  fnc_TMWB$Excess=Excess
  fnc_TMWB$dP=dP
  fnc_TMWB$ET=ET
  detach(fnc_TMWB) # IMPORTANT TO DETACH
  rm(list=c("AW", "dP", "ET", "Excess"))
  
  fnc_TMWB$Qpred=NA
  fnc_TMWB$Qpred[1]=0
  fnc_TMWB$S=NA
  fnc_TMWB$S[1]=0
  
  fcres=fnc_fcres
  attach(fnc_TMWB)
  for (t in 2:length(Qpred)){
    S[t]=S[t-1]+Excess[t]     
    Qpred[t]=fcres*S[t]
    S[t]=S[t]-Qpred[t]
  }
  fnc_TMWB$S=S
  fnc_TMWB$Qpred=Qpred # UPDATE vector BEFORE DETACHING
  detach(fnc_TMWB) # IMPORTANT TO DETACH
  rm(list=c("Qpred", "S"))
  return(fnc_TMWB)
}



#CNModel
# 
CN_Model<-function(fnc_CNModel, CNavg = 75,IaFrac = 0.05,fnc_slope=0, 
                   fnc_aspect=0,func_DAWC=.3,func_z=1000,fnc_fcres=.3) {
  
  # Energy Balance based Snow Accumulation 
  # and Melt model from the EcoHydRology package.
  attach(fnc_CNModel)
  SNO_Energy=SnowMelt(date, P, MaxTemp-3, MinTemp-3, myflowgage$declat, 
                      slope = fnc_slope, aspect = fnc_aspect, tempHt = 1, 
                      windHt = 2, groundAlbedo = 0.25,SurfEmissiv = 0.95, windSp = 2, 
                      forest = 0, startingSnowDepth_m = 0,startingSnowDensity_kg_m3=450)
  # We will update the -3 in the above to be a lapse rate adjustment
  detach(fnc_CNModel)
  fnc_CNModel$SNO=SNO_Energy$SnowWaterEq_mm
  fnc_CNModel$SNOmlt=SNO_Energy$SnowMelt_mm
  fnc_CNModel$SnowfallWatEq_mm=SNO_Energy$SnowfallWatEq_mm
  fnc_CNModel$SnowMelt_mm=SNO_Energy$SnowMelt_mm
  attach(fnc_CNModel)
  fnc_CNModel$Albedo=.23
  fnc_CNModel$Albedo[fnc_CNModel$SNO>0]=.95
  PET=PET_fromTemp(Jday=(1+as.POSIXlt(date)$yday),
                   Tmax_C = MaxTemp,Tmin_C = MinTemp,
                   lat_radians = myflowgage$declat*pi/180) * 1000
  fnc_CNModel$PET=PET
  detach(fnc_CNModel)
  rm(list="PET")
  
  fnc_CNModel$AWC=func_DAWC*func_z
  # Oh, this we want to vary some of these around our watershed!
  fnc_CNModel$dP = 0 # Initializing Net Precipitation
  fnc_CNModel$ET = 0 # Initializing ET
  fnc_CNModel$AW = 0 # Initializing AW
  fnc_CNModel$Excess = 0 # Initializing Excess
  fnc_CNModel$S =0 # Initializing S
  fnc_CNModel$Qpred=0 # Initializing Qpred
  attach(fnc_CNModel)
  SSCNavg=(1000/CNavg-10)*25.4
  SSCN=SoilStorage(S_avg=SSCNavg, field_capacity=func_DAWC*.9,
                   soil_water_content=0.1*func_DAWC, porosity=func_DAWC)
  Ia_init=IaFrac*SSCN   
  fnc_CNModel$CNavg = CNavg
  fnc_CNModel$SSCNavg = SSCNavg
  fnc_CNModel$SSCN = SSCN
  detach(fnc_CNModel)
  rm(list=c("CNavg", "SSCN", "SSCNavg"))
  fnc_CNModel$Ia = Ia_init
  attach(fnc_CNModel)
  # Those processes that are dependant on prior days conditions, we run as a 
  # loop through each of the days.
  for (t in 2:length(AW)){
    ET[t] = AW[t-1]/AWC[t-1]*PET[t]
    # Calculating Net Precipitation which adds in slope above's Excess
    dP[t] = SNO_Energy$Rain_mm[t] - ET[t] + 
      SNO_Energy$SnowMelt_mm[t] + HillslopeAboveExcess[t]    # CN Solution
    # Is the soil saturated, and thus can't take more dP? 
    if (AW[t-1] + dP[t]>=AWC[t]){
      Excess[t]=AW[t-1] + dP[t] -AWC[t]
      AW[t]=AWC[t]
      # Otherwise, if dP is less than the initial abstraction? 
      # https://en.wikipedia.org/wiki/Runoff_curve_number#Definition
    } else if (dP[t]<=Ia[t]) {
      Excess[t]=0.0
      AW[t]=AW[t-1] + dP[t]
    } else {
      Excess[t]=(dP[t]-Ia[t])^2/(dP[t]-Ia[t]+SSCN[t])
      AW[t]=AW[t-1] + dP[t] -Excess[t]
    }
    S[t]=S[t-1]+Excess[t]
    Qpred[t]=fnc_fcres*S[t]
    S[t]=S[t]-Qpred[t]
  }
  fnc_CNModel$ET=ET
  fnc_CNModel$dP=dP
  fnc_CNModel$AW=AW
  fnc_CNModel$Excess=Excess
  fnc_CNModel$S=S
  fnc_CNModel$Qpred=Qpred # UPDATE vector BEFORE DETACHING
  rm(list=c("AW", "dP", "ET", "Excess", "Qpred", "S"))
  detach(fnc_CNModel)
  return(fnc_CNModel)
}


library(EcoHydRology)
package.skeleton("BSEHydroModels",list=c("soil_wetting_above_capacity",
                "soilwetting","soildrying","TMWB_Model","CN_Model"))


install.packages("BSEHydroModels", repos=NULL)

#HM5 


