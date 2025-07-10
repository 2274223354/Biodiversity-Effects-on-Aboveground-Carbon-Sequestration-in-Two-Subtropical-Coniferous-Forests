rm(list=ls())
readDir=c("E:/data")
setwd(readDir)

library(lmerTest)
library(piecewiseSEM)
library(data.table)
library(pROC)
library(plyr)
library(dplyr)
library(doParallel)
library(foreign)
##===============================================================================
## PD:PSEM ---------
##===============================================================================
####Cunninghamia lanceolata####
####PSEM_All:PSEM####
rm(list=ls())

data_sem2=read.csv( "data_C.csv")
colnames(data_sem2 )

log.norm.scale <- function(df) {
  
  return(log_df <-log(df+1)  )  
  
}

for (c in c("pd_Broadleaf","shan_Broadleaf",         
            "FDis","FDisB"  )) {
  data_sem2[, c] <- log.norm.scale(data_sem2[, c]) 
}

colnames(data_sem2 )

#scale
log.norm.scale <- function(log_df) {
  
  return(scale(log_df)) 
}

for (c in c(2:18)) {
  data_sem2[, c] <- log.norm.scale(data_sem2[, c]) 
}
data_sem2[,2:18] <- apply(data_sem2[,2:18], 2, as.numeric)


data_sem2 <- as.data.frame(data_sem2)
data_sem2[is.na(data_sem2)] <- 0

mod_220_All<-psem(
  lmer(AverageDBH~ Density+
         Age+DBH_CV +
         #Shannon +#FDis.y+ 
         #SR+SR_Bro+SR_Other+
         #PD+    
         #pd_Broadleaf+
         #shan_Broadleaf+ #pd +
         
         # MAT+
         #MAP+
         Altitude+ (1|RegionName)
       
       
       # SoilThickness +
       #HumusLayerThickness#LitterThickness 
       ,data=data_sem2),
  
  lmer(DBH_CV~ #AverageDBH+
         Age+#Density+ 
         
         Shannon + (1|RegionName)# +#FDis.y+ 
       #SR+SR_Bro+SR_Other+
       # PD+    
       #pd_Broadleaf+
       #shan_Broadleaf+ #pd +
       
       #MAT+MAP+
       #Altitude+
       #SoilThickness +HumusLayerThickness+LitterThickness
       ,data=data_sem2),
  
  lmer(Density~ 
        Age+
         Shannon +#FDis.y+ 
         #SR+SR_Bro+SR_Other+
         PD     + (1|RegionName)
       #pd_Broadleaf+
       #shan_Broadleaf+ #pd +
       
       #MAT+MAP+
       #Altitude+
       #SoilThickness +HumusLayerThickness+LitterThickness
       ,data=data_sem2),
  
  
  
  lmer(PD~ #MAT+
         #Age+
         #MAP+
         Shannon+
         Altitude+
         SoilThickness  +
         (1|RegionName)#HumusLayerThickness+LitterThickness
       , data=data_sem2),
  
  lmer(Shannon~# MAT+
       # MAP+
         Age + (1|RegionName)#+
       #Altitude+
       # SoilThickness #HumusLayerThickness+LitterThickness
       , data=data_sem2),
  
  lmer( AGC~ Age+AverageDBH+
          Density+ 
          DBH_CV +
          #Shannon +  
         
          PD +    
         # pd_Broadleaf+
         # shan_Broadleaf+ 
          
          #MAT+
          #MAP+
          #Altitude+
          SoilThickness  + (1|RegionName)#+HumusLayerThickness+LitterThickness
        ,data=data_sem2)#,
  #Shannon%~~%shan_Broadleaf,
  #PD%~~% pd_Broadleaf
  
)

summary(mod_220_All, .progressBar = FALSE) 

AIC(mod_220_All, AIC.type = "dsep",
    aicc = T )

fn<-paste("C_psem_mode_All.R")
save(mod_220_All,file=fn)

mod_bor=mod_220_All

# 
pz_bor <- lapply(mod_bor[1:(length(mod_bor) - 1)], function(x) { which(summary(x)$co[, grep('P', colnames(summary(x)$co))] > 0.05) })
newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "C_psem_Newmode_All.R"
save(newmod_bor, file = fn_new_bor)

####PSEM_High:PSEM #########################
rm(list=ls())

data_sem2=read.csv( "data_C.csv")
colnames(data_sem2 )

C_Plot= read.csv("data_C_HLS_LCS_forests.csv")
data_sem2  =merge(C_Plot,data_sem2,by.x="Plot", by.y="Plot" )

data_sem2 = data_sem2[data_sem2$Group == "2" , ]


log.norm.scale <- function(df) {
  
  return(log_df <-log(df+1)  )  
  
}

for (c in c("pd_Broadleaf","shan_Broadleaf",         
            "FDis","FDisB"  )) {
  data_sem2[, c] <- log.norm.scale(data_sem2[, c]) 
}

colnames(data_sem2 )

#scale
log.norm.scale <- function(log_df) {
  
  return(scale(log_df)) 
}

for (c in c(3:19)) {
  data_sem2[, c] <- log.norm.scale(data_sem2[, c]) 
}
data_sem2[,3:19] <- apply(data_sem2[,3:19], 2, as.numeric)

str(data_sem2)
data_sem2[is.na(data_sem2)] <- 0

mod_220_High<-psem(
  lmer(AverageDBH~ 
       #DBH_CV +
       #Shannon +
       #SR+SR_Bro+SR_Other+
       #PD+    
       #pd_Broadleaf+
       #shan_Broadleaf+ #pd +
       
       # MAT+
       #MAP+
       Altitude+ #SoilThickness +
         (1|RegionName)#+
     # HumusLayerThickness+LitterThickness 
     ,data=data_sem2),
  
  # lm(DBH_CV~ #AverageDBH+
  #      # Density+ 
  #      
  #      Shannon# +#FDis.y+ 
  #    #SR+SR_Bro+SR_Other+
  #    # PD+    
  #    #pd_Broadleaf+
  #    #shan_Broadleaf+ #pd +
  #    
  #    #MAT+MAP+
  #    #Altitude+
  #    #SoilThickness +HumusLayerThickness+LitterThickness
  #    ,data=data_sem2),
  
  
  # lm(PD~ #MAT+
  #     # MAP# + 
  #      Shannon+shan_Broadleaf#+
  #    #Altitude+
  #    #SoilThickness +HumusLayerThickness+LitterThickness
  #    , data=data_sem2),
  
  lmer( AGC~ AverageDBH+
        #Density+ 
        # DBH_CV +
        Shannon +#FDis.y+ 
       
        #PD+    
        #pd_Broadleaf+
        shan_Broadleaf +
      
       #MAT+
      MAP+
     
      #Altitude+
      #SoilThickness #+HumusLayerThickness+LitterThickness
      +(1|RegionName) #pd +
      ,data=data_sem2)
  
  
  
)

summary(mod_220_High, .progressBar = FALSE) 

AIC(mod_220_High, AIC.type = "dsep",
    aicc = T )

fn<-paste("C_psem_mode_High.R")
save(mod_220_High,file=fn)

mod_bor=mod_220_High

# 
pz_bor <- lapply(mod_bor[1:(length(mod_bor) - 1)], function(x) { which(summary(x)$co[, grep('P', colnames(summary(x)$co))] > 0.05) })
newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "C_psem_Newmode_High.R"
save(newmod_bor, file = fn_new_bor)


####PSEM_Low:PSEM #########################
rm(list=ls())

data_sem2=read.csv( "data_C.csv")
colnames(data_sem2 )

C_Plot= read.csv("data_C_HLS_LCS_forests.csv")
data_sem2  =merge(C_Plot,data_sem2,by.x="Plot", by.y="Plot" )

data_sem2 = data_sem2[data_sem2$Group == "1" , ]


log.norm.scale <- function(df) {
  
  return(log_df <-log(df+1)  )  
  
}

for (c in c("pd_Broadleaf","shan_Broadleaf",         
            "FDis","FDisB"  )) {
  data_sem2[, c] <- log.norm.scale(data_sem2[, c]) 
}

colnames(data_sem2 )

#scale
log.norm.scale <- function(log_df) {
  
  return(scale(log_df)) 
}

for (c in c(3:19)) {
  data_sem2[, c] <- log.norm.scale(data_sem2[, c]) 
}
data_sem2[,3:19] <- apply(data_sem2[,3:19], 2, as.numeric)

str(data_sem2)
data_sem2[is.na(data_sem2)] <- 0


mod_220_Low<-psem(
  lmer(AverageDBH~ 
       DBH_CV +(1|RegionName) #+
     #Shannon +#FDis.y+ 
     #SR+SR_Bro+SR_Other+
     #PD+    
     #pd_Broadleaf+
     #shan_Broadleaf#+ #pd +
     
     # MAT+
     #MAP+
     #Altitude#+
     # SoilThickness +HumusLayerThickness+LitterThickness 
     ,data=data_sem2),
  
  lmer(DBH_CV~ #AverageDBH+
       # Density+
       
       Shannon+(1|RegionName) # +#FDis.y+
     #SR+SR_Bro+SR_Other+
     # PD+
     #pd_Broadleaf+
     #shan_Broadleaf+ #pd +
     
     #MAT+MAP+
     #Altitude+
     #SoilThickness +HumusLayerThickness+LitterThickness
     ,data=data_sem2),
  
  lmer(PD~ #MAT+
       #MAP+
       #Altitude+
       Shannon+
       #SoilThickness +
       HumusLayerThickness+(1|RegionName) #+LitterThickness
     , data=data_sem2),
  
  lmer(Shannon~ #MAT+
       # MAP+
       
       Altitude+(1|RegionName) 
     #SoilThickness +HumusLayerThickness+LitterThickness
     , data=data_sem2),
  
  
  lmer( AGC~ AverageDBH+
        #Density+ 
        # DBH_CV +
        #Shannon# +#FDis.y+ 
        #SR+SR_Bro+SR_Other+
        PD+    
        # pd_Broadleaf+
        # shan_Broadleaf+ #pd +
        
        MAP+
        Altitude+
        SoilThickness +HumusLayerThickness+(1|RegionName) #+LitterThickness
      ,data=data_sem2)
  
  
)

summary(mod_220_Low, .progressBar = FALSE) 

AIC(mod_220_Low, AIC.type = "dsep",
    aicc = T )



fn<-paste("C_psem_mode_Low.R")
save(mod_220_Low,file=fn)

mod_bor=mod_220_Low

# 
pz_bor <- lapply(mod_bor[1:(length(mod_bor) - 1)], function(x) { which(summary(x)$co[, grep('P', colnames(summary(x)$co))] > 0.05) })
newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "C_psem_Newmode_Low.R"
save(newmod_bor, file = fn_new_bor)

#################################################################################
####Pinus massoniana####
####PSEM_All:PSEM####
rm(list=ls())

data_sem2=read.csv( "data_P.csv")
colnames(data_sem2 )

log.norm.scale <- function(df) {
  
  return(log_df <-log(df+1)  )  
  
}

for (c in c("pd_Broadleaf","shan_Broadleaf",         
            "FDis","FDisB"  )) {
  data_sem2[, c] <- log.norm.scale(data_sem2[, c]) 
}

colnames(data_sem2 )

#scale
log.norm.scale <- function(log_df) {
  
  return(scale(log_df)) 
}

for (c in c(2:18)) {
  data_sem2[, c] <- log.norm.scale(data_sem2[, c]) 
}
data_sem2[,2:18] <- apply(data_sem2[,2:18], 2, as.numeric)


data_sem2 <- as.data.frame(data_sem2)
data_sem2[is.na(data_sem2)] <- 0


mod_220_All<-psem(
  lmer(AverageDBH~ Density+
       Age+
       DBH_CV +
       Shannon +#FDis.y+ 
       #SR+SR_Bro+SR_Other+
       #PD+    
       #pd_Broadleaf+
       #shan_Broadleaf+ #pd +
       
       # MAT+
       #MAP+
       Altitude+ 
      SoilThickness +(1|RegionName)#+
     #HumusLayerThickness#LitterThickness 
     ,data=data_sem2),
  
  lmer(DBH_CV~ #AverageDBH+
       #Age+
       Density+ 
       #SR_Other+
       
       Shannon +#FDis.y+ 
       #SR+SR_Bro+SR_Other+
       PD+    
       #pd_Broadleaf+
       #shan_Broadleaf+ #pd +
       
       #MAT+MAP+
       #Altitude+
       SoilThickness+ (1|RegionName) #+HumusLayerThickness+LitterThickness
     ,data=data_sem2),
  
  lmer(Density~ Age +#+ #SR_Other +
       #Shannon +#FDis.y+ 
       #SR+SR_Bro+SR_Other+
       PD + (1|RegionName)  
     #pd_Broadleaf+
     #shan_Broadleaf+ #pd +
     
     #MAT+MAP+
     #Altitude+
     #SoilThickness #+HumusLayerThickness+LitterThickness
     ,data=data_sem2),
  
  
  lmer(PD~ #MAT+
       Age+Shannon+shan_Broadleaf+ (1|RegionName)#+#MAP+
     #Altitude+
     #SoilThickness #HumusLayerThickness+LitterThickness
     , data=data_sem2),
  
  lmer(Shannon~# MAT+
       # MAP
       Age+
       #Altitude+
       SoilThickness + (1|RegionName)#HumusLayerThickness+LitterThickness
     , data=data_sem2),
  lmer( AGC~ Age+AverageDBH+
        Density+ 
        DBH_CV +
        #Shannon +  #FDis.y+ 
       
        
       #PD +
        #Altitude+#SoilThickness +
        #pd_Broadleaf+
        shan_Broadleaf+ (1|RegionName)
      ,data=data_sem2),
  
  Shannon%~~%shan_Broadleaf
)

summary(mod_220_All, .progressBar = FALSE) 

AIC(mod_220_All, AIC.type = "dsep",
    aicc = T )
fn<-paste("P_psem_mode_All.R")
save(mod_220_All,file=fn)

mod_bor=mod_220_All

# 
pz_bor <- lapply(mod_bor[1:(length(mod_bor) - 1)], function(x) { which(summary(x)$co[, grep('P', colnames(summary(x)$co))] > 0.05) })
newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "P_psem_Newmode_All.R"
save(newmod_bor, file = fn_new_bor)

####PSEM_High:PSEM #########################
rm(list=ls())

data_sem2=read.csv( "data_P.csv")

C_Plot=read.csv("data_P_HLS_LCS_forests.csv")

C_Plot=C_Plot[,c("Plot","Group"  )]

data_sem2 = merge(C_Plot,data_sem2,by.x="Plot", by.y="Plot" )
colnames(data_sem2)

data_sem2 = data_sem2[data_sem2$Group == "2" , ]
colnames(data_sem2)


log.norm.scale <- function(df) {
  
  return(log_df <-log(df+1)  )  
  
}

for (c in c("pd_Broadleaf","shan_Broadleaf",         
            "FDis","FDisB"  )) {
  data_sem2[, c] <- log.norm.scale(data_sem2[, c]) 
}

colnames(data_sem2 )

#scale
log.norm.scale <- function(log_df) {
  
  return(scale(log_df)) 
}

for (c in c(3:19)) {
  data_sem2[, c] <- log.norm.scale(data_sem2[, c]) 
}
data_sem2[,3:19] <- apply(data_sem2[,3:19], 2, as.numeric)


data_sem2 <- as.data.frame(data_sem2)
data_sem2[is.na(data_sem2)] <- 0


mod_220_High <-psem(
  lmer(AverageDBH~ Density+
       #DBH_CV +
       #Shannon +#FDis.y+ 
       #SR+SR_Bro+SR_Other+
       PD   +(1|RegionName)
    
     ,data=data_sem2),
  
  lmer(DBH_CV~ #AverageDBH+
       #Density+ 
       #SR_Other+
       Shannon   +(1|RegionName)
     
     ,data=data_sem2),
  
  # lmer(Density~ #SR_Other +
  #     Shannon +#FDis.y+ 
  #      #SR+SR_Bro+SR_Other+
  #      #PD    
  #        +(1|RegionName)
  #    
  #    ,data=data_sem2),
  
  lmer(PD~ Shannon  +(1|RegionName)
     , data=data_sem2),
  
  
  lmer( AGC~ AverageDBH+
        #Density+ 
        DBH_CV +
        Shannon +  
          #pd_Broadleaf+
         # shan_Broadleaf+
        PD   +(1|RegionName)
      
      #MAT+
      # MAP+
      #Altitude+
      # SoilThickness #+HumusLayerThickness+LitterThickness
      ,data=data_sem2)
  
  
)

summary(mod_220_High, .progressBar = FALSE) 

AIC(mod_220_High, AIC.type = "dsep",
    aicc = T )

fn<-paste("P_psem_mode_High.R")
save(mod_220_High,file=fn)

mod_bor=mod_220_High

# 
pz_bor <- lapply(mod_bor[1:(length(mod_bor) - 1)], function(x) { which(summary(x)$co[, grep('P', colnames(summary(x)$co))] > 0.05) })
newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "P_psem_Newmode_High.R"
save(newmod_bor, file = fn_new_bor)

####PSEM_Low:PSEM #########################
rm(list=ls())

data_sem2=read.csv( "data_P.csv")

C_Plot=read.csv("data_P_HLS_LCS_forests.csv")

C_Plot=C_Plot[,c("Plot","Group"  )]

data_sem2 = merge(C_Plot,data_sem2,by.x="Plot", by.y="Plot" )
colnames(data_sem2)

data_sem2 = data_sem2[data_sem2$Group == "1" , ]
colnames(data_sem2)

log.norm.scale <- function(df) {
  
  return(log_df <-log(df+1)  )  
  
}

for (c in c("pd_Broadleaf","shan_Broadleaf",         
            "FDis","FDisB"  )) {
  data_sem2[, c] <- log.norm.scale(data_sem2[, c]) 
}

colnames(data_sem2 )

#scale
log.norm.scale <- function(log_df) {
  
  return(scale(log_df)) 
}

for (c in c(3:19)) {
  data_sem2[, c] <- log.norm.scale(data_sem2[, c]) 
}
data_sem2[,3:19] <- apply(data_sem2[,3:19], 2, as.numeric)


data_sem2 <- as.data.frame(data_sem2)
data_sem2[is.na(data_sem2)] <- 0

mod_220_Low <-psem(
  lmer(AverageDBH~ #Density+
       
       DBH_CV +
       #Shannon +#FDis.y+ 
       #SR+SR_Bro+SR_Other+
       PD +(1|RegionName)  
     #pd_Broadleaf+
     #shan_Broadleaf+ #pd +
     
     # MAT+
     #MAP+
     #Altitude#+
     # SoilThickness +
     #HumusLayerThickness#LitterThickness 
     ,data=data_sem2),
  
 
  lmer(PD~ #MAT+shan_Broadleaf+
       # MAP+
       Shannon+(1|RegionName)
    # Altitude+
    # SoilThickness+ HumusLayerThickness+LitterThickness+
     , data=data_sem2),
  
  lmer( AGC~ AverageDBH+
        Density+ 
        DBH_CV +(1|RegionName)#+
      #Shannon +  #FDis.y+ 
      #species_count+
      #SR_Bro+
      #SR_Other+
      #PD #+    
      #pd_Broadleaf+
      # shan_Broadleaf
      
      #MAT+
      # MAP+
      #Altitude+
      # SoilThickness #+HumusLayerThickness+LitterThickness
      ,data=data_sem2)
  
  
)

summary(mod_220_Low, .progressBar = FALSE) 

AIC(mod_220_Low, AIC.type = "dsep",
    aicc = T )

fn<-paste("P_psem_mode_Low.R")
save(mod_220_Low,file=fn)

mod_bor=mod_220_Low

# 
pz_bor <- lapply(mod_bor[1:(length(mod_bor) - 1)], function(x) { which(summary(x)$co[, grep('P', colnames(summary(x)$co))] > 0.05) })
newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    update(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "P_psem_Newmode_Low.R"
save(newmod_bor, file = fn_new_bor)

#################################################################################
#################################################################################
##===============================================================================
## FD:PSEM ---------
##===============================================================================
rm(list=ls())

####Cunninghamia lanceolata####
####PSEM_All:PSEM####
rm(list=ls())

data_sem2=read.csv( "data_C.csv")
colnames(data_sem2 )

log.norm.scale <- function(df) {
  
  return(log_df <-log(df+1)  )  
  
}

for (c in c("pd_Broadleaf","shan_Broadleaf",         
            "FDis","FDisB"  )) {
  data_sem2[, c] <- log.norm.scale(data_sem2[, c]) 
}

colnames(data_sem2 )

#scale
log.norm.scale <- function(log_df) {
  
  return(scale(log_df)) 
}

for (c in c(2:18)) {
  data_sem2[, c] <- log.norm.scale(data_sem2[, c]) 
}
data_sem2[,2:18] <- apply(data_sem2[,2:18], 2, as.numeric)


data_sem2 <- as.data.frame(data_sem2)
data_sem2[is.na(data_sem2)] <- 0

mod_220_All<-psem(
  lmer(AverageDBH~ Density+
       Age+DBH_CV +
       #Shannon +#FDis.y+ 
       #SR+SR_Bro+SR_Other+
       #FDis.y+    
       #FDis.y_Broadleaf+
       #shan_Broadleaf+ #FDis.y +
       
       # MAT+
       #MAP+
       Altitude+ (1|RegionName)#+
     # SoilThickness +
     #HumusLayerThickness#LitterThickness 
     ,data=data_sem2),
  
  lmer(DBH_CV~ #AverageDBH+
       Age+
         PD+#Density+
       # pd_Broadleaf+
       FDis+#FDisB+ 
     #Shannon# 
     #SR+SR_Bro+SR_Other+
     # FDis.y+    
     #FDis.y_Broadleaf+
     #shan_Broadleaf+ #FDis.y +
     
     #MAT+MAP+
     Altitude+(1|RegionName)#PD
     #SoilThickness +HumusLayerThickness+LitterThickness
     ,data=data_sem2),
  
  lmer(Density~ 
       Age+PD+#Shannon +#FDis.y+ 
       #SR+SR_Bro+SR_Other+
       FDis   + (1|RegionName) 
     #FDis.y_Broadleaf+
     #shan_Broadleaf+ #FDis.y +
     #pd_Broadleaf+
     #MAT+MAP+
     #Altitude+
     #SoilThickness #HumusLayerThickness+LitterThickness
     ,data=data_sem2),
  
  lmer(FDis~ Age +#pd_Broadleaf#+#MAT+
       #Age+#MAP+
       #Altitude+
     #Shannon+
     SoilThickness + (1|RegionName)
     #HumusLayerThickness+LitterThickness
     , data=data_sem2),
  
  lmer(PD~# MAT+
       MAP+
       Age+
       Altitude+
       #SoilThickness +
         (1|RegionName)#HumusLayerThickness+LitterThickness
     , data=data_sem2),
  
  
  
  lmer( AGC~ Age+AverageDBH+
        Density+ 
        DBH_CV +
        
      # FDis + #  FDisB+  
        #pd_Broadleaf+
        #shan_Broadleaf+
        #pd+
        PD+
        #MAT+
        #MAP+
        #Altitude+
        SoilThickness + (1|RegionName)#+HumusLayerThickness+LitterThickness
      ,data=data_sem2),
   PD%~~%FDis
 
  
)

summary(mod_220_All, .progressBar = FALSE) 

AIC(mod_220_All, AIC.type = "dsep",
    aicc = T )

fn<-paste("C_psem_mode_All_FD.R")
save(mod_220_All,file=fn)

mod_bor=mod_220_All

# 
pz_bor <- lapply(mod_bor[1:(length(mod_bor) - 1)], function(x) { which(summary(x)$co[, grep('P', colnames(summary(x)$co))] > 0.05) })
newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    uFDis.yate(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "C_psem_Newmode_All_FD.R"
save(newmod_bor, file = fn_new_bor)

####PSEM_High:PSEM #########################
rm(list=ls())

data_sem2=read.csv( "data_C.csv")
colnames(data_sem2 )

C_Plot= read.csv("data_C_HLS_LCS_forests.csv")
data_sem2  =merge(C_Plot,data_sem2,by.x="Plot", by.y="Plot" )

data_sem2 = data_sem2[data_sem2$Group == "2" , ]


log.norm.scale <- function(df) {
  
  return(log_df <-log(df+1)  )  
  
}

for (c in c("pd_Broadleaf","shan_Broadleaf",         
            "FDis","FDisB"  )) {
  data_sem2[, c] <- log.norm.scale(data_sem2[, c]) 
}

colnames(data_sem2 )

#scale
log.norm.scale <- function(log_df) {
  
  return(scale(log_df)) 
}

for (c in c(3:19)) {
  data_sem2[, c] <- log.norm.scale(data_sem2[, c]) 
}
data_sem2[,3:19] <- apply(data_sem2[,3:19], 2, as.numeric)

str(data_sem2)
data_sem2[is.na(data_sem2)] <- 0



mod_220_High<-psem(
  lmer(AverageDBH~ 
       #DBH_CV +
       #Shannon +#FDis.y+ 
       #SR+SR_Bro+SR_Other+
       #FDis.y+    
       #FDis.y_Broadleaf+
       #shan_Broadleaf+ #FDis.y +
       
       # MAT+
      # MAP+
       Altitude+ (1|RegionName)#+
     # SoilThickness +HumusLayerThickness+LitterThickness 
     ,data=data_sem2),
  
  lmer(DBH_CV~ #AverageDBH+
       # Density+ 
       
       #PD +
       FDis+ (1|RegionName)
     #SR+SR_Bro+SR_Other+
     # FDis.y+    
     #FDis.y_Broadleaf+
     #shan_Broadleaf+ #FDis.y +
     
     #MAT+MAP+
     #Altitude+
     #SoilThickness +HumusLayerThickness+LitterThickness
     ,data=data_sem2),
  
  
  
  # lmer(FDis~ MAT+
  #      MAP+
  #      #shan_Broadleaf 
  #        + (1|RegionName)
  #    
  #    , data=data_sem2),
  
  # lm(Shannon~ #MAT+
  #      # MAP
  #     
  #    #Altitude+
  #    #SoilThickness +HumusLayerThickness+LitterThickness
  #    , data=data_sem2),
  
  lmer( AGC~ AverageDBH+
        #Density+ 
        DBH_CV +
        #Shannon# +
        FDis+FDisB+ MAP+
        #SR
        # PD+
        # pd_Broadleaf
        #shan_Broadleaf+ 
          (1|RegionName)
      # MAT+
     
      #Altitude+
      # SoilThickness +HumusLayerThickness+LitterThickness
      ,data=data_sem2)
 
)

summary(mod_220_High, .progressBar = FALSE) 

AIC(mod_220_High, AIC.type = "dsep",
    aicc = T )

fn<-paste("C_psem_mode_High_FD.R")
save(mod_220_High,file=fn)

mod_bor=mod_220_High

# 
pz_bor <- lapply(mod_bor[1:(length(mod_bor) - 1)], function(x) { which(summary(x)$co[, grep('P', colnames(summary(x)$co))] > 0.05) })
newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    uFDis.yate(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "C_psem_Newmode_High_FD.R"
save(newmod_bor, file = fn_new_bor)


####PSEM_Low:PSEM #########################
rm(list=ls())

data_sem2=read.csv( "data_C.csv")
colnames(data_sem2 )

C_Plot= read.csv("data_C_HLS_LCS_forests.csv")
data_sem2  =merge(C_Plot,data_sem2,by.x="Plot", by.y="Plot" )

data_sem2 = data_sem2[data_sem2$Group == "1" , ]


log.norm.scale <- function(df) {
  
  return(log_df <-log(df+1)  )  
  
}

for (c in c("pd_Broadleaf","shan_Broadleaf",         
            "FDis","FDisB"  )) {
  data_sem2[, c] <- log.norm.scale(data_sem2[, c]) 
}

colnames(data_sem2 )

#scale
log.norm.scale <- function(log_df) {
  
  return(scale(log_df)) 
}

for (c in c(3:19)) {
  data_sem2[, c] <- log.norm.scale(data_sem2[, c]) 
}
data_sem2[,3:19] <- apply(data_sem2[,3:19], 2, as.numeric)

str(data_sem2)
data_sem2[is.na(data_sem2)] <- 0



mod_220_Low<-psem(
  lmer(AverageDBH~ 
       DBH_CV + (1|RegionName)#+
     #Shannon +#FDis.y+ 
     #SR+SR_Bro+SR_Other+
     #FDis.y+    
     #FDis.y_Broadleaf+
     #shan_Broadleaf+ #FDis.y +
     
     # MAT+
     #MAP+
     #Altitude#+
     # SoilThickness +HumusLayerThickness+LitterThickness 
     ,data=data_sem2),
  
  # lm(FDis.y~ #MAT+
  #     #MAP+
  #     #Altitude+
  #      Shannon
  #   #SoilThickness +HumusLayerThickness+LitterThickness
  #   , data=data_sem2),
  
  
  lmer( AGC~ AverageDBH+
        #Density+ 
        # DBH_CV +
        #Shannon +#FDis.y+ 
        #SR+SR_Bro+SR_Other+
        #FDis+  FDisB+  
        #FDis.y_Broadleaf+
        #pd_Broadleaf+
        #shan_Broadleaf+
        #FDis.y +
        PD+
        # MAT+
          MAP+
        #Altitude+
        SoilThickness +HumusLayerThickness+ (1|RegionName)#+LitterThickness
      ,data=data_sem2)
  
  
)

summary(mod_220_Low, .progressBar = FALSE) 

AIC(mod_220_Low, AIC.type = "dsep",
    aicc = T )



fn<-paste("C_psem_mode_Low_FD.R")
save(mod_220_Low,file=fn)

mod_bor=mod_220_Low

# 
pz_bor <- lapply(mod_bor[1:(length(mod_bor) - 1)], function(x) { which(summary(x)$co[, grep('P', colnames(summary(x)$co))] > 0.05) })
newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    uFDis.yate(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "C_psem_Newmode_Low_FD.R"
save(newmod_bor, file = fn_new_bor)

#################################################################################
####Pinus massoniana####
####PSEM_All:PSEM####
rm(list=ls())

data_sem2=read.csv( "data_P.csv")
colnames(data_sem2 )


log.norm.scale <- function(df) {
  
  return(log_df <-log(df+1)  )  
  
}

for (c in c("pd_Broadleaf","shan_Broadleaf",         
            "FDis","FDisB"  )) {
  data_sem2[, c] <- log.norm.scale(data_sem2[, c]) 
}

colnames(data_sem2 )

#scale
log.norm.scale <- function(log_df) {
  
  return(scale(log_df)) 
}

for (c in c(2:18)) {
  data_sem2[, c] <- log.norm.scale(data_sem2[, c]) 
}
data_sem2[,2:18] <- apply(data_sem2[,2:18], 2, as.numeric)


data_sem2 <- as.data.frame(data_sem2)
data_sem2[is.na(data_sem2)] <- 0


mod_220_All<-psem(
  lmer(AverageDBH~ Density+
       Age+#PD+
       #DBH_CV +
      # Shannon +
         #FDisB+   
       FDis+    
      
       # pd_Broadleaf+ 
       
       #MAT+HumusLayerThickness+LitterThickness +
       MAP+
      # Altitude+
      SoilThickness+ (1|RegionName) #+
     
     ,data=data_sem2),
  
  lmer(DBH_CV~ #AverageDBH+
       #shan_Broadleaf+
       #Age+
       Density+ 
       #SR_Other+
       
       #Shannon +#FDis.y+ 
       #SR+SR_Bro+SR_Other+
       FDis+    
       #FDis.y_Broadleaf+
       #shan_Broadleaf+ #FDis.y +
       
       #MAT+MAP+
       #Altitude#+
       SoilThickness+ (1|RegionName) #+HumusLayerThickness+LitterThickness
     ,data=data_sem2),
  
  lmer(Density~ #Age+
         PD+(1|RegionName)#+ #SR_Other +
     #Shannon# +#FDis.y+ 
     #SR+SR_Bro+SR_Other+
     # FDis.y +   
     #FDis.y_Broadleaf+
     #shan_Broadleaf+ #FDis.y +
     
     #MAT+MAP+
     #Altitude+
     #SoilThickness #+HumusLayerThickness+LitterThickness
     ,data=data_sem2),
  
  lmer(FDis~ Age+#MAT+
       #Shannon+
      # shan_Broadleaf+#MAP+
       #Altitude+
       SoilThickness + (1|RegionName)#HumusLayerThickness+LitterThickness
     , data=data_sem2),
  #   lm(Shannon~# MAT+
  #        # MAP
  # Age+
  #        #Altitude+
  #        SoilThickness #HumusLayerThickness+LitterThickness
  #      , data=data_sem2),
  #   
  
  
  lmer( AGC~ Age+AverageDBH+
        Density+ 
        DBH_CV +
        
        PD+
        #FDisB+
          (1|RegionName)
      ,data=data_sem2),PD%~~%FDis
 
)

summary(mod_220_All, .progressBar = FALSE) 

AIC(mod_220_All, AIC.type = "dsep",
    aicc = T )
fn<-paste("P_psem_mode_All_FD.R")
save(mod_220_All,file=fn)

mod_bor=mod_220_All

# 
pz_bor <- lapply(mod_bor[1:(length(mod_bor) - 1)], function(x) { which(summary(x)$co[, grep('P', colnames(summary(x)$co))] > 0.05) })
newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    uFDis.yate(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "P_psem_Newmode_All_FD.R"
save(newmod_bor, file = fn_new_bor)

####PSEM_High:PSEM #########################
rm(list=ls())

data_sem2=read.csv( "data_P.csv")

C_Plot=read.csv("data_P_HLS_LCS_forests.csv")

C_Plot=C_Plot[,c("Plot","Group"  )]

data_sem2 = merge(C_Plot,data_sem2,by.x="Plot", by.y="Plot" )
colnames(data_sem2)

data_sem2 = data_sem2[data_sem2$Group == "2" , ]
colnames(data_sem2)


log.norm.scale <- function(df) {
  
  return(log_df <-log(df+1)  )  
  
}

for (c in c("pd_Broadleaf","shan_Broadleaf",         
            "FDis","FDisB"  )) {
  data_sem2[, c] <- log.norm.scale(data_sem2[, c]) 
}

colnames(data_sem2 )

#scale
log.norm.scale <- function(log_df) {
  
  return(scale(log_df)) 
}

for (c in c(3:19)) {
  data_sem2[, c] <- log.norm.scale(data_sem2[, c]) 
}
data_sem2[,3:19] <- apply(data_sem2[,3:19], 2, as.numeric)


data_sem2 <- as.data.frame(data_sem2)
data_sem2[is.na(data_sem2)] <- 0

mod_220_High <-psem(
  lmer(AverageDBH~ Density+
       
       #DBH_CV +
         PD+
       #Shannon +#FDis.y+ 
       #SR+SR_Bro+SR_Other+
       #FDis  + 
       (1|RegionName) 
     
     ,data=data_sem2),
  
  lmer(DBH_CV~ #AverageDBH+
       
       #Density+ 
       #SR_Other+
       
       #Shannon +#+FDis.y+ 
       #SR+SR_Bro+SR_Other+
       FDis + 
       #FDis.y_Broadleaf+
       #shan_Broadleaf+ #FDis.y +
       
       #MAT+MAP+
       #Altitude+
       SoilThickness + 
         (1|RegionName)#+HumusLayerThickness+LitterThickness
     ,data=data_sem2),
  
  # lmer(Density~ 
  #      
  #      #DBH_CV +
  #      #Shannon +#FDis.y+ 
  #      #SR+SR_Bro+SR_Other+
  #      FDis  +PD + 
  #      (1|RegionName)
  #    
  #    ,data=data_sem2),
  # 
 
  
  lmer(PD~ #MAT+
       
       
       #MAP+
       # Altitude+
       SoilThickness + 
       (1|RegionName)#+HumusLayerThickness+LitterThickness
     , data=data_sem2),
  
  
  lmer( AGC~ AverageDBH+
       # Density+ 
       #  DBH_CV +
        # Shannon + 
        FDis + #FDisB+   
       
        #pd_Broadleaf+
        PD+ 
          (1|RegionName)
      #MAT+
      # MAP+
      #Altitude+
      # SoilThickness #+HumusLayerThickness+LitterThickness
      ,data=data_sem2),PD%~~% FDis
  
 
 
)


summary(mod_220_High, .progressBar = FALSE) 

AIC(mod_220_High, AIC.type = "dsep",
    aicc = T )

fn<-paste("P_psem_mode_High_FD.R")
save(mod_220_High,file=fn)

mod_bor=mod_220_High

# 
pz_bor <- lapply(mod_bor[1:(length(mod_bor) - 1)], function(x) { which(summary(x)$co[, grep('P', colnames(summary(x)$co))] > 0.05) })
newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    uFDis.yate(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "P_psem_Newmode_High_FD.R"
save(newmod_bor, file = fn_new_bor)

####PSEM_Low:PSEM #########################
rm(list=ls())

data_sem2=read.csv( "data_P.csv")

C_Plot=read.csv("data_P_HLS_LCS_forests.csv")

C_Plot=C_Plot[,c("Plot","Group"  )]

data_sem2 = merge(C_Plot,data_sem2,by.x="Plot", by.y="Plot" )
colnames(data_sem2)

data_sem2 = data_sem2[data_sem2$Group == "1" , ]
colnames(data_sem2)

log.norm.scale <- function(df) {
  
  return(log_df <-log(df+1)  )  
  
}

for (c in c("pd_Broadleaf","shan_Broadleaf",         
            "FDis","FDisB"  )) {
  data_sem2[, c] <- log.norm.scale(data_sem2[, c]) 
}

colnames(data_sem2 )

#scale
log.norm.scale <- function(log_df) {
  
  return(scale(log_df)) 
}

for (c in c(3:19)) {
  data_sem2[, c] <- log.norm.scale(data_sem2[, c]) 
}
data_sem2[,3:19] <- apply(data_sem2[,3:19], 2, as.numeric)


data_sem2 <- as.data.frame(data_sem2)
data_sem2[is.na(data_sem2)] <- 0

mod_220_Low <-psem(
  lmer(AverageDBH~ #Density+
       
       DBH_CV +
       PD+ 
         (1|RegionName)
     # Shannon +FDis.y 
     #SR+SR_Bro+SR_Other+
     #FDis.y   
     #FDis.y_Broadleaf+
     #shan_Broadleaf+ #FDis.y +
     
     # MAT+
     #MAP+
     #Altitude#+
     # SoilThickness +
     #HumusLayerThickness#
     #LitterThickness 
     ,data=data_sem2),
  
  
  
  lmer( AGC~ AverageDBH+
        Density+ 
        DBH_CV+  
          (1|RegionName)
      ,data=data_sem2)
  
  
)

summary(mod_220_Low, .progressBar = FALSE) 

AIC(mod_220_Low, AIC.type = "dsep",
    aicc = T )

fn<-paste("P_psem_mode_Low_FD.R")
save(mod_220_Low,file=fn)

mod_bor=mod_220_Low

# 
pz_bor <- lapply(mod_bor[1:(length(mod_bor) - 1)], function(x) { which(summary(x)$co[, grep('P', colnames(summary(x)$co))] > 0.05) })
newmod_bor <- lapply(1:length(pz_bor), function(i) {
  dropcoeff <- vector()
  if (length(pz_bor[[i]][!names(pz_bor[[i]]) %in% '(Intercept)']) > 0) {
    dropcoeff <- names(pz_bor[[i]])[!names(pz_bor[[i]]) %in% '(Intercept)']
    uFDis.yate(mod_bor[[i]], as.formula(paste('. ~ . -', paste(dropcoeff, collapse = '-'))), data = mod_bor[[i]]@frame)
  } else {
    mod_bor[[i]]
  }
})
newmod_bor
fn_new_bor <- "P_psem_Newmode_Low_FD.R"
save(newmod_bor, file = fn_new_bor)

#################################################################################
#################################################################################