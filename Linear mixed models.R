##===============================================================================
## PD:LMM ---------
##===============================================================================
rm(list=ls())

readDir=c("E:/MYData/Compare_P_C/result")
setwd(readDir)
library(lmerTest)
library(piecewiseSEM)
library(data.table)
library(pROC)
library(plyr)
library(dplyr)
library(doParallel)
library(foreign)
library(car)
library(MuMIn)
####Cunninghamia lanceolata PD：LMM_All####
rm(list=ls())

data_sem2=read.csv( "E:/MYData/Compare_P_C/投稿/Ecological Applications/code/data_C.csv")
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

modall= lmer( AGC~Age+ AverageDBH+
              Density+ 
              DBH_CV +
              Shannon + 
              PD +    
             # pd_Broadleaf+
              shan_Broadleaf+
              MAT+
              MAP+
              Altitude+
              SoilThickness +HumusLayerThickness+LitterThickness+(1|RegionName)
            ,data=data_sem2, REML = FALSE)


summary(modall)
r.squaredGLMM(modall)
round(vif(modall), 2)
####Model Selection####
ddall<-dredge(modall,options(na.action = "na.fail"))
subset(ddall,delta<2)
deall<-model.avg(ddall, subset = delta < 2)
summary(deall) 

####Extract model data####
library(broom)

model_avg <- summary(deall) 

coef_full_avg <- model_avg$coefmat.full

coef_table_full <- as.data.frame(coef_full_avg)


colnames(coef_table_full) <- c("Estimate", "Std.Error","Adjusted.SE", "z.value","Pr(>|z|)")
coef_table_full


write.csv(coef_table_full, "PD_full_average_All_310.csv")
#################################################################################
rm(list=ls())
####Cunninghamia lanceolata PD：LMM_High####
rm(list=ls())

data_sem2=read.csv( "E:/MYData/Compare_P_C/投稿/Ecological Applications/code/data_C.csv")
colnames(data_sem2 )

C_Plot= read.csv("E:/MYData/Compare_P_C/投稿/Ecological Applications/code/data_C_HLS_LCS_forests.csv")
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


modall= lmer( AGC~ AverageDBH+
              Density+ 
              DBH_CV +
              Shannon + 
              PD +    
              
              #pd_Broadleaf+
              shan_Broadleaf+
              
              MAT+
              MAP+
              Altitude+
              SoilThickness +HumusLayerThickness+LitterThickness+(1|RegionName)
            ,data=data_sem2, REML = FALSE)


summary(modall)
r.squaredGLMM(modall)
vif(modall) 
round(vif(modall), 2)
####Model Selection####
library(MuMIn)
ddall<-dredge(modall,options(na.action = "na.fail"))
subset(ddall,delta<2)
deall<-model.avg(ddall, subset = delta < 2)
summary(deall) 

####Extract model data####
library(broom)

model_avg <- summary(deall) 


coef_full_avg <- model_avg$coefmat.full


coef_table_full <- as.data.frame(coef_full_avg)


colnames(coef_table_full) <- c("Estimate", "Std.Error","Adjusted.SE", "z.value","Pr(>|z|)")
coef_table_full


write.csv(coef_table_full, "PD_full_average_High_310.csv")

#################################################################################
rm(list=ls())
####Cunninghamia lanceolata PD：LMM_Low####
rm(list=ls())

data_sem2=read.csv( "E:/MYData/Compare_P_C/投稿/Ecological Applications/code/data_C.csv")
colnames(data_sem2 )

C_Plot= read.csv("E:/MYData/Compare_P_C/投稿/Ecological Applications/code/data_C_HLS_LCS_forests.csv")
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





modall= lmer( AGC~ AverageDBH+
              Density+ 
              DBH_CV +
              Shannon + 
             
              
              PD +    
              
              #pd_Broadleaf+
              shan_Broadleaf+
              
              MAT+
              MAP+
              Altitude+
              SoilThickness +HumusLayerThickness+LitterThickness+(1|RegionName)
            ,data=data_sem2, REML = FALSE)
round(vif(modall), 2)

summary(modall)
r.squaredGLMM(modall)
vif(modall) 
####Model Selection####
library(MuMIn)
ddall<-dredge(modall,options(na.action = "na.fail"))
subset(ddall,delta<2)
deall<-model.avg(ddall, subset = delta < 2)
summary(deall) 

####Extract model data####
library(broom)

model_avg <- summary(deall) 


coef_full_avg <- model_avg$coefmat.full


coef_table_full <- as.data.frame(coef_full_avg)


colnames(coef_table_full) <- c("Estimate", "Std.Error","Adjusted.SE", "z.value","Pr(>|z|)")
coef_table_full


write.csv(coef_table_full, "PD_full_average_Low_310.csv")


#################################################################################
#################################################################################
####Pinus massoniana PD：LMM_All####
rm(list=ls())

data_sem2=read.csv( "E:/MYData/Compare_P_C/投稿/Ecological Applications/code/data_P.csv")
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



modall=  lmer( AGC~ Age+
               AverageDBH+
               Density+ 
               DBH_CV +
               Shannon + 
              
              
               PD +    
              # pd_Broadleaf+
               shan_Broadleaf+
               
               MAT+
               MAP+
               Altitude+
               SoilThickness +HumusLayerThickness+LitterThickness+ (1|RegionName)
             ,data=data_sem2, REML = FALSE)
round(vif(modall), 2)

summary(modall)
r.squaredGLMM(modall)
vif(modall) 
####Model Selection####
library(MuMIn)
ddall<-dredge(modall,options(na.action = "na.fail"))
subset(ddall,delta<2)
deall<-model.avg(ddall, subset = delta < 2)
summary(deall) 

####Extract model data####
library(broom)

model_avg <- summary(deall) 


coef_full_avg <- model_avg$coefmat.full


coef_table_full <- as.data.frame(coef_full_avg)


colnames(coef_table_full) <- c("Estimate", "Std.Error","Adjusted.SE", "z.value","Pr(>|z|)")
coef_table_full


write.csv(coef_table_full, "PD_full_average_All_220.csv")

#################################################################################
####Pinus massoniana PD：LMM_High####
rm(list=ls())

data_sem2=read.csv( "E:/MYData/Compare_P_C/投稿/Ecological Applications/code/data_P.csv")

C_Plot=read.csv("E:/MYData/Compare_P_C/投稿/Ecological Applications/code/data_P_HLS_LCS_forests.csv")

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
  
  return((log_df - min(log_df)) / (max(log_df) - min(log_df)))   
}

for (c in c(3:19)) {
  data_sem2[, c] <- log.norm.scale(data_sem2[, c]) 
}
data_sem2[,3:19] <- apply(data_sem2[,3:19], 2, as.numeric)


data_sem2 <- as.data.frame(data_sem2)
data_sem2[is.na(data_sem2)] <- 0


modall=  lmer( AGC~ 
               AverageDBH+
               Density+ 
               DBH_CV +
               Shannon + 
               
               PD +    
               #pd_Broadleaf+
               shan_Broadleaf+ 
               
               MAT+
               MAP+
               Altitude+
               SoilThickness +HumusLayerThickness+LitterThickness+ (1|RegionName)
             ,data=data_sem2, REML = FALSE)
round(vif(modall), 2)

summary(modall)
r.squaredGLMM(modall)
vif(modall) 
####Model Selection####
library(MuMIn)
ddall<-dredge(modall,options(na.action = "na.fail"))
subset(ddall,delta<2)
deall<-model.avg(ddall, subset = delta < 2)
summary(deall) 

####Extract model data####
library(broom)

model_avg <- summary(deall) 


coef_full_avg <- model_avg$coefmat.full


coef_table_full <- as.data.frame(coef_full_avg)


colnames(coef_table_full) <- c("Estimate", "Std.Error","Adjusted.SE", "z.value","Pr(>|z|)")
coef_table_full


write.csv(coef_table_full, "PD_full_average_High_220.csv")

#################################################################################

####Pinus massoniana PD：LMM_Low####
rm(list=ls())

data_sem2=read.csv( "E:/MYData/Compare_P_C/投稿/Ecological Applications/code/data_P.csv")

C_Plot=read.csv("E:/MYData/Compare_P_C/投稿/Ecological Applications/code/data_P_HLS_LCS_forests.csv")

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


modall=  lmer( AGC~ 
               AverageDBH+
               Density+ 
               DBH_CV +
               Shannon + 
               # FDis.y+ 
               #SR+
               #SR_Bro+
               #SR_Other+
               PD +    
               #pd_Broadleaf+
               shan_Broadleaf+ #pd +
               
               MAT+
               MAP+
               Altitude+
               SoilThickness +HumusLayerThickness+LitterThickness+ (1|RegionName)
             ,data=data_sem2, REML = FALSE)
round(vif(modall), 2)

summary(modall)
r.squaredGLMM(modall)
vif(modall) 
####Model Selection####
library(MuMIn)
ddall<-dredge(modall,options(na.action = "na.fail"))
subset(ddall,delta<2)
deall<-model.avg(ddall, subset = delta < 2)
summary(deall) 

####Extract model data####
library(broom)

model_avg <- summary(deall) 


coef_full_avg <- model_avg$coefmat.full


coef_table_full <- as.data.frame(coef_full_avg)


colnames(coef_table_full) <- c("Estimate", "Std.Error","Adjusted.SE", "z.value","Pr(>|z|)")
coef_table_full


write.csv(coef_table_full, "PD_full_average_Low_220.csv")

#################################################################################
##===============================================================================
## FD:LMM ---------
##===============================================================================
rm(list=ls())
####Cunninghamia lanceolata FD：LMM_All####
rm(list=ls())

data_sem2=read.csv( "E:/MYData/Compare_P_C/投稿/Ecological Applications/code/data_C.csv")
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



modall= lmer( AGC~Age+ AverageDBH+
              Density+ 
              DBH_CV +
              #Shannon + 
              FDis+ FDisB+
              
              PD +    
              
             # pd_Broadleaf+
             # shan_Broadleaf+
              
              MAT+
              MAP+
              Altitude+
              SoilThickness +HumusLayerThickness+LitterThickness+ (1|RegionName)
            ,data=data_sem2, REML = FALSE)
round(vif(modall), 2)
summary(modall)
r.squaredGLMM(modall)
vif(modall) 
####Model Selection####
library(MuMIn)
ddall<-dredge(modall,options(na.action = "na.fail"))
subset(ddall,delta<2)
deall<-model.avg(ddall, subset = delta < 2)
summary(deall) 

####Extract model data####
library(broom)

model_avg <- summary(deall) 


coef_full_avg <- model_avg$coefmat.full


coef_table_full <- as.data.frame(coef_full_avg)


colnames(coef_table_full) <- c("Estimate", "Std.Error","Adjusted.SE", "z.value","Pr(>|z|)")
coef_table_full


write.csv(coef_table_full, "FD_full_average_All_310.csv")




#################################################################################
rm(list=ls())
####Cunninghamia lanceolata FD：LMM_High####
rm(list=ls())

data_sem2=read.csv( "E:/MYData/Compare_P_C/投稿/Ecological Applications/code/data_C.csv")
colnames(data_sem2 )

C_Plot= read.csv("E:/MYData/Compare_P_C/投稿/Ecological Applications/code/data_C_HLS_LCS_forests.csv")
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

modall= lmer( AGC~ AverageDBH+
              Density+ 
              DBH_CV +
              #Shannon + 
              FDis+ FDisB+
              
              PD +    
              
              #pd_Broadleaf+
              #shan_Broadleaf+
              
              MAT+
              MAP+
              Altitude+
              SoilThickness +HumusLayerThickness+LitterThickness+(1|RegionName)
            ,data=data_sem2, REML = FALSE)

round(vif(modall), 2)
summary(modall)
r.squaredGLMM(modall)
vif(modall) 
####Model Selection####
library(MuMIn)
ddall<-dredge(modall,options(na.action = "na.fail"))
subset(ddall,delta<2)
deall<-model.avg(ddall, subset = delta < 2)
summary(deall) 

####Extract model data####
library(broom)

model_avg <- summary(deall) 


coef_full_avg <- model_avg$coefmat.full


coef_table_full <- as.data.frame(coef_full_avg)


colnames(coef_table_full) <- c("Estimate", "Std.Error","Adjusted.SE", "z.value","Pr(>|z|)")
coef_table_full


write.csv(coef_table_full, "FD_full_average_High_310.csv")

#################################################################################
rm(list=ls())
####Cunninghamia lanceolata FD：LMM_Low####
rm(list=ls())

data_sem2=read.csv( "E:/MYData/Compare_P_C/投稿/Ecological Applications/code/data_C.csv")
colnames(data_sem2 )

C_Plot= read.csv("E:/MYData/Compare_P_C/投稿/Ecological Applications/code/data_C_HLS_LCS_forests.csv")
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

modall= lmer( AGC~ AverageDBH+
              Density+ 
              DBH_CV +
              #Shannon + 
              FDis+ FDisB+
              
              PD +    
              
              #pd_Broadleaf+
             # shan_Broadleaf+
              
              MAT+
              MAP+
              Altitude+
              SoilThickness +HumusLayerThickness+LitterThickness+(1|RegionName)
            ,data=data_sem2, REML = FALSE)
round(vif(modall), 2)

summary(modall)
r.squaredGLMM(modall)
vif(modall) 
####Model Selection####
library(MuMIn)
ddall<-dredge(modall,options(na.action = "na.fail"))
subset(ddall,delta<2)
deall<-model.avg(ddall, subset = delta < 2)
summary(deall) 

####Extract model data####
library(broom)

model_avg <- summary(deall) 


coef_full_avg <- model_avg$coefmat.full


coef_table_full <- as.data.frame(coef_full_avg)


colnames(coef_table_full) <- c("Estimate", "Std.Error","Adjusted.SE", "z.value","Pr(>|z|)")
coef_table_full


write.csv(coef_table_full, "FD_full_average_Low_310.csv")


#################################################################################
#################################################################################
####Pinus massoniana FD：LMM_All####
rm(list=ls())

data_sem2=read.csv( "E:/MYData/Compare_P_C/投稿/Ecological Applications/code/data_P.csv")
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



modall=  lmer( AGC~ Age+
               AverageDBH+
               Density+ 
               DBH_CV +
              
               FDis+ FDisB+ 
             
               PD +    
              # pd_Broadleaf+
               #shan_Broadleaf+ #pd +
               
               MAT+
               MAP+
               Altitude+
               SoilThickness +HumusLayerThickness+LitterThickness+(1|RegionName)
             ,data=data_sem2, REML = FALSE)
round(vif(modall), 2)

summary(modall)
r.squaredGLMM(modall)
vif(modall) 
####Model Selection####
library(MuMIn)
ddall<-dredge(modall,options(na.action = "na.fail"))
subset(ddall,delta<2)
deall<-model.avg(ddall, subset = delta < 2)
summary(deall) 

####Extract model data####
library(broom)

model_avg <- summary(deall) 


coef_full_avg <- model_avg$coefmat.full


coef_table_full <- as.data.frame(coef_full_avg)


colnames(coef_table_full) <- c("Estimate", "Std.Error","Adjusted.SE", "z.value","Pr(>|z|)")
coef_table_full


write.csv(coef_table_full, "FD_full_average_All_220.csv")

#################################################################################










####Pinus massoniana FD：LMM_High####
rm(list=ls())

data_sem2=read.csv( "E:/MYData/Compare_P_C/投稿/Ecological Applications/code/data_P.csv")

C_Plot=read.csv("E:/MYData/Compare_P_C/投稿/Ecological Applications/code/data_P_HLS_LCS_forests.csv")

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
  
  return((log_df - min(log_df)) / (max(log_df) - min(log_df)))  
}

for (c in c(3:19)) {
  data_sem2[, c] <- log.norm.scale(data_sem2[, c]) 
}
data_sem2[,3:19] <- apply(data_sem2[,3:19], 2, as.numeric)


data_sem2 <- as.data.frame(data_sem2)
data_sem2[is.na(data_sem2)] <- 0

modall=  lmer( AGC~ 
               AverageDBH+
               Density+ 
               DBH_CV +
               
               
               FDis+ FDisB+
               
               PD +    
               #pd_Broadleaf+
               
               
               MAT+
               MAT+
               Altitude+
               SoilThickness +HumusLayerThickness+LitterThickness+(1|RegionName) 
             ,data=data_sem2, REML = FALSE )
library(car)
round(vif(modall), 2)
summary(modall)
r.squaredGLMM(modall)
vif(modall) 
####Model Selection####
library(MuMIn)
ddall<-dredge(modall,options(na.action = "na.fail"))
subset(ddall,delta<2)
deall<-model.avg(ddall, subset = delta < 2)
summary(deall) 

####Extract model data####
library(broom)

model_avg <- summary(deall) 


coef_full_avg <- model_avg$coefmat.full


coef_table_full <- as.data.frame(coef_full_avg)


colnames(coef_table_full) <- c("Estimate", "Std.Error","Adjusted.SE", "z.value","Pr(>|z|)")
coef_table_full


write.csv(coef_table_full, "FD_full_average_High_220.csv")

#################################################################################










####Pinus massoniana FD：LMM_Low####
rm(list=ls())

data_sem2=read.csv( "E:/MYData/Compare_P_C/投稿/Ecological Applications/code/data_P.csv")

C_Plot=read.csv("E:/MYData/Compare_P_C/投稿/Ecological Applications/code/data_P_HLS_LCS_forests.csv")

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
modall=  lmer( AGC~ 
               AverageDBH+
               Density+ 
               DBH_CV +
               # Shannon + 
               FDis+ FDisB+
               #SR+
               #SR_Bro+
               #SR_Other+
               PD +    
               #pd_Broadleaf+
             
               MAT+
               MAT+
               Altitude+
               SoilThickness +HumusLayerThickness+LitterThickness+(1|RegionName)
             ,data=data_sem2, REML = FALSE)
round(vif(modall), 2)

summary(modall)
r.squaredGLMM(modall)
vif(modall) 
####Model Selection####
library(MuMIn)
ddall<-dredge(modall,options(na.action = "na.fail"))
subset(ddall,delta<2)
deall<-model.avg(ddall, subset = delta < 2)
summary(deall) 

####Extract model data####
library(broom)

model_avg <- summary(deall) 


coef_full_avg <- model_avg$coefmat.full


coef_table_full <- as.data.frame(coef_full_avg)


colnames(coef_table_full) <- c("Estimate", "Std.Error","Adjusted.SE", "z.value","Pr(>|z|)")
coef_table_full


write.csv(coef_table_full, "FD_full_average_Low_220.csv")

#################################################################################
##===============================================================================
## FD:LMMPlot ---------
##===============================================================================
#######Cunninghamia lanceolata：plot_all######
rm(list=ls())
library(tidyverse)
library(ggplot2)
library(ggh4x)
library(dplyr)

df_All<-read.csv("FD_full_average_All_310.csv",check.names = F)
df_High<-read.csv("FD_full_average_High_310.csv",check.names = F)
df_Low<-read.csv("FD_full_average_Low_310.csv",check.names = F)

df_All$Grid<-"General"
df_High$Grid<-"High"
df_Low$Grid<-"Low"

df=rbind( df_All,df_High,df_Low )
colnames(df) <- c("var" ,"Estimate", "Std. Error","Adjusted.SE", "z.value","Pr(>|z|)", "Grid"   )

unique(df$var)

df <- df %>%
  filter(!( var %in% c("(Intercept)") )  ) %>%
  dplyr::mutate(group = case_when(
    var %in% c("Age","AverageDBH" ,"DBH_CV"  ,"Density") ~ "Stand structural",
    var %in% c( "PD" ,"pd_Broadleaf"  ,"FDis" ,"FDisB") ~ "Diversity",
    var %in% c("Altitude",
               "MAT", "MAP" , "SoilThickness" ,"HumusLayerThickness", "LitterThickness") ~ "Geographic Location",
    TRUE ~ NA_character_  # In case there are any terms not listed above
  ))



cols<-c("#ADD9E6" , "#d2de61" ,"#dd7500")


####To calculate the percentage for each category########

B=df %>% 
  group_by(Grid,var
  ) %>%
  dplyr::mutate(groups=fct_relevel(var,
                                   rev(c()))) %>% 
  group_by(Grid,
           #group,
           var
  ) %>% 
  dplyr::summarise(sum_value = sum(abs(Estimate), na.rm = TRUE)) %>% 
  dplyr::mutate(new_col = sum_value / sum(sum_value, na.rm = TRUE))  %>% 
  dplyr::mutate(new_col = round(new_col, 4))  
B
####Plot1########

df <- df %>%
  dplyr::mutate(
    sign = if_else(
      Grid == "Low", Estimate +2*abs(`Std. Error`) ,
      if_else(Grid == "High", Estimate + 2*abs( `Std. Error`),
              if_else(Grid == "General", Estimate +2*abs( `Std. Error`), NA_real_))
    )
  )
custom_order <- c("AverageDBH" ,"DBH_CV"  ,"Density",
                  "FDisB", "FDis"  ,#"pd_Broadleaf"  ,
                  "PD" ,"Age",
                  "Altitude",
                  "MAT", "MAP" , "SoilThickness" ,"HumusLayerThickness", "LitterThickness"  )

unique(df$var)
p_all310= df %>%
  dplyr::mutate(Grid = fct_relevel(Grid, c("General", "Low","High")))  %>% 
  dplyr::mutate(group1=fct_relevel(var,custom_order))  %>% 
  group_by(group1, Grid) %>% 
  arrange(Estimate) %>% 
  dplyr::mutate(var = fct_relevel(var, var)) %>% 
  dplyr::mutate(signi = case_when(
    `Pr(>|z|)` < 0.1 & `Pr(>|z|)` >= 0.05 ~ '',
    `Pr(>|z|)` > 0.1 ~ '',
    `Pr(>|z|)` > 0.05 ~ '',
    `Pr(>|z|)` < 0.05 & `Pr(>|z|)` >= 0.01 ~ '*',
    `Pr(>|z|)` < 0.01 & `Pr(>|z|)` >= 0.001 ~ '**',
    `Pr(>|z|)` < 0.001 ~ '***'
  )) %>% 
  
  arrange(group) %>%  
  ggplot(aes(x = Estimate, y = fct_reorder(var, desc(group1)), color = Grid, group = Grid)) +
  geom_vline(xintercept=c(0), linetype="longdash")+
  geom_linerange(
    aes(xmin = Estimate - `Std. Error`- `Std. Error`, xmax = Estimate +  `Std. Error`+ `Std. Error`, color = Grid, group = Grid),
    size =0.1, 
    color = "#040505",
    alpha =1, position = position_dodge(width = 0.6)
  ) +
  geom_point(
    aes(color = Grid, group = Grid),
    show.legend = F,
    size = 3, alpha = 1, 
    #shape=21,
    #color = "black",
    #fill = "black",
    position = position_dodge(width = 0.6)
  )+
  scale_shape_manual(values = c(21,21,21,21,21)) +
  
  #scale_color_manual(values = c("General"="#a1a1a4", "Low"="#B7D982","High"="#1572CD" )) +
  #scale_fill_manual(values =  c("General"="#a1a1a4", "Low"="#B7D982","High"="#1572CD" )) +
  

  geom_text(aes(x = sign , y = var,label = signi),
            # size =3, family = "serif",  color ="black",fontface = "bold", vjust = -0.3,position = position_dodge(width = 0.6), hjust = 0.5) +
            size =3, family = "serif",  color ="black",
            fontface = "bold", position = position_dodge(width = 0.6), hjust = 0) +
  
  geom_hline(yintercept = seq(0.5, length(custom_order) - 0.5), linetype = "solid", color = "#f3f4f7") +  # 添加横线
  
  theme_bw( )+
  
  labs(title = expression((a)~italic(Cunninghamia~lanceolata)))+
  theme(plot.title = element_text(size=16,family="serif"))+
  
  labs(x="Parameter Estimate",y=NULL)+
  theme(axis.title.x = element_text(size = 14, family =  "serif"))+
  #labs(x=NULL,y=NULL)+
  
  #theme(axis.title.x = element_text(size = 12, family =  "serif"))+
  
  scale_x_continuous(limits = c(-0.9,1),
                     breaks = c(-0.5,0,0.5,1))+
  
  
  scale_y_discrete(
    labels = c(
      
      expression(Average~DBH),
      expression(DBH~variation),
      expression(Stand~density),
      
      expression(FD[Broadleaves]),
      expression(FD),
      #expression(PD[Broadleaf]),
      expression(PD),
      expression(Stand~age),
      
      #expression(Species~diversity),
      
      expression(Elevation),
      expression(MAT),
      expression(MAP),
      expression(ST),
      expression(HLT),
      expression(LT)
    ),
    limits = rev#,position = "right"
  )+
  theme(#axis.text.x = element_blank(),
    axis.ticks = element_line(),
    axis.text.x = element_text(size = 12,color="black",family="serif"),
    axis.text.y = element_text(size = 14,color="black",face="bold",family="serif"),
    #axis.text.y =  element_blank(),
    axis.ticks.y = element_line()
  )+
  theme(
    panel.grid = element_blank(),
    # panel.background = element_rect(fill = "#ffffff"),  
    # panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # 设置面板边框线
    
    # panel.grid.major = element_line(color = "#f3f4f7" ),  
    #panel.grid.minor = element_line(color = "#f3f4f7"),  # 设置次要网格线
    
    axis.line = element_line(color = "black", size = 0.5),  # 设置坐标轴线
    
    axis.title = element_text(color = "black"),  # 设置坐标轴标题颜色
    axis.text = element_text(color = "black"),  # 设置坐标轴文本颜色
    axis.line.x = element_line(color = "black", size = 0.5),  # 设置x轴线颜色和大小
    axis.line.y = element_line(color = "black", size = 0.5) , # 设置y轴线颜色和大小
    #axis.line.y = element_blank() , # 设置y轴线颜色和大小
    legend.position = ""
  )+
  
  #guides(fill = guide_legend(reverse = TRUE))+
  # theme(legend.position = "right", legend.title = element_text(size = 20), 
  #        legend.text = element_text(size = 12,family="serif"#,face = "italic"
  #                                   ))+
  
  scale_color_manual(values = cols)


p_all310


################
#######Pinus massoniana：plot_all######
#rm(list=ls())
library(tidyverse)
library(ggplot2)
library(ggh4x)
library(dplyr)

df_All<-read.csv("FD_full_average_All_220.csv",check.names = F)
df_High<-read.csv("FD_full_average_High_220.csv",check.names = F)
df_Low<-read.csv("FD_full_average_Low_220.csv",check.names = F)

df_All$Grid<-"General"
df_High$Grid<-"High"
df_Low$Grid<-"Low"

df=rbind( df_All,df_High,df_Low )
colnames(df) <- c("var" ,"Estimate", "Std. Error","Adjusted.SE", "z.value","Pr(>|z|)", "Grid"   )

df <- df %>%
  filter(!( var %in% c("(Intercept)") )  ) %>%
  dplyr::mutate(group = case_when(
    var %in% c("Age","AverageDBH"  ,"DBH_CV"  ,"Density") ~ "Stand structural",
    var %in% c(  "PD" ,"pd_Broadleaf"  ,"FDis" ,"FDisB") ~ "Diversity",
    var %in% c("Altitude",
               "MAT" , "MAP"  , "SoilThickness" ,"HumusLayerThickness", "LitterThickness") ~ "Geographic Location",
    TRUE ~ NA_character_  # In case there are any terms not listed above
  ))

unique(df$var)

cols<-c("#ADD9E6" , "#d2de61" ,"#dd7500")


####To calculate the percentage for each category########

B=df %>% 
  group_by(Grid,var
  ) %>%
  dplyr::mutate(groups=fct_relevel(var,
                                   rev(c()))) %>% 
  group_by(Grid,
           #group,
           var
  ) %>% 
  dplyr::summarise(sum_value = sum(abs(Estimate), na.rm = TRUE)) %>% 
  dplyr::mutate(new_col = sum_value / sum(sum_value, na.rm = TRUE))  %>% 
  dplyr::mutate(new_col = round(new_col, 4))  
B
####Plot1########

df <- df %>%
  dplyr::mutate(
    sign = if_else(
      Grid == "Low", Estimate +2*abs(`Std. Error`) ,
      if_else(Grid == "High", Estimate + 2*abs( `Std. Error`),
              if_else(Grid == "General", Estimate + 2*abs( `Std. Error`), NA_real_))
    )
  )
custom_order <- c("AverageDBH"  ,"DBH_CV"  ,"Density","FDisB",
                  "FDis"  ,#"pd_Broadleaf"  ,
                  "PD" ,"Age",
                  "Altitude",
                  "MAT" , "MAP"  , "SoilThickness" ,"HumusLayerThickness", "LitterThickness"  )

unique(df$var)
p_all220= df %>%
  dplyr::mutate(Grid = fct_relevel(Grid, c("General", "Low","High")))  %>% 
  dplyr::mutate(group1=fct_relevel(var,custom_order))  %>% 
  group_by(group1, Grid) %>% 
  arrange(Estimate) %>% 
  dplyr::mutate(var = fct_relevel(var, var)) %>% 
  dplyr::mutate(signi = case_when(
    `Pr(>|z|)` < 0.1 & `Pr(>|z|)` >= 0.05 ~ '',
    `Pr(>|z|)` > 0.1 ~ '',
    `Pr(>|z|)` > 0.05 ~ '',
    `Pr(>|z|)` < 0.05 & `Pr(>|z|)` >= 0.01 ~ '*',
    `Pr(>|z|)` < 0.01 & `Pr(>|z|)` >= 0.001 ~ '**',
    `Pr(>|z|)` < 0.001 ~ '***'
  )) %>% 
  
  arrange(group) %>%  
  ggplot(aes(x = Estimate, y = fct_reorder(var, desc(group1)), color = Grid, group = Grid)) +
  geom_vline(xintercept=c(0), linetype="longdash")+
  geom_linerange(
    aes(xmin = Estimate - `Std. Error`- `Std. Error`, xmax = Estimate + `Std. Error`+ `Std. Error`, color = Grid, group = Grid),
    size =0.1, 
    color = "#040505",
    alpha =1, position = position_dodge(width = 0.6)
  ) +
  geom_point(
    aes(color = Grid, group = Grid),
    show.legend = F,
    size = 3, alpha = 1, 
    #shape=21,
    #color = "black",
    #fill = "black",
    position = position_dodge(width = 0.6)
  )+
  scale_shape_manual(values = c(21,21,21,21,21)) +
  
  #scale_color_manual(values = c("General"="#a1a1a4", "Low"="#B7D982","High"="#1572CD" )) +
  #scale_fill_manual(values =  c("General"="#a1a1a4", "Low"="#B7D982","High"="#1572CD" )) +
  

  geom_text(aes(x = sign , y = var,label = signi),
            # size =3, family = "serif",  color ="black",fontface = "bold", vjust = -0.3,position = position_dodge(width = 0.6), hjust = 0.5) +
            size =3, family = "serif",  color ="black",
            fontface = "bold", position = position_dodge(width = 0.6), hjust = 0) +
  
  geom_hline(yintercept = seq(0.5, length(custom_order) - 0.5), linetype = "solid", color = "#f3f4f7") +  # 添加横线
  
  theme_bw( )+
  
  labs(title = expression((b)~italic(Pinus~massoniana)))+
  theme(plot.title = element_text(size=16,family="serif"))+
  
  labs(x="Parameter Estimate",y=NULL)+
  theme(axis.title.x = element_text(size = 14, family =  "serif"))+
  #labs(x=NULL,y=NULL)+
  
  #theme(axis.title.x = element_text(size = 12, family =  "serif"))+
  
  scale_x_continuous(limits = c(-0.8,1),
                     breaks = c(-0.5,0,0.5,1))+
  
  
  scale_y_discrete(
    labels = c(
     
      expression(Average~DBH),
      expression(DBH~variation),
      expression(Stand~density),
      expression(FD[Broadleaves]),
      
      expression(FD),
      
      #expression(PD[Broadleaf]),
      expression(PD), expression(Stand~age),
      #expression(Shannon[Broadleaf]),
      #expression(Species~diversity),
      
      expression(Elevation),
      expression(MAT),
      expression(MAP),
      expression(ST),
      expression(HLT),
      expression(LT)
    ),
    limits = rev#,position = "right"
  )+
  theme(#axis.text.x = element_blank(),
    axis.ticks = element_line(),
    axis.text.x = element_text(size = 12,color="black",family="serif"),
    axis.text.y = element_text(size = 14,color="black",face="bold",family="serif"),
    #axis.text.y =  element_blank(),
    axis.ticks.y = element_line()
  )+
  theme(
    panel.grid = element_blank(),
    # panel.background = element_rect(fill = "#ffffff"),  
    # panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # 设置面板边框线
    
    # panel.grid.major = element_line(color = "#f3f4f7" ),  
    #panel.grid.minor = element_line(color = "#f3f4f7"),  # 设置次要网格线
    
    axis.line = element_line(color = "black", size = 0.5),  # 设置坐标轴线
    
    axis.title = element_text(color = "black"),  # 设置坐标轴标题颜色
    axis.text = element_text(color = "black"),  # 设置坐标轴文本颜色
    axis.line.x = element_line(color = "black", size = 0.5),  # 设置x轴线颜色和大小
    axis.line.y = element_line(color = "black", size = 0.5) , # 设置y轴线颜色和大小
    #axis.line.y = element_blank() , # 设置y轴线颜色和大小
    legend.position = ""
  )+
  
  #guides(fill = guide_legend(reverse = TRUE))+
  # theme(legend.position = "right", legend.title = element_text(size = 20), 
  #        legend.text = element_text(size = 12,family="serif"#,face = "italic"
  #                                   ))+
  
  scale_color_manual(values =cols)


p_all220


####ggarrange####
library(patchwork)
library(ggpubr)
Plot_Liear1=ggarrange(p_all310,p_all220,
                      ncol = 1,nrow =2,
                      heights = c(1,1),
                      widths = c(1,1),common.legend = FALSE, align = c("hv"))

Plot_Liear1
ggsave("Plot_Liear_FD250328.png", path = "E:/MYData/Compare_P_C/figure/Plot_Liear",width =5, height =8,dpi=600, plot=Plot_Liear1)


#################################################################################
#################################################################################
####Cunninghamia lanceolata：Plot2：100%Plot############
df<-read.csv("FD_full_average_All_310.csv",check.names = F)

colnames(df) <- c("var" ,"Estimate", "Std. Error","Adjusted.SE", "z.value","Pr(>|z|)"  )

unique(df$var)

df <- df %>%
  dplyr::mutate(across(where(is.numeric), round, 2))


df <- df %>%
  filter(!( var %in% c("(Intercept)") )  ) %>%
  dplyr::mutate(group = case_when(
    var %in% c("Age","AverageDBH" ,"DBH_CV"  ,"Density") ~ "Stand structural",
    var %in% c( "PD" ,"pd_Broadleaf"  ,"FDis" ,"FDisB") ~ "Diversity",
    var %in% c("Altitude",
               "MAT", "MAP" , "SoilThickness" ,"HumusLayerThickness", "LitterThickness") ~ "Geographic Location",
    TRUE ~ NA_character_  # In case there are any terms not listed above
  ))

#cols1=c("#DDE9F7", "#B5D3E9", "#82BADA","#539DCD","#004b79","#436f8e"  "#7E9E8D")
cols1=c(
  "#9AB560" ,
  "#463D09" ,
  #"#D7D1C2",
  "#fff9ea"
)
df1=df

p_all2_G=df1 %>% 
  dplyr::mutate(group=fct_relevel(group,
                           rev(c("Stand structural", "Diversity",
                                 #"CWM",
                                 "Geographic Location"
                           )))) %>% 
  group_by(group) %>% 
  dplyr::summarise(sum_value=sum(abs(Estimate))) %>% 
   dplyr::mutate(
    new_col=sum_value/sum(sum_value)*100) %>% 
  ggplot(aes(x = 1, y = new_col, label = group)) +
  geom_col(aes(fill = group, color = group), alpha =1, show.legend = FALSE) +
  scale_fill_manual(values = cols1) +
  scale_color_manual(values =cols1    ) +
  #scale_color_manual(values = c("#C0BDAB", "#C0BDAB","#C0BDAB","#C0BDAB"  )    ) + 
  coord_flip()+
  
  scale_fill_manual(values = cols1)+
  scale_y_continuous(expand = c(0,0), limits = c(0,100), breaks = c(0,20,40,60,80,100))+
  
  
  theme_minimal()+
  theme(
    axis.text.x = element_text(size = 12, family = "serif"),
    axis.title.y = element_blank(),  
    axis.line.y = element_blank(),   
    axis.line.x = element_line(),
    axis.ticks.x = element_line(),
    axis.text.y = element_blank(),
    plot.title = element_blank(),
    axis.title.x = element_text(size = 14, 
                                family =  "serif")
  ) +
  
  
  labs(y="Relative effect of estimates (%)")+
  
  labs(subtitle  = expression(Adj.~ italic(R)^2==61~ '%' ))+
  #labs(x="Pinus massoniana")+
  labs(title = expression(General))+
  
  #ggtitle(" ") +
  
  theme(plot.title = element_text(hjust = 0, face = "bold"))+
  
  theme(plot.subtitle = element_text(size = 14, family = "serif") 
  )+
  theme(plot.title = element_text(size = 14, 
                                  family =  "serif"))
p_all2_G
#################
df2<-read.csv("FD_full_average_Low_310.csv",check.names = F)

colnames(df2) <- c("var" ,"Estimate", "Std. Error","Adjusted.SE", "z.value","Pr(>|z|)"  )

unique(df2$var)

df2 <- df2 %>%
  dplyr::mutate(across(where(is.numeric), round, 2))


df2 <- df2 %>%
  filter(!( var %in% c("(Intercept)") )  ) %>%
  dplyr::mutate(group = case_when(
    var %in% c("Age","AverageDBH" ,"DBH_CV"  ,"Density") ~ "Stand structural",
    var %in% c( "PD" ,"pd_Broadleaf"  ,"FDis" ,"FDisB") ~ "Diversity",
    var %in% c("Altitude",
               "MAT", "MAP" , "SoilThickness" ,"HumusLayerThickness", "LitterThickness") ~ "Geographic Location",
    TRUE ~ NA_character_  # In case there are any terms not listed above
  ))

#cols1=c("#DDE9F7", "#B5D3E9", "#82BADA","#539DCD","#004b79","#436f8e"  "#7E9E8D")
cols1=c(
  "#9AB560" ,
  "#463D09" ,
  #"#D7D1C2",
  "#fff9ea"
)
df1=df
cols1=c( #"#DA561D" , 
  "#9AB560" ,
  "#463D09",
  #"#D7D1C2",
  "#fff9ea"
)

p_all2_L=df2 %>% 
  dplyr::mutate(group=fct_relevel(group,
                           rev(c("Stand structural", "Diversity",#"CWM" ,
                                 "Geographic Location"
                           )))) %>% 
  group_by(group) %>% 
  dplyr::summarise(sum_value=sum(abs(Estimate))) %>% 
  dplyr::mutate(new_col=sum_value/sum(sum_value)*100) %>% 
  ggplot(aes(x = 1, y = new_col, label = group)) +
  geom_col(aes(fill = group, color = group), alpha =1, show.legend = FALSE) +
  scale_fill_manual(values = cols1) +
  #scale_color_manual(values = c("#B7D982" ,"#B7D982" , "#B7D982" ,"#B7D982" )    ) +
  scale_color_manual(values =cols1  ) + 
  coord_flip()+
  
  scale_fill_manual(values = cols1)+
  scale_y_continuous(expand = c(0,0), limits = c(0,100), breaks = c(0,20,40,60,80,100))+
  
  
  theme_minimal()+
  theme(
    axis.text.x = element_blank(),  
    #axis.text.x = element_text(size = 10, family =  "serif"),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),  
    axis.line.y = element_blank(),
    axis.ticks.x = element_blank( ), 
    axis.text.y = element_blank(),  
    plot.title = element_blank()
  )+
  
  
  labs(y="Relative effect of estimates (%)")+
  
  labs(subtitle  = expression(Adj.~ italic(R)^2==84~ '%' ))+
  #labs(x="Pinus massoniana")+
  labs(title = expression(Low-carbon-sequestration~forests))+
  
  #ggtitle(" ") +
  
  theme(plot.title = element_text(hjust = 0, face = "bold"))+
  
  theme(plot.subtitle = element_text(size = 14, family = "serif") 
  )+
  theme(plot.title = element_text(size = 14, 
                                  family =  "serif"))
p_all2_L
#################
df3<-read.csv("FD_full_average_High_310.csv",check.names = F)

colnames(df3) <- c("var" ,"Estimate", "Std. Error","Adjusted.SE", "z.value","Pr(>|z|)"  )

unique(df3$var)

df3 <- df3 %>%
  dplyr::mutate(across(where(is.numeric), round, 2))


df3 <- df3 %>%
  filter(!( var %in% c("(Intercept)") )  ) %>%
  dplyr::mutate(group = case_when(
    var %in% c("Age","AverageDBH" ,"DBH_CV"  ,"Density") ~ "Stand structural",
    var %in% c( "PD" ,"pd_Broadleaf"  ,"FDis" ,"FDisB") ~ "Diversity",
    var %in% c("Altitude",
               "MAT", "MAP" , "SoilThickness" ,"HumusLayerThickness", "LitterThickness") ~ "Geographic Location",
    TRUE ~ NA_character_  # In case there are any terms not listed above
  ))
df3 %>% 
  dplyr::mutate(group=fct_relevel(var,
                           rev(c()))) %>% 
  group_by(group) %>% 
  dplyr::summarise(sum_value=sum(abs(Estimate))) %>% 
  dplyr::mutate(new_col=sum_value/sum(sum_value))
df3
df3 %>% 
  # dplyr::mutate(group=fct_relevel(var,
  # rev(c()))) %>% #求每类变量
  group_by(group) %>% 
  dplyr::summarise(sum_value=sum(abs(Estimate))) %>% 
  dplyr::mutate(new_col=sum_value/sum(sum_value))

df3 %>% 
  # dplyr::mutate(group=fct_relevel(var,
  # rev(c()))) %>% #求每类变量
  group_by(group) %>% 
  dplyr::summarise(sum_value=sum(abs(Estimate))) %>% 
  dplyr::mutate(new_col=sum_value/sum(sum_value))

#cols1=c("#DDE9F7", "#B5D3E9", "#82BADA","#539DCD","#004b79","#436f8e"  "#7E9E8D")
cols1=c( #"#DA561D" , 
  "#9AB560" ,
  "#463D09" ,
  #"#D7D1C2",
  "#fff9ea"
)


p_all2_H=df3 %>% 
  dplyr::mutate(group=fct_relevel(group,
                           rev(c("Stand structural", "Diversity",#"CWM",
                                 "Geographic Location")))) %>% 
  group_by(group) %>% 
  dplyr::summarise(sum_value=sum(abs(Estimate))) %>% 
  dplyr::mutate(new_col=sum_value/sum(sum_value)*100) %>% 
  ggplot(aes(x = 1, y = new_col, label = group)) +
  geom_col(aes(fill = group, color = group), alpha =1, show.legend = FALSE) +
  scale_fill_manual(values = cols1) +
  #scale_color_manual(values = c("#1572CD"  ,"#1572CD" ,"#1572CD" ,"#1572CD"      )    ) + 
  scale_color_manual(values = cols1    ) + 
  coord_flip()+
  
  scale_fill_manual(values = cols1)+
  scale_y_continuous(expand = c(0,0), limits = c(0,100), breaks = c(0,20,40,60,80,100))+
  
  
  theme_minimal()+
  theme(
    axis.text.x = element_blank(),  
    #axis.text.x = element_text(size = 10, family =  "serif"),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),  
    axis.line.y = element_blank(),
    axis.ticks.x = element_blank( ), 
    axis.text.y = element_blank(),  
    plot.title = element_blank()
  )+
 
  labs(y="Relative effect of estimates (%)")+
  
  labs(subtitle  = expression(Adj.~ italic(R)^2==72~ '%' ))+
  #labs(x="Pinus massoniana")+
  labs(title = expression(High-carbon-sequestration~forests))+
  
  #ggtitle(" ") +
  
  theme(plot.title = element_text(hjust = 0, face = "bold"))+
  
  theme(plot.subtitle = element_text(size = 14, family = "serif") 
  )+
  theme(plot.title = element_text(size = 14, 
                                  family =  "serif"))

p_all2_H

####ggarrange######################
library(patchwork)
library(patchwork)
pG1=ggarrange(p_all2_H,p_all2_L,p_all2_G ,
              ncol =1,nrow = 3,
              heights = c(1,1),
              widths = c(1,1),common.legend = FALSE, align = c("hv"))
pG1
##################################################################################
####Pinus massoniana：Plot2：100%Plot############
dfP<-read.csv("FD_full_average_All_220.csv",check.names = F)

colnames(dfP) <- c("var" ,"Estimate", "Std. Error","Adjusted.SE", "z.value","Pr(>|z|)"  )

unique(dfP$var)

dfP <- dfP %>%
  dplyr::mutate(across(where(is.numeric), round, 2))


dfP <- dfP %>%
  
  filter(!( var %in% c("(Intercept)") )  ) %>%
  dplyr::mutate(group = case_when(
    var %in% c("Age","AverageDBH"  ,"DBH_CV"  ,"Density") ~ "Stand structural",
    var %in% c(  "PD" ,"pd_Broadleaf"  ,"FDis" ,"FDisB") ~ "Diversity",
    var %in% c("Altitude",
               "MAT" , "MAP"  , "SoilThickness" ,"HumusLayerThickness", "LitterThickness") ~ "Geographic Location",
    TRUE ~ NA_character_  # In case there are any terms not listed above
  ))

cols1=c( 
  "#9AB560" ,
  "#463D09" ,
  #"#D7D1C2",
  "#fff9ea"
)

p_all2_PG=dfP %>% 
  dplyr::mutate(group=fct_relevel(group,
                           rev(c("Stand structural", "Diversity",#"CWM",
                                 "Geographic Location")))) %>% 
  group_by(group) %>% 
  dplyr::summarise(sum_value=sum(abs(Estimate))) %>% 
  dplyr::mutate(new_col=sum_value/sum(sum_value)*100) %>% 
  ggplot(aes(x = 1, y = new_col, label = group)) +
  geom_col(aes(fill = group, color = group), alpha =1, show.legend = FALSE) +
  scale_fill_manual(values = cols1) +
  scale_color_manual(values = cols1    ) +
  #scale_color_manual(values = c("#a1a1a4","#a1a1a4","#a1a1a4" ,"#a1a1a4" )    ) +
  coord_flip()+
  
  scale_fill_manual(values = cols1)+
  scale_y_continuous(expand = c(0,0), limits = c(0,100), breaks = c(0,20,40,60,80,100))+
  
  
  theme_minimal()+
  theme(
    axis.text.x = element_text(size = 12, family = "serif"),
    axis.title.y = element_blank(),  
    axis.line.y = element_blank(),   
    axis.line.x = element_line(),
    axis.ticks.x = element_line(),
    axis.text.y = element_blank(),
    plot.title = element_blank(),
    axis.title.x = element_text(size = 14, 
                                family =  "serif")
  ) +
  labs(y="Relative effect of estimates (%)")+
  
  labs(subtitle  = expression(Adj.~ italic(R)^2==78~ '%' ))+
  #labs(x="Pinus massoniana")+
  labs(title = expression(General))+
  
  #ggtitle(" ") +
  
  theme(plot.title = element_text(hjust = 0, face = "bold")
        
        
  )+
  
  theme(plot.subtitle = element_text(size = 14, family = "serif") 
  )+
  theme(plot.title = element_text(size = 14, 
                                  family =  "serif"))
p_all2_PG
#################
dfP2<-read.csv("FD_full_average_Low_220.csv",check.names = F)

colnames(dfP2) <- c("var" ,"Estimate", "Std. Error","Adjusted.SE", "z.value","Pr(>|z|)"  )

unique(dfP2$var)

dfP2 <- dfP2 %>%
  dplyr::mutate(across(where(is.numeric), round, 2))


dfP2 <- dfP2 %>%
  
  filter(!( var %in% c("(Intercept)") )  ) %>%
  dplyr::mutate(group = case_when(
    var %in% c("Age","AverageDBH"  ,"DBH_CV"  ,"Density") ~ "Stand structural",
    var %in% c(  "PD" ,"pd_Broadleaf"  ,"FDis" ,"FDisB") ~ "Diversity",
    var %in% c("Altitude",
               "MAT" , "MAP"  , "SoilThickness" ,"HumusLayerThickness", "LitterThickness") ~ "Geographic Location",
    TRUE ~ NA_character_  # In case there are any terms not listed above
  ))

cols1=c( 
  "#9AB560" ,
 # "#463D09" ,
  #"#D7D1C2",
  "#fff9ea"
)
p_all2_PL=dfP2 %>% 
  dplyr::mutate(group=fct_relevel(group,
                           rev(c("Stand structural", #"Diversity",
                                 "Geographic Location")))) %>% 
  group_by(group) %>% 
  dplyr::summarise(sum_value=sum(abs(Estimate))) %>% 
  dplyr::mutate(new_col=sum_value/sum(sum_value)*100) %>% 
  ggplot(aes(x = 1, y = new_col, label = group)) +
  geom_col(aes(fill = group, color = group), alpha =1, show.legend = FALSE) +
  scale_fill_manual(values = cols1) +
  #scale_color_manual(values = cols1    ) + 
  scale_color_manual(values = cols1 ) +
  coord_flip()+
  
  scale_fill_manual(values = cols1)+
  scale_y_continuous(expand = c(0,0), limits = c(0,100), breaks = c(0,20,40,60,80,100))+
  
  
  theme_minimal()+
  theme(
    axis.text.x = element_blank(),  
    #axis.text.x = element_text(size = 10, family =  "serif"),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),  
    axis.line.y = element_blank(),
    axis.ticks.x = element_blank( ), 
    axis.text.y = element_blank(),  
    plot.title = element_blank()
  )+
  
  
  labs(y="Relative effect of estimates (%)")+
  
  labs(subtitle  = expression(Adj.~ italic(R)^2==91~ '%' ))+
  #labs(x="Pinus massoniana")+
  labs(title = expression(Low-carbon-sequestration~forests))+
  
  #ggtitle(" ") +
  
  theme(plot.title = element_text(hjust = 0, face = "bold"))+
  
  theme(plot.subtitle = element_text(size = 14, family = "serif") 
  )+
  theme(plot.title = element_text(size = 14, 
                                  family =  "serif"))
p_all2_PL
#################
dfP3<-read.csv("FD_full_average_High_220.csv",check.names = F)

colnames(dfP3) <- c("var" ,"Estimate", "Std. Error","Adjusted.SE", "z.value","Pr(>|z|)"  )

unique(dfP3$var)

dfP3 <- dfP3 %>%
  dplyr::mutate(across(where(is.numeric), round, 2))


dfP3 <- dfP3 %>%
  
  filter(!( var %in% c("(Intercept)") )  ) %>%
  dplyr::mutate(group = case_when(
    var %in% c("Age","AverageDBH"  ,"DBH_CV"  ,"Density") ~ "Stand structural",
    var %in% c(  "PD" ,"pd_Broadleaf"  ,"FDis" ,"FDisB") ~ "Diversity",
    var %in% c("Altitude",
               "MAT" , "MAP"  , "SoilThickness" ,"HumusLayerThickness", "LitterThickness") ~ "Geographic Location",
    TRUE ~ NA_character_  # In case there are any terms not listed above
  ))




cols1=c( 
  "#9AB560" ,
  "#463D09" ,
  #"#D7D1C2",
  "#fff9ea"
)

p_all2_PH=dfP3 %>% 
  dplyr::mutate(group=fct_relevel(group,
                           rev(c("Stand structural","Diversity", "Geographic Location")))) %>% 
  group_by(group) %>% 
  dplyr::summarise(sum_value=sum(abs(Estimate))) %>% 
  dplyr::mutate(new_col=sum_value/sum(sum_value)*100) %>% 
  ggplot(aes(x = 1, y = new_col, label = group)) +
  geom_col(aes(fill = group, color = group), alpha =1, show.legend = FALSE) +
  scale_fill_manual(values = cols1) +
  scale_color_manual(values = cols1) +
  #scale_color_manual(values = c("#1572CD"  ,"#1572CD" ,"#1572CD" ,"#1572CD"      )    ) + 
  coord_flip()+
  
  scale_fill_manual(values = cols1)+
  scale_y_continuous(expand = c(0,0), limits = c(0,100), breaks = c(0,20,40,60,80,100))+
  
  
  theme_minimal()+
  theme(
    axis.text.x = element_blank(),  
    #axis.text.x = element_text(size = 10, family =  "serif"),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),  
    axis.line.y = element_blank(),
    axis.ticks.x = element_blank( ), 
    axis.text.y = element_blank(),  
    plot.title = element_blank()
  )+
 
  labs(y="Relative effect of estimates (%)")+
  
  labs(subtitle  = expression(Adj.~ italic(R)^2==75~ '%' ))+
  #labs(x="Pinus massoniana")+
  labs(title = expression(High-carbon-sequestration~forests))+
  
  #ggtitle(" ") +
  
  theme(plot.title = element_text(hjust = 0, face = "bold"))+
  
  theme(plot.subtitle = element_text(size = 14, family = "serif") 
  )+
  theme(plot.title = element_text(size = 14, 
                                  family =  "serif"))

p_all2_PH

##########################
library(patchwork)
library(patchwork)
pGP1=ggarrange(p_all2_PH,p_all2_PL,p_all2_PG ,
               ncol =1,nrow = 3,
               heights = c(1,1),
               widths = c(1,1),common.legend = FALSE, align = c("hv"))
pGP1

######ggarrange##########
pG=ggarrange(#p_all1,
  pG1,
  #p_all1_P,
  pGP1,
  ncol =2,nrow =1,
  heights = c(1,1),
  #widths = c(1.6,0.4),
  widths = c(1,1),
  common.legend = FALSE, align = c("hv"))
pG

ggsave("Plot_Liear2_FD_250329.png", path = "E:/MYData/Compare_P_C/figure/Plot_Liear",width =6, height =4.2,dpi=1200, plot=pG)

##===============================================================================
## PD:LMMPlot ---------
##===============================================================================
#######Cunninghamia lanceolata：plot_all######
rm(list=ls())
library(tidyverse)
library(ggplot2)
library(ggh4x)
library(dplyr)

df_All<-read.csv("PD_full_average_All_310.csv",check.names = F)
df_High<-read.csv("PD_full_average_High_310.csv",check.names = F)
df_Low<-read.csv("PD_full_average_Low_310.csv",check.names = F)

df_All$Grid<-"General"
df_High$Grid<-"High"
df_Low$Grid<-"Low"

df=rbind( df_All,df_High,df_Low )
colnames(df) <- c("var" ,"Estimate", "Std. Error","Adjusted.SE", "z.value","Pr(>|z|)", "Grid"   )

df <- df %>%
  filter(!( var %in% c("(Intercept)") )  ) %>%
  dplyr::mutate(group = case_when(
    var %in% c("Age","AverageDBH" ,"DBH_CV"  ,"Density") ~ "Stand structural",
    var %in% c( "PD" ,"pd_Broadleaf"  ,"Shannon" ,"shan_Broadleaf") ~ "Diversity",
    var %in% c("Altitude",
               "MAT", "MAP" , "SoilThickness" ,"HumusLayerThickness", "LitterThickness") ~ "Geographic Location",
    TRUE ~ NA_character_  # In case there are any terms not listed above
  ))



cols<-c("#ADD9E6" , "#d2de61" ,"#dd7500")


####To calculate the percentage for each category########

B=df %>% 
  group_by(Grid,var
           ) %>%
  dplyr::mutate(groups=fct_relevel(var,
                            rev(c()))) %>% 
  group_by(Grid,
           #group,
           var
           ) %>% 
  dplyr::summarise(sum_value = sum(abs(Estimate), na.rm = TRUE)) %>% 
  dplyr::mutate(new_col = sum_value / sum(sum_value, na.rm = TRUE))  %>% 
  dplyr::mutate(new_col = round(new_col, 4))  
B
####Plot1########

df <- df %>%
  dplyr::mutate(
    sign = if_else(
      Grid == "Low", Estimate +2*abs(`Std. Error`) ,
      if_else(Grid == "High", Estimate + 2*abs( `Std. Error`),
              if_else(Grid == "General", Estimate + 2*abs( `Std. Error`), NA_real_))
    )
  )
custom_order <- c("AverageDBH" ,"DBH_CV"  ,"Density",
                 # "pd_Broadleaf"  ,
                  "PD" ,"shan_Broadleaf","Shannon" ,"Age",
                  "Altitude",
                  "MAT", "MAP" , "SoilThickness" ,"HumusLayerThickness", "LitterThickness"  )

unique(df$var)
p_all310= df %>%
  dplyr::mutate(Grid = fct_relevel(Grid, c("General", "Low","High")))  %>% 
  dplyr::mutate(group1=fct_relevel(var,custom_order))  %>% 
  group_by(group1, Grid) %>% 
  arrange(Estimate) %>% 
  dplyr::mutate(var = fct_relevel(var, var)) %>% 
  dplyr::mutate(signi = case_when(
    `Pr(>|z|)` < 0.1 & `Pr(>|z|)` >= 0.05 ~ '',
    `Pr(>|z|)` > 0.1 ~ '',
    `Pr(>|z|)` > 0.05 ~ '',
    `Pr(>|z|)` < 0.05 & `Pr(>|z|)` >= 0.01 ~ '*',
    `Pr(>|z|)` < 0.01 & `Pr(>|z|)` >= 0.001 ~ '**',
    `Pr(>|z|)` < 0.001 ~ '***'
  )) %>% 
  
  arrange(group) %>%  
  ggplot(aes(x = Estimate, y = fct_reorder(var, desc(group1)), color = Grid, group = Grid)) +
  geom_vline(xintercept=c(0), linetype="longdash")+
  geom_linerange(
    aes(xmin = Estimate - `Std. Error`- `Std. Error`, xmax = Estimate + `Std. Error`+ `Std. Error`, 
        color = Grid, group = Grid),
    size =0.1, 
    color = "#040505",
    alpha =1, position = position_dodge(width = 0.6)
  ) +
  geom_point(
    aes(color = Grid, group = Grid),
    show.legend = T,
    size = 3, alpha = 1, 
    #shape=21,
    #color = "black",
    #fill = "black",
    position = position_dodge(width = 0.6)
  )+
  scale_shape_manual(values = c(21,21,21,21,21)) +
  
  #scale_color_manual(values = c("General"="#a1a1a4", "Low"="#B7D982","High"="#1572CD" )) +
  #scale_fill_manual(values =  c("General"="#a1a1a4", "Low"="#B7D982","High"="#1572CD" )) +
  

  geom_text(aes(x = sign , y = var,label = signi),
            # size =3, family = "serif",  color ="black",fontface = "bold", vjust = -0.3,position = position_dodge(width = 0.6), hjust = 0.5) +
            size =3, family = "serif",  color ="black",
            fontface = "bold", position = position_dodge(width = 0.6), hjust = 0) +
  
  geom_hline(yintercept = seq(0.5, length(custom_order) - 0.5), linetype = "solid", color = "#f3f4f7") +  # 添加横线
  
  theme_bw( )+
  
  labs(title = expression((a)~italic(Cunninghamia~lanceolata)))+
  theme(plot.title = element_text(size=16,family="serif"))+
  
  labs(x="Parameter Estimate",y=NULL)+
  theme(axis.title.x = element_text(size = 14, family =  "serif"))+
  #labs(x=NULL,y=NULL)+
  
  #theme(axis.title.x = element_text(size = 12, family =  "serif"))+
  scale_color_manual(values = cols)+
  scale_x_continuous(limits = c(-0.8,1),
                     breaks = c(-0.5,0,0.5,1))+
  
  
  scale_y_discrete(
    labels = c(
      
      expression(Average~DBH),
      expression(DBH~variation),
      expression(Stand~density),
      
      #expression(PD[Broadleaf]),
      expression(PD),
      #expression(Shannon[Broadleaf]),
      #expression(Shannon),
      expression(H*"'"[Broadleaves]), expression(H*"'"),
      
      expression(Stand~age),
      
      
      expression(Elevation),
      expression(MAT),
      expression(MAP),
      expression(ST),
      expression(HLT),
      expression(LT)
    ),
    limits = rev#,position = "right"
  )+
  theme(#axis.text.x = element_blank(),
    axis.ticks = element_line(),
    axis.text.x = element_text(size = 12,color="black",family="serif"),
    axis.text.y = element_text(size = 14,color="black",face="bold",family="serif"),
    #axis.text.y =  element_blank(),
    axis.ticks.y = element_line()
  )+
  theme(
    panel.grid = element_blank(),
    # panel.background = element_rect(fill = "#ffffff"),  
    # panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # 设置面板边框线
    
    # panel.grid.major = element_line(color = "#f3f4f7" ),  
    #panel.grid.minor = element_line(color = "#f3f4f7"),  # 设置次要网格线
    
    axis.line = element_line(color = "black", size = 0.5),  # 设置坐标轴线
    
    axis.title = element_text(color = "black"),  # 设置坐标轴标题颜色
    axis.text = element_text(color = "black"),  # 设置坐标轴文本颜色
    axis.line.x = element_line(color = "black", size = 0.5),  # 设置x轴线颜色和大小
    axis.line.y = element_line(color = "black", size = 0.5) , # 设置y轴线颜色和大小
    #axis.line.y = element_blank() , # 设置y轴线颜色和大小
    legend.position = ""
  )
  
 


p_all310

#ggsave("Plot_LiearPlot例.png", path = "E:/MYData/Compare_P_C/figure/Plot_Liear",width =6, height =4.2,dpi=1200, plot=p_all1)
################
#######Pinus massoniana：plot_all######

library(tidyverse)
library(ggplot2)
library(ggh4x)
library(dplyr)

df_All<-read.csv("PD_full_average_All_220.csv",check.names = F)
df_High<-read.csv("PD_full_average_High_220.csv",check.names = F)
df_Low<-read.csv("PD_full_average_Low_220.csv",check.names = F)

df_All$Grid<-"General"
df_High$Grid<-"High"
df_Low$Grid<-"Low"

df=rbind( df_All,df_High,df_Low )
colnames(df) <- c("var" ,"Estimate", "Std. Error","Adjusted.SE", "z.value","Pr(>|z|)", "Grid"   )
unique(df$var)
df <- df %>%
  filter(!( var %in% c("(Intercept)") )  ) %>%
  dplyr::mutate(group = case_when(
    var %in% c("Age","AverageDBH"  ,"DBH_CV"  ,"Density") ~ "Stand structural",
    var %in% c( "PD" ,"pd_Broadleaf"  ,"Shannon" ,"shan_Broadleaf") ~ "Diversity",
    var %in% c("Altitude",
               "MAT" , "MAP"  , "SoilThickness" ,"HumusLayerThickness", "LitterThickness") ~ "Geographic Location",
    TRUE ~ NA_character_  # In case there are any terms not listed above
  ))

unique(df$var)

cols<-c("#ADD9E6" , "#d2de61" ,"#dd7500")


####To calculate the percentage for each category########

B=df %>% 
  group_by(Grid,var
  ) %>%
  dplyr::mutate(groups=fct_relevel(var,
                                   rev(c()))) %>% 
  group_by(Grid,
           #group,
           var
  ) %>% 
  dplyr::summarise(sum_value = sum(abs(Estimate), na.rm = TRUE)) %>% 
  dplyr::mutate(new_col = sum_value / sum(sum_value, na.rm = TRUE))  %>% 
  dplyr::mutate(new_col = round(new_col, 4))  
B
####Plot1########

df <- df %>%
  dplyr::mutate(
    sign = if_else(
      Grid == "Low", Estimate +2*abs(`Std. Error`) ,
      if_else(Grid == "High", Estimate + 2*abs( `Std. Error`),
              if_else(Grid == "General", Estimate + 2*abs( `Std. Error`), NA_real_))
    )
  )



custom_order <- c("AverageDBH"  ,"DBH_CV"  ,"Density",
                 # "pd_Broadleaf"  ,
                  "PD" ,"shan_Broadleaf","Shannon" ,"Age",
                  "Altitude",
                  "MAT" , 
                  "MAP"  , "SoilThickness" ,"HumusLayerThickness", "LitterThickness" 
                  )

unique(df$var)
p_all220= df %>%
  dplyr::mutate(Grid = fct_relevel(Grid, c("General", "Low","High")))  %>%  
  dplyr::mutate(group1=fct_relevel(var,custom_order))  %>% 
  group_by(group1, Grid) %>% 
  arrange(Estimate) %>% 
  dplyr::mutate(var = fct_relevel(var, var)) %>% 
  dplyr::mutate(signi = case_when(
    `Pr(>|z|)` < 0.1 & `Pr(>|z|)` >= 0.05 ~ '',
    `Pr(>|z|)` > 0.1 ~ '',
    `Pr(>|z|)` > 0.05 ~ '',
    `Pr(>|z|)` < 0.05 & `Pr(>|z|)` >= 0.01 ~ '*',
    `Pr(>|z|)` < 0.01 & `Pr(>|z|)` >= 0.001 ~ '**',
    `Pr(>|z|)` < 0.001 ~ '***'
  )) %>% 
  
  arrange(group) %>%  
  ggplot(aes(x = Estimate, y = fct_reorder(var, desc(group1)), color = Grid, group = Grid)) +
  geom_vline(xintercept=c(0), linetype="longdash")+
  geom_linerange(
    aes(xmin = Estimate - `Std. Error`- `Std. Error`, xmax = Estimate + `Std. Error`+ `Std. Error`, color = Grid, group = Grid),
    size =0.1, 
    color = "#040505",
    alpha =1, position = position_dodge(width = 0.6)
  ) +
  geom_point(
    aes(color = Grid, group = Grid),
    show.legend = F,
    size = 3, alpha = 1, 
    #shape=21,
    #color = "black",
    #fill = "black",
    position = position_dodge(width = 0.6)
  )+
  scale_shape_manual(values = c(21,21,21,21,21)) +
  
  #scale_color_manual(values = c("General"="#a1a1a4", "Low"="#B7D982","High"="#1572CD" )) +
  #scale_fill_manual(values =  c("General"="#a1a1a4", "Low"="#B7D982","High"="#1572CD" )) +
  

  geom_text(aes(x = sign , y = var,label = signi),
            # size =3, family = "serif",  color ="black",fontface = "bold", vjust = -0.3,position = position_dodge(width = 0.6), hjust = 0.5) +
            size =3, family = "serif",  color ="black",
            fontface = "bold", position = position_dodge(width = 0.6), hjust = 0) +
  
  geom_hline(yintercept = seq(0.5, length(custom_order) - 0.5), linetype = "solid", color = "#f3f4f7") +  # 添加横线
  
  theme_bw( )+
  
  labs(title = expression((b)~italic(Pinus~massoniana)))+
  theme(plot.title = element_text(size=16,family="serif"))+
  
  labs(x="Parameter Estimate",y=NULL)+
  theme(axis.title.x = element_text(size = 14, family =  "serif"))+
  #labs(x=NULL,y=NULL)+
  
  #theme(axis.title.x = element_text(size = 12, family =  "serif"))+
  
  scale_x_continuous(limits = c(-0.8,1.1),
                      breaks = c(-0.5,0,0.5,1))+
  
  
  scale_y_discrete(
    labels = c(
     
      expression(Average~DBH),
      expression(DBH~variation),
      expression(Stand~density),
      
      #expression(PD[Broadleaf]),
      expression(PD),
     # expression(Shannon[Broadleaf]),
      #expression(Shannon), 
     expression(H*"'"[Broadleaf]), expression(H*"'"),
      expression(Stand~age),
      
      expression(Elevation),
      expression(MAT),
      expression(MAP),
      expression(ST),
      expression(HLT),
      expression(LT)
    ),
    limits = rev#,position = "right"
  )+
  theme(#axis.text.x = element_blank(),
    axis.ticks = element_line(),
    axis.text.x = element_text(size = 12,color="black",family="serif"),
    axis.text.y = element_text(size = 14,color="black",face="bold",family="serif"),
    #axis.text.y =  element_blank(),
    axis.ticks.y = element_line()
  )+
  theme(
    panel.grid = element_blank(),
    # panel.background = element_rect(fill = "#ffffff"),  
    # panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # 设置面板边框线
    
    # panel.grid.major = element_line(color = "#f3f4f7" ),  
    #panel.grid.minor = element_line(color = "#f3f4f7"),  # 设置次要网格线
    
    axis.line = element_line(color = "black", size = 0.5),  # 设置坐标轴线
    
    axis.title = element_text(color = "black"),  # 设置坐标轴标题颜色
    axis.text = element_text(color = "black"),  # 设置坐标轴文本颜色
    axis.line.x = element_line(color = "black", size = 0.5),  # 设置x轴线颜色和大小
    axis.line.y = element_line(color = "black", size = 0.5) , # 设置y轴线颜色和大小
    #axis.line.y = element_blank() , # 设置y轴线颜色和大小
    legend.position = ""
  )+
  
  guides(fill = guide_legend(reverse = TRUE))+
  theme(legend.position = "right", legend.title = element_text(size = 20),
         legend.text = element_text(size = 12,family="serif"#,face = "italic"
                                    ))+
  
  scale_color_manual(values = cols)


p_all220

#ggsave("Plot_LiearPlot例.png", path = "E:/MYData/Compare_P_C/figure/Plot_Liear",width =6, height =4.2,dpi=1200, plot=p_all1)

####ggarrange####
library(patchwork)
library(ggpubr)
Plot_Liear1=ggarrange(p_all310,p_all220,
                      ncol = 1,nrow =2,
                      heights = c(1,1),
                      widths = c(1,1),common.legend = FALSE, align = c("hv"))

Plot_Liear1
ggsave("Plot_Liear_PD250415.png", path = "E:/MYData/Compare_P_C/figure/Plot_Liear",width =5, height =8,dpi=600, plot=Plot_Liear1)


#################################################################################
#################################################################################
####Cunninghamia lanceolata：Plot2：100%Plot############
df<-read.csv("PD_full_average_All_310.csv",check.names = F)

colnames(df) <- c("var" ,"Estimate", "Std. Error","Adjusted.SE", "z.value","Pr(>|z|)"  )

unique(df$var)

# df <- df %>%
#   dplyr::mutate(across(where(is.numeric), round,2))



df <- df %>%
  filter(!( var %in% c("(Intercept)") )  ) %>%
  dplyr::mutate(group = case_when(
    var %in% c("Age","AverageDBH" ,"DBH_CV"  ,"Density") ~ "Stand structural",
    var %in% c( "PD" ,"pd_Broadleaf"  ,"Shannon" ,"shan_Broadleaf") ~ "Diversity",
    var %in% c("Altitude",
               "MAT", "MAP" , "SoilThickness" ,"HumusLayerThickness", "LitterThickness") ~ "Geographic Location",
    TRUE ~ NA_character_  # In case there are any terms not listed above
  ))



#cols1=c("#DDE9F7", "#B5D3E9", "#82BADA","#539DCD","#004b79","#436f8e"  "#7E9E8D")
cols1=c(
  "#9AB560" ,
  "#463D09" ,
  #"#D7D1C2",
  "#fff9ea"
)
df1=df

p_all2_G=df1 %>% 
  dplyr::mutate(group=fct_relevel(group,
                           rev(c("Stand structural", "Diversity",
                                 #"CWM",
                                 "Geographic Location"
                           )))) %>% 
  group_by(group) %>% 
  dplyr::summarise(sum_value=sum(abs(Estimate))) %>% 
  dplyr::mutate(
    new_col=sum_value/sum(sum_value)*100) %>% 
  ggplot(aes(x = 1, y = new_col, label = group)) +
  geom_col(aes(fill = group, color = group), alpha =1, show.legend = FALSE) +
  scale_fill_manual(values = cols1) +
  scale_color_manual(values =cols1    ) +
  #scale_color_manual(values = c("#C0BDAB", "#C0BDAB","#C0BDAB","#C0BDAB"  )    ) + 
  coord_flip()+
  
  scale_fill_manual(values = cols1)+
  scale_y_continuous(expand = c(0,0), limits = c(0,100), breaks = c(0,20,40,60,80,100))+
  
  
  theme_minimal()+
  theme(
    axis.text.x = element_text(size = 12, family = "serif"),
    axis.title.y = element_blank(),  
    axis.line.y = element_blank(),   
    axis.line.x = element_line(),
    axis.ticks.x = element_line(),
    axis.text.y = element_blank(),
    plot.title = element_blank(),
    axis.title.x = element_text(size = 14, 
                                family =  "serif")
  ) +
  
 
  
  labs(y="Relative effect of estimates (%)")+
  
  labs(subtitle  = expression(Adj.~ italic(R)^2==61~ '%' ))+
  #labs(x="Pinus massoniana")+
  labs(title = expression(General))+
  
  #ggtitle(" ") +
  
  theme(plot.title = element_text(hjust = 0, face = "bold"))+
  
  theme(plot.subtitle = element_text(size = 14, family = "serif") 
  )+
  theme(plot.title = element_text(size = 14, 
                                  family =  "serif"))
p_all2_G
#################
df2<-read.csv("PD_full_average_Low_310.csv",check.names = F)

colnames(df2) <- c("var" ,"Estimate", "Std. Error","Adjusted.SE", "z.value","Pr(>|z|)"  )

unique(df2$var)

df2 <- df2 %>%
  dplyr::mutate(across(where(is.numeric), round, 2))



df2 <- df2 %>%
  filter(!( var %in% c("(Intercept)") )  ) %>%
  dplyr::mutate(group = case_when(
    var %in% c("Age","AverageDBH" ,"DBH_CV"  ,"Density") ~ "Stand structural",
    var %in% c( "PD" ,"pd_Broadleaf"  ,"Shannon" ,"shan_Broadleaf") ~ "Diversity",
    var %in% c("Altitude",
               "MAT", "MAP" , "SoilThickness" ,"HumusLayerThickness", "LitterThickness") ~ "Geographic Location",
    TRUE ~ NA_character_  # In case there are any terms not listed above
  ))



#cols1=c("#DDE9F7", "#B5D3E9", "#82BADA","#539DCD","#004b79","#436f8e"  "#7E9E8D")
cols1=c(
  "#9AB560" ,
  "#463D09" ,
  #"#D7D1C2",
  "#fff9ea"
)
df1=df
cols1=c( #"#DA561D" , 
  "#9AB560" ,
  "#463D09",
  #"#D7D1C2",
  "#fff9ea"
)

p_all2_L=df2 %>% 
  dplyr::mutate(group=fct_relevel(group,
                           rev(c("Stand structural", "Diversity",#"CWM" ,
                                 "Geographic Location"
                           )))) %>% 
  group_by(group) %>% 
  dplyr::summarise(sum_value=sum(abs(Estimate))) %>% 
  dplyr::mutate(new_col=sum_value/sum(sum_value)*100) %>% 
  ggplot(aes(x = 1, y = new_col, label = group)) +
  geom_col(aes(fill = group, color = group), alpha =1, show.legend = FALSE) +
  scale_fill_manual(values = cols1) +
  #scale_color_manual(values = c("#B7D982" ,"#B7D982" , "#B7D982" ,"#B7D982" )    ) +
  scale_color_manual(values =cols1  ) + 
  coord_flip()+
  
  scale_fill_manual(values = cols1)+
  scale_y_continuous(expand = c(0,0), limits = c(0,100), breaks = c(0,20,40,60,80,100))+
  
  
  theme_minimal()+
  theme(
    axis.text.x = element_blank(),  
    #axis.text.x = element_text(size = 10, family =  "serif"),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),  
    axis.line.y = element_blank(),
    axis.ticks.x = element_blank( ), 
    axis.text.y = element_blank(),  
    plot.title = element_blank()
  )+
  
  
  labs(y="Relative effect of estimates (%)")+
  
  labs(subtitle  = expression(Adj.~ italic(R)^2==83~ '%' ))+
  #labs(x="Pinus massoniana")+
  labs(title = expression(Low-carbon-sequestration~forests))+
  
  #ggtitle(" ") +
  
  theme(plot.title = element_text(hjust = 0, face = "bold"))+
  
  theme(plot.subtitle = element_text(size = 14, family = "serif") 
  )+
  theme(plot.title = element_text(size = 14, 
                                  family =  "serif"))
p_all2_L
#################
df3<-read.csv("PD_full_average_High_310.csv",check.names = F)

colnames(df3) <- c("var" ,"Estimate", "Std. Error","Adjusted.SE", "z.value","Pr(>|z|)"  )

unique(df3$var)

df3 <- df3 %>%
  dplyr::mutate(across(where(is.numeric), round, 2))



df3 <- df3 %>%
  filter(!( var %in% c("(Intercept)") )  ) %>%
  dplyr::mutate(group = case_when(
    var %in% c("Age","AverageDBH" ,"DBH_CV"  ,"Density") ~ "Stand structural",
    var %in% c( "PD" ,"pd_Broadleaf"  ,"Shannon" ,"shan_Broadleaf") ~ "Diversity",
    var %in% c("Altitude",
               "MAT", "MAP" , "SoilThickness" ,"HumusLayerThickness", "LitterThickness") ~ "Geographic Location",
    TRUE ~ NA_character_  # In case there are any terms not listed above
  ))


df3 %>% 
  dplyr::mutate(group=fct_relevel(var,
                           rev(c()))) %>% 
  group_by(group) %>% 
  dplyr::summarise(sum_value=sum(abs(Estimate))) %>% 
  dplyr::mutate(new_col=sum_value/sum(sum_value))
df3
df3 %>% 
  # dplyr::mutate(group=fct_relevel(var,
  # rev(c()))) %>% #求每类变量
  group_by(group) %>% 
  dplyr::summarise(sum_value=sum(abs(Estimate))) %>% 
  dplyr::mutate(new_col=sum_value/sum(sum_value))

df3 %>% 
  # dplyr::mutate(group=fct_relevel(var,
  # rev(c()))) %>% #求每类变量
  group_by(group) %>% 
  dplyr::summarise(sum_value=sum(abs(Estimate))) %>% 
  dplyr::mutate(new_col=sum_value/sum(sum_value))

#cols1=c("#DDE9F7", "#B5D3E9", "#82BADA","#539DCD","#004b79","#436f8e"  "#7E9E8D")
cols1=c( #"#DA561D" , 
  "#9AB560" ,
  "#463D09" ,
  #"#D7D1C2",
  "#fff9ea"
)


p_all2_H=df3 %>% 
  dplyr::mutate(group=fct_relevel(group,
                           rev(c("Stand structural", "Diversity",#"CWM",
                                 "Geographic Location")))) %>% 
  group_by(group) %>% 
  dplyr::summarise(sum_value=sum(abs(Estimate))) %>% 
  dplyr::mutate(new_col=sum_value/sum(sum_value)*100) %>% 
  ggplot(aes(x = 1, y = new_col, label = group)) +
  geom_col(aes(fill = group, color = group), alpha =1, show.legend = FALSE) +
  scale_fill_manual(values = cols1) +
  #scale_color_manual(values = c("#1572CD"  ,"#1572CD" ,"#1572CD" ,"#1572CD"      )    ) + 
  scale_color_manual(values = cols1    ) + 
  coord_flip()+
  
  scale_fill_manual(values = cols1)+
  scale_y_continuous(expand = c(0,0), limits = c(0,100), breaks = c(0,20,40,60,80,100))+
  
  
  theme_minimal()+
  theme(
    axis.text.x = element_blank(),  
    #axis.text.x = element_text(size = 10, family =  "serif"),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),  
    axis.line.y = element_blank(),
    axis.ticks.x = element_blank( ), 
    axis.text.y = element_blank(),  
    plot.title = element_blank()
  )+
  # theme(
  #   axis.text.x = element_text(size = 12, family = "serif"),
  #   axis.title.y = element_blank(),  
  #   axis.line.y = element_blank(),   
  #   axis.line.x = element_line(),
  #   axis.ticks.x = element_line(),
  #   axis.text.y = element_blank(),
  #   plot.title = element_blank(),
  #   axis.title.x = element_text(size = 14, 
  #                               family =  "serif")
  # ) +
  labs(y="Relative effect of estimates (%)")+
  
  labs(subtitle  = expression(Adj.~ italic(R)^2==70~ '%' ))+
  #labs(x="Pinus massoniana")+
  labs(title = expression(High-carbon-sequestration~forests))+
  
  #ggtitle(" ") +
  
  theme(plot.title = element_text(hjust = 0, face = "bold"))+
  
  theme(plot.subtitle = element_text(size = 14, family = "serif") 
  )+
  theme(plot.title = element_text(size = 14, 
                                  family =  "serif"))

p_all2_H

##########################
library(patchwork)
library(patchwork)
pG1=ggarrange(p_all2_H,p_all2_L, p_all2_G ,
              ncol =1,nrow = 3,
              heights = c(1,1),
              widths = c(1,1),common.legend = FALSE, align = c("hv"))
pG1
##################################################################################
####Pinus massoniana：Plot2：100%Plot############
dfP<-read.csv("PD_full_average_All_220.csv",check.names = F)

colnames(dfP) <- c("var" ,"Estimate", "Std. Error","Adjusted.SE", "z.value","Pr(>|z|)"  )

unique(dfP$var)

dfP <- dfP %>%
  dplyr::mutate(across(where(is.numeric), round, 2))


dfP <- dfP %>%
  filter(!( var %in% c("(Intercept)") )  ) %>%
  dplyr::mutate(group = case_when(
    var %in% c("Age","AverageDBH"  ,"DBH_CV"  ,"Density") ~ "Stand structural",
    var %in% c( "PD" ,"pd_Broadleaf"  ,"Shannon" ,"shan_Broadleaf") ~ "Diversity",
    var %in% c("Altitude",
               "MAT" , "MAP"  , "SoilThickness" ,"HumusLayerThickness", "LitterThickness") ~ "Geographic Location",
    TRUE ~ NA_character_  # In case there are any terms not listed above
  ))

cols1=c( 
  "#9AB560" ,
  "#463D09" ,
  #"#D7D1C2",
  "#fff9ea"
)

p_all2_PG=dfP %>% 
  dplyr::mutate(group=fct_relevel(group,
                           rev(c("Stand structural", "Diversity",#"CWM",
                                 "Geographic Location")))) %>% 
  group_by(group) %>% 
  dplyr::summarise(sum_value=sum(abs(Estimate))) %>% 
  dplyr::mutate(new_col=sum_value/sum(sum_value)*100) %>% 
  ggplot(aes(x = 1, y = new_col, label = group)) +
  geom_col(aes(fill = group, color = group), alpha =1, show.legend = FALSE) +
  scale_fill_manual(values = cols1) +
  scale_color_manual(values = cols1    ) +
  #scale_color_manual(values = c("#a1a1a4","#a1a1a4","#a1a1a4" ,"#a1a1a4" )    ) +
  coord_flip()+
  
  scale_fill_manual(values = cols1)+
  scale_y_continuous(expand = c(0,0), limits = c(0,100), breaks = c(0,20,40,60,80,100))+
  
  
  theme_minimal()+
  theme(
    axis.text.x = element_text(size = 12, family = "serif"),
    axis.title.y = element_blank(),  
    axis.line.y = element_blank(),   
    axis.line.x = element_line(),
    axis.ticks.x = element_line(),
    axis.text.y = element_blank(),
    plot.title = element_blank(),
    axis.title.x = element_text(size = 14, 
                                family =  "serif")
  ) +
  
 
  
  labs(y="Relative effect of estimates (%)")+
  
  labs(subtitle  = expression(Adj.~ italic(R)^2==78~ '%' ))+
  #labs(x="Pinus massoniana")+
  labs(title = expression(General))+
  
  #ggtitle(" ") +
  
  theme(plot.title = element_text(hjust = 0, face = "bold")
        
        
  )+
  
  theme(plot.subtitle = element_text(size = 14, family = "serif") 
  )+
  theme(plot.title = element_text(size = 14, 
                                  family =  "serif"))
p_all2_PG
#################
dfP2<-read.csv("PD_full_average_Low_220.csv",check.names = F)

colnames(dfP2) <- c("var" ,"Estimate", "Std. Error","Adjusted.SE", "z.value","Pr(>|z|)"  )

unique(dfP2$var)

dfP2 <- dfP2 %>%
  dplyr::mutate(across(where(is.numeric), round, 2))


dfP2 <- dfP2 %>%
  filter(!( var %in% c("(Intercept)") )  ) %>%
  dplyr::mutate(group = case_when(
    var %in% c("Age","AverageDBH"  ,"DBH_CV"  ,"Density") ~ "Stand structural",
    var %in% c( "PD" ,"pd_Broadleaf"  ,"Shannon" ,"shan_Broadleaf") ~ "Diversity",
    var %in% c("Altitude",
               "MAT" , "MAP"  , "SoilThickness" ,"HumusLayerThickness", "LitterThickness") ~ "Geographic Location",
    TRUE ~ NA_character_  # In case there are any terms not listed above
  ))

cols1=c( 
  "#9AB560" ,
  "#463D09" ,
  #"#D7D1C2",
  "#fff9ea"
)
p_all2_PL=dfP2 %>% 
  dplyr::mutate(group=fct_relevel(group,
                           rev(c("Stand structural", "Diversity",
                                 "Geographic Location")))) %>% 
  group_by(group) %>% 
  dplyr::summarise(sum_value=sum(abs(Estimate))) %>% 
  dplyr::mutate(new_col=sum_value/sum(sum_value)*100) %>% 
  ggplot(aes(x = 1, y = new_col, label = group)) +
  geom_col(aes(fill = group, color = group), alpha =1, show.legend = FALSE) +
  scale_fill_manual(values = cols1) +
  #scale_color_manual(values = cols1    ) + 
  scale_color_manual(values = cols1 ) +
  coord_flip()+
  
  scale_fill_manual(values = cols1)+
  scale_y_continuous(expand = c(0,0), limits = c(0,100), breaks = c(0,20,40,60,80,100))+
  
  
  theme_minimal()+
  theme(
    axis.text.x = element_blank(),  
    #axis.text.x = element_text(size = 10, family =  "serif"),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),  
    axis.line.y = element_blank(),
    axis.ticks.x = element_blank( ), 
    axis.text.y = element_blank(),  
    plot.title = element_blank()
  )+
  
  labs(y="Relative effect of estimates (%)")+
  
  labs(subtitle  = expression(Adj.~ italic(R)^2==92~ '%' ))+
  #labs(x="Pinus massoniana")+
  labs(title = expression(Low-carbon-sequestration~forests))+
  
  #ggtitle(" ") +
  
  theme(plot.title = element_text(hjust = 0, face = "bold"))+
  
  theme(plot.subtitle = element_text(size = 14, family = "serif") 
  )+
  theme(plot.title = element_text(size = 14, 
                                  family =  "serif"))
p_all2_PL
#################
dfP3<-read.csv("PD_full_average_High_220.csv",check.names = F)

colnames(dfP3) <- c("var" ,"Estimate", "Std. Error","Adjusted.SE", "z.value","Pr(>|z|)"  )

unique(dfP3$var)

dfP3 <- dfP3 %>%
  dplyr::mutate(across(where(is.numeric), round, 2))


dfP3 <- dfP3 %>%
  
  filter(!( var %in% c("(Intercept)") )  ) %>%
  dplyr::mutate(group = case_when(
    var %in% c("Age","AverageDBH"  ,"DBH_CV"  ,"Density") ~ "Stand structural",
    var %in% c( "PD" ,"pd_Broadleaf"  ,"Shannon" ,"shan_Broadleaf") ~ "Diversity",
    var %in% c("Altitude",
               "MAT" , "MAP"  , "SoilThickness" ,"HumusLayerThickness", "LitterThickness") ~ "Geographic Location",
    TRUE ~ NA_character_  # In case there are any terms not listed above
  ))




cols1=c( 
  #"#9AB560" ,
  "#463D09" ,
  #"#D7D1C2",
  "#fff9ea"
)

p_all2_PH=dfP3 %>% 
  dplyr::mutate(group=fct_relevel(group,
                           rev(c("Stand structural","Diversity"#, "Geographic Location"
                                 )))) %>% 
  group_by(group) %>% 
  dplyr::summarise(sum_value=sum(abs(Estimate))) %>% 
  dplyr::mutate(new_col=sum_value/sum(sum_value)*100) %>% 
  ggplot(aes(x = 1, y = new_col, label = group)) +
  geom_col(aes(fill = group, color = group), alpha =1, show.legend = FALSE) +
  scale_fill_manual(values = cols1) +
  scale_color_manual(values = cols1) +
  #scale_color_manual(values = c("#1572CD"  ,"#1572CD" ,"#1572CD" ,"#1572CD"      )    ) + 
  coord_flip()+
  
  scale_fill_manual(values = cols1)+
  scale_y_continuous(expand = c(0,0), limits = c(0,100), breaks = c(0,20,40,60,80,100))+
  
  
  theme_minimal()+
  theme(
    axis.text.x = element_blank(),  
    #axis.text.x = element_text(size = 10, family =  "serif"),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),  
    axis.line.y = element_blank(),
    axis.ticks.x = element_blank( ), 
    axis.text.y = element_blank(),  
    plot.title = element_blank()
  )+
  
  labs(y="Relative effect of estimates (%)")+
  
  labs(subtitle  = expression(Adj.~ italic(R)^2==77~ '%' ))+
  #labs(x="Pinus massoniana")+
  labs(title = expression(High-carbon-sequestration~forests))+
  
  #ggtitle(" ") +
  
  theme(plot.title = element_text(hjust = 0, face = "bold"))+
  
  theme(plot.subtitle = element_text(size = 14, family = "serif") 
  )+
  theme(plot.title = element_text(size = 14, 
                                  family =  "serif"))

p_all2_PH

##########################
library(patchwork)
library(patchwork)
pGP1=ggarrange(p_all2_PH,p_all2_PL,p_all2_PG ,
               ncol =1,nrow = 3,
               heights = c(1,1),
               widths = c(1,1),common.legend = FALSE, align = c("hv"))
pGP1

#######################

######ggarrange##########
pG=ggarrange(#p_all1,
  pG1,
  #p_all1_P,
  pGP1,
  ncol =2,nrow =1,
  heights = c(1,1),
  #widths = c(1.6,0.4),
  widths = c(1,1),
  common.legend = FALSE, align = c("hv"))
pG

ggsave("Plot_Liear2_PD_250415.png", path = "E:/MYData/Compare_P_C/figure/Plot_Liear",width =6, height =4.2,dpi=1200, plot=pG)
