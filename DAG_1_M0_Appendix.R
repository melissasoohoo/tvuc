##############################################################
#Purpose: Simulations 
#         This is the 2-time point DAG simulated on a logit scale  
#         This is DAG 1- where U is associated with Y: Model 0

##############################################################


##############################################################
#                                                            #
#               Prepare for Simulation                       #
#               & G-computation                              #
#                                                            #
##############################################################
setwd(" ")

library(dplyr)
library(msm)

rm(list = ls())
cat("\014")

set.seed(8921)
Kreps <-1000
id= as.numeric(1:50000)

#set the parameters
#independent variables 
pz1=0.35
pus=0.25

#intercept
intz =log(0.50)
intx1=log(0.20)
intx2=log(0.25)
inty =log(0.025)

#Relation of var
ORz  =log(1.20)
ORx  =log(1.50)
ORz1y=log(1.20)
ORz2y=log(1.50)

#Relation to Y 
ORxy1=log(1.60)
ORxy2=log(2.00)

#Parameters for U 
ORu1=(log(2))
ORu2=(log(2))

simmed=data.frame(id)
obs =nrow(simmed)

# to store for simulation 
colnames=c("Krep_num", 
          "est_0", "est_x1", "est_x2",
          "est_z1", "est_z2", 
          "est_u1", "est_u2",  
          "SE_0",  "SE_x1", "SE_x2",
          "SE_z1", "SE_z2", 
          "SE_u1", "SE_u2",  
          "P_x1", "P_x2",
          "P_z1", "P_z2",
          "P_u1", "P_u2", 
          "P_Y")

P1_OR_sim <- matrix(NA, nrow=Kreps,ncol=22)
colnames(P1_OR_sim)=colnames

#Store Final Estimates  
colnames0=c("Krep_num", 
            "Y00_est_t","Y10_est_t","Y01_est_t","Y11_est_t",
            "JDE_est_t","dox1_t","dox2_t",
            "Y00_est_b","Y10_est_b","Y01_est_b","Y11_est_b",
            "JDE_est_b","dox1_b","dox2_b")
P1_OR_Est <- matrix(NA, nrow=Kreps,ncol=15)
colnames(P1_OR_Est)=colnames0

#Store True Outcome model 
colnames1=c("Krep_num", "delta_0", 
            "delta_x1", "delta_x2",
            "delta_z1", "delta_z2", 
            "delta_u1", "delta_u2")
delta_mod <- matrix(NA, nrow=Kreps, ncol=8)
colnames(delta_mod)=colnames1

#Store True Intermediate models : Z2
colnames2=c("Krep_num", "z2_beta_0", "z2_beta_x1", "z2_beta_z1")
z2_beta_mod <- matrix(NA, nrow=Kreps, ncol=4)
colnames(z2_beta_mod)=colnames2

#Store Biased Outcome model 
colnames5=c("Krep_num",    "bi_delta_0", 
            "bi_delta_x1", "bi_delta_x2", 
            "bi_delta_z1", "bi_delta_z2" )
bi_delta_mod <- matrix(NA, nrow=Kreps, ncol=6)
colnames(bi_delta_mod)=colnames5

#Store Biased Intermediate models : Z2
colnames6=c("Krep_num", "z2_bi_beta_0", "z2_bi_beta_x1", "z2_bi_beta_z1")
z2_bi_beta_mod <- matrix(NA, nrow=Kreps, ncol=4)
colnames(z2_bi_beta_mod)=colnames6

rm (colnames, colnames0, colnames1, colnames2,  colnames5, colnames6 )

index= 1


##############################################################
#                                                            #
#               Simulation Start:                            #
#               & G-computation models                       #
#                                                            #
##############################################################


for(j in 1:Kreps  ){ 
  
  simmed$z1= rbinom(obs, 1, pz1)
  
  simmed$u1= rbinom(obs, 1, pus)
  simmed$u2= rbinom(obs, 1, pus)

  simmed$x1= rbinom(obs, 1, 
                      plogis(intx1 + (ORz*simmed$z1) + (ORu1*simmed$u1)))
  
  simmed$z2= rbinom(obs, 1, 
                      plogis(intz + (ORz*simmed$z1) + (ORx*simmed$x1)))
  
  simmed$x2= rbinom(obs, 1, 
                      plogis(intx2 + (ORx*simmed$x1) + (ORz*simmed$z1)+
                              (ORz*simmed$z2) + (ORu2*simmed$u2))) 

  simmed$y=  rbinom(obs, 1,   
                      plogis(inty + 
                               (ORxy1*simmed$x1) + (ORxy2*simmed$x2) +
                               (ORz1y*simmed$z1) + (ORz2y*simmed$z2) +
                               (ORu1*simmed$u1) + (ORu2*simmed$u2)  ))

  #Check the convergence of the model:
  esty_mod <- glm(y~x1 + x2 + z1 + z2  + u1 + u2 , 
                  data=simmed,
                  family = binomial (link="logit"))
  
  vcv <- vcov(esty_mod)
  
  P1_OR_sim[j,1]=j  
  P1_OR_sim[j,2]=(coef(esty_mod)[1])  #int
  P1_OR_sim[j,3]=(coef(esty_mod)[2])  #x1
  P1_OR_sim[j,4]=(coef(esty_mod)[3])  #x2
  P1_OR_sim[j,5]=(coef(esty_mod)[4])  #z1
  P1_OR_sim[j,6]=(coef(esty_mod)[5])  #z2
  P1_OR_sim[j,7]=(coef(esty_mod)[6])  #u1
  P1_OR_sim[j,8]=(coef(esty_mod)[7])  #u2

  P1_OR_sim[j,9]=sqrt(diag(vcv)[1])  #SE of intercept
  P1_OR_sim[j,10]=sqrt(diag(vcv)[2]) #SE of estx1
  P1_OR_sim[j,11]=sqrt(diag(vcv)[3]) #SE of estx2
  P1_OR_sim[j,12]=sqrt(diag(vcv)[4]) #SE of estz1
  P1_OR_sim[j,13]=sqrt(diag(vcv)[5]) #SE of estz2
  P1_OR_sim[j,14]=sqrt(diag(vcv)[6]) #SE of estu1
  P1_OR_sim[j,15]=sqrt(diag(vcv)[7]) #SE of estu2

  P1_OR_sim[j,16]=mean(simmed$x1) #p(x1)
  P1_OR_sim[j,17]=mean(simmed$x2) #p(x2)
  P1_OR_sim[j,18]=mean(simmed$z1) #p(z1)
  P1_OR_sim[j,19]=mean(simmed$z2) #p(z2)
  P1_OR_sim[j,20]=mean(simmed$u1) #p(u1)
  P1_OR_sim[j,21]=mean(simmed$u2) #p(u2)
  P1_OR_sim[j,22]=mean(simmed$y)  #p(y)
  

  ##############################################################
  #                                                            #
  #              True Model, if we had U                       #
  #                                                            #
  ##############################################################
  
  delta=glm(y ~ x1 + x2 + z1 + z2  + u1 + u2  , 
            data= simmed, 
            family= binomial (link="logit"))
  
  delta_mod[j,1]=j  
  delta_mod[j,2]=(coef(delta)[1])#bo
  delta_mod[j,3]=(coef(delta)[2])#x1
  delta_mod[j,4]=(coef(delta)[3])#x2
  delta_mod[j,5]=(coef(delta)[4])#z1
  delta_mod[j,6]=(coef(delta)[5])#z1
  delta_mod[j,7]=(coef(delta)[6])#u1
  delta_mod[j,8]=(coef(delta)[7])#u2
  
  #modeling the descendants - for the true model: 
  z2_beta = glm(z2 ~ x1 + z1 ,
                data= simmed, 
                family= binomial (link="logit"))

  z2_beta_mod[j,1]=j  
  z2_beta_mod[j,2]=(coef(z2_beta)[1])#int
  z2_beta_mod[j,3]=(coef(z2_beta)[2])#x1
  z2_beta_mod[j,4]=(coef(z2_beta)[3])#z1
  
  ##############################################################
  #                                                            #
  #              Biased Model, if we didn't have U             #
  #                                                            #
  ##############################################################
  
  bi_delta=glm(y ~ x1 + x2 + z1 + z2 ,
              data=simmed, family= binomial (link="logit"))
  
  bi_delta_mod[j,1]=j  
  bi_delta_mod[j,2]=(coef(bi_delta)[1]) #int
  bi_delta_mod[j,3]=(coef(bi_delta)[2]) #x1
  bi_delta_mod[j,4]=(coef(bi_delta)[3]) #x2
  bi_delta_mod[j,5]=(coef(bi_delta)[4]) #z1
  bi_delta_mod[j,6]=(coef(bi_delta)[5]) #z2

  #modeling the descendants - for the biased model: 
  z2_bi_beta = glm(z2 ~ x1 + z1 ,
                   data= simmed, 
                   family= binomial (link="logit"))

  z2_bi_beta_mod[j,1]=j  
  z2_bi_beta_mod[j,2]=(coef(z2_bi_beta)[1])#int
  z2_bi_beta_mod[j,3]=(coef(z2_bi_beta)[2])#x1
  z2_bi_beta_mod[j,4]=(coef(z2_bi_beta)[3])#z1
  

  ##############################################################
  #                                                            #
  #               Calculate the PO:                            #
  #                                                            #
  ##############################################################
  
  #Z2 for the True Model
  
  simmed$z2_0 = 1/(1+exp(-(z2_beta_mod[j,2] + z2_beta_mod[j,3]*0 + z2_beta_mod[j,4]*simmed$z1)))
  simmed$z2_1 = 1/(1+exp(-(z2_beta_mod[j,2] + z2_beta_mod[j,3]*1 + z2_beta_mod[j,4]*simmed$z1)))

  #True Model PO
  
  simmed$y00_t = 1/(1+exp(-(delta_mod[j,2] + delta_mod[j,3]*0 + delta_mod[j,4]*0 +
                            delta_mod[j,5]*simmed$z1 + delta_mod[j,6]*simmed$z2_0 + 
                            delta_mod[j,7]*simmed$u1 + delta_mod[j,8]*simmed$u2 ))) 

  simmed$y10_t = 1/(1+exp(-(delta_mod[j,2] + delta_mod[j,3]*1 + delta_mod[j,4]*0 +
                            delta_mod[j,5]*simmed$z1 + delta_mod[j,6]*simmed$z2_1 + 
                            delta_mod[j,7]*simmed$u1 + delta_mod[j,8]*simmed$u2 ))) 

  simmed$y01_t = 1/(1+exp(-(delta_mod[j,2] + delta_mod[j,3]*0 + delta_mod[j,4]*1 +
                            delta_mod[j,5]*simmed$z1 + delta_mod[j,6]*simmed$z2_0 + 
                            delta_mod[j,7]*simmed$u1 + delta_mod[j,8]*simmed$u2 ))) 

  simmed$y11_t = 1/(1+exp(-(delta_mod[j,2] + delta_mod[j,3]*1 + delta_mod[j,4]*1 +
                            delta_mod[j,5]*simmed$z1 + delta_mod[j,6]*simmed$z2_1 + 
                            delta_mod[j,7]*simmed$u1 + delta_mod[j,8]*simmed$u2 ))) 

simmed$y_jde_t   =  (simmed$y11_t/(1-simmed$y11_t)) /(simmed$y00_t/(1-simmed$y00_t))
  
  #Calculate the effect estimate 
  simmed$y_dox1_t = (simmed$y10_t/(1-simmed$y10_t)) /(simmed$y00_t/(1-simmed$y00_t))
  simmed$y_dox2_t = (simmed$y01_t/(1-simmed$y01_t)) /(simmed$y00_t/(1-simmed$y00_t))
  
  #Z2 for the Bias Model
  
  simmed$z2_bi_0 = 1/(1+exp(-(z2_bi_beta_mod[j,2] + z2_bi_beta_mod[j,3]*0 + z2_bi_beta_mod[j,4]*simmed$z1)))
  simmed$z2_bi_1 = 1/(1+exp(-(z2_bi_beta_mod[j,2] + z2_bi_beta_mod[j,3]*1 + z2_bi_beta_mod[j,4]*simmed$z1)))
  
  #Bias Model PO
  
  simmed$y00_b = 1/(1+exp(-(bi_delta_mod[j,2] + bi_delta_mod[j,3]*0 + bi_delta_mod[j,4]*0 +
                                bi_delta_mod[j,5]*simmed$z1 + bi_delta_mod[j,6]*simmed$z2_bi_0))) 
  
  simmed$y10_b = 1/(1+exp(-(bi_delta_mod[j,2] + bi_delta_mod[j,3]*1 + bi_delta_mod[j,4]*0 +
                              bi_delta_mod[j,5]*simmed$z1 + bi_delta_mod[j,6]*simmed$z2_bi_1))) 
  
  simmed$y01_b = 1/(1+exp(-(bi_delta_mod[j,2] + bi_delta_mod[j,3]*0 + bi_delta_mod[j,4]*1 +
                              bi_delta_mod[j,5]*simmed$z1 + bi_delta_mod[j,6]*simmed$z2_bi_0))) 
  
  simmed$y11_b = 1/(1+exp(-(bi_delta_mod[j,2] + bi_delta_mod[j,3]*1 + bi_delta_mod[j,4]*1 +
                              bi_delta_mod[j,5]*simmed$z1 + bi_delta_mod[j,6]*simmed$z2_bi_1))) 
  
  simmed$y_jde_b   =  (simmed$y11_b/(1-simmed$y11_b)) /(simmed$y00_b/(1-simmed$y00_b))
  
  #Compute each estimate 
  simmed$y_dox1_b = (simmed$y10_b/(1-simmed$y10_b)) /(simmed$y00_b/(1-simmed$y00_b))
  simmed$y_dox2_b = (simmed$y01_b/(1-simmed$y01_b)) /(simmed$y00_b/(1-simmed$y00_b))

  #save all mean estimates per K replication
  
  P1_OR_Est[j,1]=j
  P1_OR_Est[j,2]=mean(simmed$y00_t)
  P1_OR_Est[j,3]=mean(simmed$y10_t)
  P1_OR_Est[j,4]=mean(simmed$y01_t)
  P1_OR_Est[j,5]=mean(simmed$y11_t)
  P1_OR_Est[j,6]=mean(simmed$y_jde_t)
  P1_OR_Est[j,7]=mean(simmed$y_dox1_t)
  P1_OR_Est[j,8]=mean(simmed$y_dox2_t)
  
  P1_OR_Est[j,9]=mean(simmed$y00_b)
  P1_OR_Est[j,10]=mean(simmed$y10_b)
  P1_OR_Est[j,11]=mean(simmed$y01_b)
  P1_OR_Est[j,12]=mean(simmed$y11_b)
  P1_OR_Est[j,13]=mean(simmed$y_jde_b)
  P1_OR_Est[j,14]=mean(simmed$y_dox1_b)
  P1_OR_Est[j,15]=mean(simmed$y_dox2_b)

  cat(j)
  index=index+1
  
}

##############################################################
#                                                            #
#               save and summarize estimates                 #
#               across all K=1000 reps                       #
##############################################################

Y00_T_est=median(P1_OR_Est[,2])
Y00_T_LCL=quantile(P1_OR_Est[,2], c(0.025))
Y00_T_UCL=quantile(P1_OR_Est[,2], c(0.975))

Y10_T_est=median(P1_OR_Est[,3])
Y10_T_LCL=quantile(P1_OR_Est[,3], c(0.025))
Y10_T_UCL=quantile(P1_OR_Est[,3], c(0.975))

Y01_T_est=median(P1_OR_Est[,4])
Y01_T_LCL=quantile(P1_OR_Est[,4], c(0.025))
Y01_T_UCL=quantile(P1_OR_Est[,4], c(0.975))

Y11_T_est=median(P1_OR_Est[,5])
Y11_T_LCL=quantile(P1_OR_Est[,5], c(0.025))
Y11_T_UCL=quantile(P1_OR_Est[,5], c(0.975))

OR_True_P50=median(P1_OR_Est[,6])
OR_True_LCL=quantile(P1_OR_Est[,6], c(0.025))
OR_True_UCL=quantile(P1_OR_Est[,6], c(0.975))

YDOX1_T_est=median(P1_OR_Est[,7])
YDOX1_T_LCL=quantile(P1_OR_Est[,7], c(0.025))
YDOX1_T_UCL=quantile(P1_OR_Est[,7], c(0.975))

YDOX2_T_est=median(P1_OR_Est[,8])
YDOX2_T_LCL=quantile(P1_OR_Est[,8], c(0.025))
YDOX2_T_UCL=quantile(P1_OR_Est[,8], c(0.975))

#Biased
Y00_B_est=median(P1_OR_Est[,9])
Y00_B_LCL=quantile(P1_OR_Est[,9], c(0.025))
Y00_B_UCL=quantile(P1_OR_Est[,9], c(0.975))

Y10_B_est=median(P1_OR_Est[,10])
Y10_B_LCL=quantile(P1_OR_Est[,10], c(0.025))
Y10_B_UCL=quantile(P1_OR_Est[,10], c(0.975))

Y01_B_est=median(P1_OR_Est[,11])
Y01_B_LCL=quantile(P1_OR_Est[,11], c(0.025))
Y01_B_UCL=quantile(P1_OR_Est[,11], c(0.975))

Y11_B_est=median(P1_OR_Est[,12])
Y11_B_LCL=quantile(P1_OR_Est[,12], c(0.025))
Y11_B_UCL=quantile(P1_OR_Est[,12], c(0.975))

OR_Bias_P50=median(P1_OR_Est[,13])
OR_Bias_LCL=quantile(P1_OR_Est[,13], c(0.025))
OR_Bias_UCL=quantile(P1_OR_Est[,13], c(0.975))

YDOX1_B_est=median(P1_OR_Est[,14])
YDOX1_B_LCL=quantile(P1_OR_Est[,14], c(0.025))
YDOX1_B_UCL=quantile(P1_OR_Est[,14], c(0.975))

YDOX2_B_est=median(P1_OR_Est[,15])
YDOX2_B_LCL=quantile(P1_OR_Est[,15], c(0.025))
YDOX2_B_UCL=quantile(P1_OR_Est[,15], c(0.975))

fin_OR_True=data.frame("OR2T", "DAG 1 MODEL 0", Y00_T_est,Y00_T_LCL,Y00_T_UCL,
                                        Y10_T_est,Y10_T_LCL,Y10_T_UCL,
                                        Y01_T_est,Y01_T_LCL,Y01_T_UCL,
                                        Y11_T_est,Y11_T_LCL,Y11_T_UCL,
                                        OR_True_P50,OR_True_LCL,OR_True_UCL,
                                        YDOX1_T_est,YDOX1_T_LCL,YDOX1_T_UCL, 
                                        YDOX2_T_est,YDOX2_T_LCL,YDOX2_T_UCL
                       )

fin_OR_Bias=data.frame("OR2B", "DAG 1 MODEL 0", Y00_B_est,Y00_B_LCL,Y00_B_UCL,
                                        Y10_B_est,Y10_B_LCL,Y10_B_UCL,
                                        Y01_B_est,Y01_B_LCL,Y01_B_UCL,
                                        Y11_B_est,Y11_B_LCL,Y11_B_UCL,
                                        OR_Bias_P50,OR_Bias_LCL,OR_Bias_UCL,
                                        YDOX1_B_est,YDOX1_B_LCL,YDOX1_B_UCL, 
                                        YDOX2_B_est,YDOX2_B_LCL,YDOX2_B_UCL
                       )

#QAQC on the simulation

ab.est0 <- mean(abs(P1_OR_sim[ , 2] -  (inty)))
ab.estx1 <- mean(abs(P1_OR_sim[ , 3] -  (ORxy1)))
ab.estx2 <- mean(abs(P1_OR_sim[ , 4] -  (ORxy2)))
ab.estz1 <- mean(abs(P1_OR_sim[ , 5] -  (ORz1y)))
ab.estz2 <- mean(abs(P1_OR_sim[ , 6] -  (ORz2y))) 
ab.estu1 <- mean(abs(P1_OR_sim[ , 7] -  (ORu1))) 
ab.estu2 <- mean(abs(P1_OR_sim[ , 8] -  (ORu2))) 

mse.est0  <- mean((P1_OR_sim[ , 2] -  (inty))^2)
mse.estx1 <- mean((P1_OR_sim[ , 3] -  (ORxy1))^2) 
mse.estx2 <- mean((P1_OR_sim[ , 4] -  (ORxy2))^2) 
mse.estz1 <- mean((P1_OR_sim[ , 5] -  (ORz1y))^2)
mse.estz2 <- mean((P1_OR_sim[ , 6] -  (ORz2y))^2)  
mse.estu1 <- mean((P1_OR_sim[ , 7] -  (ORu1))^2) 
mse.estu2 <- mean((P1_OR_sim[ , 8] -  (ORu2))^2) 

p_x1 =mean(P1_OR_sim[,16])
p_x2 =mean(P1_OR_sim[,17])
p_z1 =mean(P1_OR_sim[,18])
p_z2 =mean(P1_OR_sim[,19])
p_u1 =mean(P1_OR_sim[,20])
p_u2 =mean(P1_OR_sim[,21])
p_y =mean(P1_OR_sim[,22])

absdta=data.frame("OR2T", "DAG 1 MODEL 0", ab.est0,ab.estx1,ab.estx2,
                  ab.estz1,ab.estz2,ab.estu1,ab.estu2)
msedta=data.frame("OR2T", "DAG 1 MODEL 0", mse.est0,mse.estx1,mse.estx2,
                  mse.estz1,mse.estz2,mse.estu1,mse.estu2)

pdta=data.frame("OR2", "DAG 1 MODEL 0", p_x1,p_x2,
                p_z1,p_z2,p_u1,p_u2,p_y)

# Merge columns and SAVE 
P1_OR_Est_DAG_1_M0=P1_OR_Est

#save the estimates as a Rdata file for future reference 
save(P1_OR_Est_DAG_1_M0, file=" .rda")

#summarize the estimates into a Text File

library(knitr)
sink(" .txt")

#Final estimates 

knitr::kable("DAG 1 MODEL 0")
knitr::kable(fin_OR_True,digits=5)
knitr::kable(fin_OR_Bias,digits=5)

knitr::kable(absdta,digits=5)
knitr::kable(msedta,digits=5)

knitr::kable(pdta,digits=5)

sink()


