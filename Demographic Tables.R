library(lsr)
library(effsize)
library(confintr)

###############Table 1 Deomgraphics#####################

#Bringing in phenotype data
data = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/aim2analysis_final2.csv', header=T)
data2 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/OldFiles/aim2analysis_demo.csv', header=T)

test=merge(data,data2)
sbp=test[na.omit(test$SBP_new),]

table(data$gender)
prop.table(table(data$gender))

table(data$race)
prop.table(table(data$race))

table(data$edu)
prop.table(table(data$edu))

table(data$smoking)
prop.table(table(data$smoking))

mean(data$age)
sd(data$age)


#Medications
table(data$bpmeds)
prop.table(table(data$bpmeds))
sum(is.na(data$bpmeds))

table(data$diabmeds)
prop.table(table(data$diabmeds))
sum(is.na(data$diabmeds))

table(data$cholmeds)
prop.table(table(data$cholmeds))
sum(is.na(data$cholmeds))


#Cardiometabolic Risk Factors
mean(data$BMI_impute,na.rm=T)
sd(data$BMI_impute, na.rm=T)
sum(is.na(data$BMI_impute))

mean(data$WC_new,na.rm=T)
sd(data$WC_new,na.rm=T)
sum(is.na(data$WC_new))

mean(data$SBP,na.rm=T)
sd(data$SBP,na.rm=T)
sum(is.na(data$SBP))

mean(data$DBP,na.rm=T)
sd(data$DBP,na.rm=T)
sum(is.na(data$DBP))

mean(data$PGLUFF,na.rm=T)
sd(data$PGLUFF,na.rm=T)
sum(is.na(data$PGLUFF))

mean(data2$PCHOL,na.rm=T)
sd(data2$PCHOL,na.rm=T)
sum(is.na(data2$PCHOL))

mean(test$PLDLC,na.rm=T)
sd(test$PLDLC,na.rm=T)
sum(is.na(test$PLDLC))

mean(data$PHDLD,na.rm=T)
sd(data$PHDLD,na.rm=T)
sum(is.na(data$PHDLD))

mean(data$PCRP,na.rm=T)
sd(data$PCRP,na.rm=T)
sum(is.na(data$PCRP))


mean(data2$TC_new,na.rm=T)
sd(data2$TC_new,na.rm=T)
sum(is.na(data2$TC_new))

mean(data2$LDL_new,na.rm=T)
sd(data2$LDL_new,na.rm=T)
sum(is.na(data2$LDL_new))

mean(test$PCHOL,na.rm=T)
sd(test$PCHOL,na.rm=T)
sum(is.na(test$PCHOL))

mean(test$DBP_new,na.rm=T)
sd(test$DBP_new,na.rm=T)
sum(is.na(test$DBP_new))


###################Supplemental Table 2 - Comparing participants included vs.excluded in primary analysis from DNA methylation sample##### 
###Bringing in full DNA methylation sample and finsl analysis with demographics
fullsample = read.csv('//orion.sph.umich.edu/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/missing.csv', header=T)
ptype=read.csv('//orion.sph.umich.edu/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/aim2analysis_final2.csv', header=T)
kept=ptype[,1:2]

###Merging datasests and creating an indicated for included vs. excluded 
mergedsample=merge(fullsample, kept, by="id", all.x=T)
mergedsample$excluded=ifelse(mergedsample$indicator2010 >= 0,0,1) 
mergedsample$excluded[is.na(mergedsample$excluded)] = 1 


#Creating seaparate datasets for included vs. excluded participants
included=mergedsample[mergedsample$excluded==0,]
excluded=mergedsample[mergedsample$excluded==1,]

#Age - means, SDs, t-test, and Cohen's D
mean(excluded$age)
sd(excluded$age)
mean(included$age)
sd(included$age)
t.test(age~excluded, alternative = c("two.sided"),data=mergedsample)
cohensD(included$age, excluded$age)

#Gender - N(%), chi-sq test, and Cramers's V
table(excluded$gender)
prop.table(table(excluded$gender))
table(included$gender)
prop.table(table(included$gender))
chisq.test(mergedsample$gender, mergedsample$excluded, correct=FALSE)
gender = matrix(c(774,576,1575,1093), nrow = 2)
cramersv(gender)

#Race - N(%), chi-sq test, and Cramers's V
table(excluded$race)
prop.table(table(excluded$race))
table(included$race)
prop.table(table(included$race))
chisq.test(mergedsample$race, mergedsample$excluded,correct=FALSE)
race = matrix(c(298,666,309,54,269,1981,349,68), nrow = 4)
cramersv(race)

table(excluded$edu)
prop.table(table(excluded$edu))
table(included$edu)
prop.table(table(included$edu))
chisq.test(mergedsample$edu, mergedsample$excluded, correct=FALSE)
edu = matrix(c(329,758,262,346,1617,705), nrow = 3)
chisq.test(edu, correct=FALSE)
cramersv(edu)

#Smoking - N(%), chi-sq test, and Cramers's V
table(excluded$smoking)
prop.table(table(excluded$smoking))
table(included$smoking)
prop.table(table(included$smoking))
smoking = matrix(c(540,527,218,1207,1134,327), nrow = 3)
chisq.test(missing$smoking, missing$missing, correct=FALSE)
chisq.test(smoking, correct=FALSE)
cramersv(smoking)


######EXTRA

table(miss$workpay)
prop.table(table(miss$workpay))
table(there$workpay)
prop.table(table(there$workpay))
workpay = matrix(c(365,917,1514,1175), nrow = 2)
chisq.test(missing$workpay, missing$missing, correct=FALSE)
chisq.test(workpay, correct=FALSE)
cramersv(workpay)


table(miss$partner_spouse)
prop.table(table(miss$partner_spouse))
table(there$partner_spouse)
prop.table(table(there$partner_spouse))
chisq.test(missing$partner_spouse, missing$missing, correct=FALSE)
partner = matrix(c(473,812,784,1905), nrow = 2)
chisq.test(partner, correct=FALSE)
cramersv(partner)

table(miss$haschild)
prop.table(table(miss$haschild))
table(there$haschild)
prop.table(table(there$haschild))
chisq.test(missing$haschild, missing$missing, correct=FALSE)
child = matrix(c(159,1170,298,2390), nrow = 2)
chisq.test(child, correct=FALSE)
cramersv(partner)





