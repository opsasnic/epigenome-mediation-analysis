ptype=read.csv('//orion.sph.umich.edu/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/aim2analysis_final2.csv', header=T)
trig=read.csv('//orion.sph.umich.edu/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/trig.csv', header=T)
trig2=read.csv('//orion.sph.umich.edu/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/trig2.csv', header=T)



#aim2 analysis plan 
ptype$crp3=log(ptype$PCRP)
ptype$glucose=log(ptype$PGLUFF2)


xwalk = read.csv('//orion.sph.umich.edu/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Methylation_Apr_2021/transfer2/xwalk.csv', header=T)
plateinfo = read.csv('//orion.sph.umich.edu/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Methylation_Apr_2021/transfer04.13.2021/finalPD.csv', header=T)[,c(23, 3, 4, 9, 8)]
cells = read.csv('//orion.sph.umich.edu/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Methylation_Apr_2021/transfer2/DNAm_cellpropest_HRSn4018.csv', header=T)
info = merge(merge(xwalk, plateinfo, by='FID'), cells, by.x='id', by.y='sample')
info$row = as.factor(substr(info$BCD_Well, 1, 1))
info$col = as.factor(substr(info$BCD_Well, 2, 3))

ptype2 = merge(ptype, info[,-c(4,5)], by='id')
ptype2$Slide = as.factor(ptype2$Slide)

#PCs
PC.comb = read.csv('//orion.sph.umich.edu/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Genotype_Jan2015_Processing/PCA/Phase1234/Top_PCs.csv')
colnames(PC.comb) = c('local_id', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')

data = merge(ptype2, PC.comb, by='local_id')

data2 = merge(data, trig, by='id')
data3 = merge(data, trig2, by='id')

data3$trig=log(data3$PTGF2)



nomeds=data[data$diabmeds==0,]
nostatin=data[data$cholmeds==0,]
fasting=data[data$PFASTYN==1,]
fasting2=data2[data2$PFASTYN==1,]
fasting3=data3[data3$PFASTYN==1,]
fastnomeds=nomeds[nomeds$PFASTYN==1,]
nostatin=data3[data3$cholmeds==0,]


#Total Cholesterol
mod1a=lm(TC_new ~ cumbur6_zscore2 + gender + age +  edu  + indicator2010 + PFASTYN + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10, data=data)
mod1b=lm(TC_new~PFASTYN+cumbur6_zscore2  + age + gender + edu + indicator2010 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + smoking , data=data)

summary(mod1a)
summary(mod1b)

#LDL
mod2a=lm(LDL_new ~ cumbur6_zscore2 + gender + age +  edu  + indicator2010 + PFASTYN + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10, data=data)

mod2b=lm(LDL_new~PFASTYN+cumbur6_zscore2  + age + gender + edu + indicator2010 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + smoking , data=data)
mod2c=lm(LDL_new~PFASTYN+cumbur6_zscore  + age + gender + edu + indicator2010 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + smoking + BMI.x , data=data)
mod2c=lm(LDL_new~cumbur6_zscore  + age + gender + edu + indicator2010 + PFASTYN + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + smoking , data=data)

summary(mod2a)
summary(mod2b)
summary(mod2c)

#HDL
mod3a=lm(PHDLD ~ cumbur6_zscore2 + gender + age +  edu  + indicator2010 + PFASTYN + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10, data=data)
mod3b=lm(PHDLD~PFASTYN+cumbur6_zscore2  + age + gender + edu + indicator2010 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + smoking , data=data)
mod3c=lm(PHDLD ~ cumbur6_zscore + gender + age +  edu  + indicator2010 + PFASTYN + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10 + smoking + BMI.x, data=data)
mod3c=lm(PHDLD ~ cumbur6_zscore + gender + age +  edu  + indicator2010 + PFASTYN + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10 + smoking, data=data)

summary(mod3a)
summary(mod3b)
summary(mod3c)

#CRP
mod4a=lm(crp3 ~ cumbur6_zscore2 + gender + age +  edu  + indicator2010 + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10, data=data)
mod4b=lm(crp3~cumbur6_zscore2  + age + gender + edu + indicator2010 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10+smoking, data=data)
summary(mod4a) 
summary(mod4b) 
summary(mod4c) 

#Glucose
mod5a=lm(glucose ~ cumbur6_zscore2 + gender + age +  edu  + indicator2010 + PFASTYN + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10, data=nomeds)
mod5b=lm(glucose~cumbur6_zscore2 + PFASTYN + age + gender + edu + indicator2010 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +smoking, data=nomeds)
mod5c=lm(glucose ~ cumbur6_zscore + gender + age +  edu  + indicator2010 + PFASTYN + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10 + smoking + BMI.x, data=nomeds)
mod5c=lm(glucose ~ cumbur6_zscore2 + gender + age +  edu  + indicator2010 + PFASTYN + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10 + smoking, data=nomeds)

summary(mod5a)
summary(mod5b)
summary(mod5c)


#SBP
mod6a=lm(SBP_new ~ cumbur6_zscore + gender + age +  edu  + indicator2010 + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10, data=data)

mod6b=lm(SBP_new~cumbur6_zscore  + age + gender + edu + indicator2010 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10+smoking, data=data)
mod6c=lm(SBP_new ~ cumbur6_zscore + gender + age +  edu  + indicator2010 + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10 + smoking + BMI.x, data=data)
mod6c=lm(SBP_new ~ cumbur6_zscore + gender + age +  edu  + indicator2010 + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10 + smoking, data=data)

summary(mod6a)
summary(mod6b)
summary(mod6c)


#DBP
mod7=lm(DBP_new ~ cumbur6_zscore + gender + age +  edu  + indicator2010 + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10, data=data)

mod7b=lm(DBP_new~cumbur6_zscore  + age + gender + edu + indicator2010 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +smoking, data=data)
mod7c=lm(DBP_new~cumbur6_zscore  + age + gender + edu + indicator2010 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + BMI.x +smoking, data=data)
mod7c=lm(DBP_new ~ cumbur6_zscore + gender + age +  edu  + indicator2010 + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10 + smoking + BMI.x, data=data)
mod7c=lm(DBP_new ~ cumbur6_zscore + gender + age +  edu  + indicator2010 + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10 + smoking, data=data)

summary(mod7)
summary(mod7a)
summary(mod7c)


#BMI
mod8a=lm(BMI_impute ~ cumbur6_zscore + gender + age +  edu  + indicator2010 + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10, data=data)
mod8b=lm(BMI_impute ~ cumbur6_zscore + gender + age +  edu  + indicator2010 + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10 + smoking, data=data)

summary(mod8a)
summary(mod8b)


#WC
mod9a=lm(WC_new ~ cumbur6_zscore2 + gender + age +  edu  + indicator2010 + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10, data=data)

mod9b=lm(WC_new~cumbur6_zscore2  + age + gender + edu + indicator2010 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + smoking, data=data)
mod9c=lm(WC_new ~ cumbur6_zscore + gender + age +  edu  + indicator2010 + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10 + smoking, data=data)


summary(mod9a)
summary(mod9b)
summary(mod9c)



#Fasting status 
mod1=lm(TC_new ~ cumbur6_zscore2 + gender + age +  edu  + indicator2010  + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10, data=fasting3)
mod2=lm(PHDLD ~ cumbur6_zscore2 + gender + age +  edu  + indicator2010 + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10, data=fasting3)
mod3=lm(LDL_new ~ cumbur6_zscore2 + gender + age +  edu  + indicator2010 + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10, data=fasting3)
mod4=lm(glucose ~ cumbur6_zscore2 + gender + age +  edu  + indicator2010 + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10, data=fastnomeds)
mod5=lm(PTGF2 ~ cumbur6_zscore2 + gender + age +  edu  + indicator2010 + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10, data=fasting3)

summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)

#Fasting status - smoking 
mod1=lm(TC_new ~ cumbur6_zscore2 + gender + age +  edu  + indicator2010  + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10 + smoking, data=fasting3)
mod2=lm(PHDLD ~ cumbur6_zscore2 + gender + age +  edu  + indicator2010 + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10 + smoking, data=fasting3)
mod3=lm(LDL_new ~ cumbur6_zscore2 + gender + age +  edu  + indicator2010 + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10+ smoking, data=fasting3)
mod4=lm(glucose ~ cumbur6_zscore2 + gender + age +  edu  + indicator2010 + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10 + smoking, data=fastnomeds)
mod5=lm(PTGF2 ~ cumbur6_zscore2 + gender + age +  edu  + indicator2010 + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10 + smoking, data=fasting3)

summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)

mean(data$PCRP, na.rm=TRUE)
sd(nomeds$PGLUFF, na.rm=TRUE)
mean(data$PGLUFF, na.rm=TRUE)
sd(data$TC_new, na.rm=TRUE)
sd(data$TC_new, na.rm=TRUE)

sd(nomeds$PGLUFF, na.rm=TRUE)
mean(data$PGLUFF, na.rm=TRUE)

sd(data$PGLUFF2, na.rm=TRUE)
mean(data$PGLUFF2, na.rm=TRUE)

sd(nomeds$PGLUFF2, na.rm=TRUE)
mean(nomeds$PGLUFF2, na.rm=TRUE)



sum(is.na(nomeds$PGLUFF2))
sum(is.na(data$PGLUFF2))

