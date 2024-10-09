####Read in phenotype file comtaining psychosocial stress measure, demographic covariates, and cardiometabolic risk factors#####
ptype=read.csv('//orion.sph.umich.edu/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/aim2analysis_final2.csv', header=T)
trig2=read.csv('//orion.sph.umich.edu/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/trig2.csv', header=T)
ptype2 = merge(ptype, trig2, by='id')


###Log transform CRP, glucose to achieve a more normal distribution 
ptype2$crp3=log(ptype2$PCRP)
ptype2$glucose=log(ptype2$PGLUFF2)

###reading in the crosswalk file for later merge 
xwalk = read.csv('//orion.sph.umich.edu/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Methylation_Apr_2021/transfer2/xwalk.csv', header=T)
plateinfo = read.csv('//orion.sph.umich.edu/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Methylation_Apr_2021/transfer04.13.2021/finalPD.csv', header=T)[,c(23, 3, 4, 9, 8)]
cells = read.csv('//orion.sph.umich.edu/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Methylation_Apr_2021/transfer2/DNAm_cellpropest_HRSn4018.csv', header=T)
info = merge(merge(xwalk, plateinfo, by='FID'), cells, by.x='id', by.y='sample')
info$row = as.factor(substr(info$BCD_Well, 1, 1))
info$col = as.factor(substr(info$BCD_Well, 2, 3))
ptype3 = merge(ptype2, info[,-c(4,5)], by='id')
ptype3$Slide = as.factor(ptype3$Slide)

###reading in the top 10 genetic ancestry PCs that are covariates in the models
PC.comb = read.csv('//orion.sph.umich.edu/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Genotype_Jan2015_Processing/PCA/Phase1234/Top_PCs.csv')
colnames(PC.comb) = c('local_id', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')

toteffectsdata = merge(ptype3, PC.comb, by='local_id')


##Restricting the data to individuals who do not take diabetes medication - use  this dataset for glucose associations 
nomeds=toteffectsdata[toteffectsdata$diabmeds==0,]

###Total effect models with cumulative psychosocial stress (cumbur6_zscore2) and releveant cardiometabolic risk factor
###Model A is the primary model and Model B includes the addition of smoking status as a covarate

#Total Cholesterol
mod1a=lm(TC_new ~ cumbur6_zscore2 + gender + age +  edu  + indicator2010 + PFASTYN + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10, data=toteffectsdata)
mod1b=lm(TC_new~PFASTYN+cumbur6_zscore2  + age + gender + edu + indicator2010 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + smoking , data=toteffectsdata)

summary(mod1a)
summary(mod1b)

#LDL
mod2a=lm(LDL_new ~ cumbur6_zscore2 + gender + age +  edu  + indicator2010 + PFASTYN + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10, data=toteffectsdata)
mod2b=lm(LDL_new~PFASTYN+cumbur6_zscore2  + age + gender + edu + indicator2010 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + smoking , data=toteffectsdata)

summary(mod2a)
summary(mod2b)

#HDL
mod3a=lm(PHDLD ~ cumbur6_zscore2 + gender + age +  edu  + indicator2010 + PFASTYN + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10, data=toteffectsdata)
mod3b=lm(PHDLD~PFASTYN+cumbur6_zscore2  + age + gender + edu + indicator2010 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + smoking , data=toteffectsdata)

summary(mod3a)
summary(mod3b)

#CRP
mod4a=lm(crp3 ~ cumbur6_zscore2 + gender + age +  edu  + indicator2010 + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10, data=toteffectsdata)
mod4b=lm(crp3~cumbur6_zscore2  + age + gender + edu + indicator2010 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10+smoking, data=toteffectsdata)
summary(mod4a) 
summary(mod4b) 

#Glucose
mod5a=lm(glucose ~ cumbur6_zscore2 + gender + age +  edu  + indicator2010 + PFASTYN + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10, data=nomeds)
mod5b=lm(glucose~cumbur6_zscore2 + PFASTYN + age + gender + edu + indicator2010 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +smoking, data=nomeds)

summary(mod5a)
summary(mod5b)


#SBP -- note that "SBP_new" has already pre-adjusted for hypertension medication use as described in the methods
mod6a=lm(SBP_new ~ cumbur6_zscore2 + gender + age +  edu  + indicator2010 + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10, data=toteffectsdata)
mod6b=lm(SBP_new~cumbur6_zscore2  + age + gender + edu + indicator2010 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10+smoking, data=toteffectsdata)

summary(mod6a)
summary(mod6b)


#DBP -- note that "DBP_new" has already pre-adjusted for hypertension medication use as described in the methods
mod7a=lm(DBP_new ~ cumbur6_zscore2 + gender + age +  edu  + indicator2010 + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10, data=toteffectsdata)
mod7b=lm(DBP_new~cumbur6_zscore2 + age + gender + edu + indicator2010 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +smoking, data=toteffectsdata)

summary(mod7a)
summary(mod7b)


#BMI
mod8a=lm(BMI_impute ~ cumbur6_zscore2 + gender + age +  edu  + indicator2010 + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10, data=toteffectsdata)
mod8b=lm(BMI_impute ~ cumbur6_zscore2 + gender + age +  edu  + indicator2010 + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10 + smoking, data=toteffectsdata)

summary(mod8a)
summary(mod8b)


#Waist Circumference
mod9a=lm(WC_new ~ cumbur6_zscore2 + gender + age +  edu  + indicator2010 + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10, data=toteffectsdata)
mod9b=lm(WC_new~cumbur6_zscore2  + age + gender + edu + indicator2010 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + smoking, data=toteffectsdata)

summary(mod9a)
summary(mod9b)

#Triglycerides
mod10a=lm(PTGF2 ~ cumbur6_zscore2 + gender + age +  edu  + indicator2010 + PFASTYN + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10, data=toteffectsdata)
mod10b=lm(PTGF2~cumbur6_zscore2  + age + gender + edu + indicator2010 + PFASTYN + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + smoking , data=toteffectsdata)

summary(mod10a)
summary(mod10b)



