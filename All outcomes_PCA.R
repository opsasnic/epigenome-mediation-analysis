library(MASS)
library(hdi)
library(HDMT)
library(HIMA)
library(zoo)
library(mediation)

########################################BMI################################################################

#Reading in BMI phenotype, CpG mediators, and merging data
bmi = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/HIMA/BMIdata.csv', header=T)
mediators = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/HIMA/cpgs/BMImediators.csv', header=F, skip=1)
header = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/HIMA/cpgs/BMImediators.csv', header=F, nrows=1)

methy = data.frame(id=unlist(header)[-1], t(mediators[,-c(1)]), row.names = NULL)
colnames(methy)[-1] = as.vector(mediators$V1)
rm(mediators)

xwalk = read.csv('/net/orion/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Methylation_Apr_2021/transfer2/xwalk.csv', header=T)
plateinfo = read.csv('/net/orion/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Methylation_Apr_2021/transfer04.13.2021/finalPD.csv', header=T)[,c(23, 3, 4, 9, 8)]
cells = read.csv('/net/orion/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Methylation_Apr_2021/transfer2/DNAm_cellpropest_HRSn4018.csv', header=T)
info = merge(merge(xwalk, plateinfo, by='FID'), cells, by.x='id', by.y='sample')
info$row = as.factor(substr(info$BCD_Well, 1, 1))
info$col = as.factor(substr(info$BCD_Well, 2, 3))

ptype2 = merge(bmi, info[,-c(4,5)], by='id')
ptype2$Slide = as.factor(ptype2$Slide)
#PCs
PC.comb = read.csv('/net/orion/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Genotype_Jan2015_Processing/PCA/Phase1234/Top_PCs.csv')
colnames(PC.comb) = c('local_id', 'PC1.x', 'PC2.x', 'PC3.x', 'PC4.x', 'PC5.x', 'PC6.x', 'PC7.x', 'PC8.x', 'PC9.x', 'PC10.x')

hrs_bmi = merge(ptype2, PC.comb, by='local_id')
hrs_bmi=hrs_bmi[,-1]
hrs_bmi=hrs_bmi[,-10:-22]
hrs_bmi=hrs_bmi[,-6]

bmimediators=merge(hrs_bmi,methy, by="id")
bmimediators=na.omit(bmimediators)

bmimediators=merge(hrs_bmi,methy, by="id")
bmimediators1=bmimediators[,-1]
bmimediators2=na.omit(bmimediators1)
meds=bmimediators2[18:67]

#PC analysis 
PC_bmi = prcomp(meds, retx = T, center = T, scale. = F, rank.=10)
PC_bminew=as.data.frame(PC_bmi$x)

PC1=cbind(bmimediators2,PC_bminew)
PC1_1=PC1[,-18:-67]

summary(PC_bmi)

####Mediation Adjusting for smoking#####

#Mediate PC1
out.PC1=lm(BMI_impute~cumbur6_zscore2 + age + gender + edu + smoking  + PC1 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC1=lm(PC1~cumbur6_zscore2 + age + gender + edu + smoking  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC1_1<- mediate(med.PC1, out.PC1, treat = "cumbur6_zscore2", mediator = "PC1", robustSE = TRUE, sims = 10000)
summary(med.out_PC1_1)

#Mediate PC2
out.PC2=lm(BMI_impute~cumbur6_zscore2 + age + gender + edu + smoking + PC2 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC2=lm(PC2~cumbur6_zscore2 + age + gender + edu + smoking  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC2<- mediate(med.PC2, out.PC2, treat = "cumbur6_zscore2", mediator = "PC2", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC2)

#Mediate PC3
out.PC3=lm(BMI_impute~cumbur6_zscore2 + age + gender + edu + smoking + PC3 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC3=lm(PC3~cumbur6_zscore2 + age + gender + edu + smoking  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC3<- mediate(med.PC3, out.PC3, treat = "cumbur6_zscore2", mediator = "PC3", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC3)

#Mediate PC4
out.PC4=lm(BMI_impute~cumbur6_zscore2 + age + gender + edu + smoking + PC4 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC4=lm(PC4~cumbur6_zscore2 + age + gender + edu + smoking  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC4<- mediate(med.PC4, out.PC4, treat = "cumbur6_zscore2", mediator = "PC4", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC4)

#Mediate PC5
out.PC5=lm(BMI_impute~cumbur6_zscore2 + age + gender + edu + smoking  + PC5 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC5=lm(PC5~cumbur6_zscore2 + age + gender + edu + smoking  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC5<- mediate(med.PC5, out.PC5, treat = "cumbur6_zscore2", mediator = "PC5", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC5)

#Mediate PC6
out.PC6=lm(BMI_impute~cumbur6_zscore2 + age + gender + edu + smoking  + PC6 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC6=lm(PC6~cumbur6_zscore2 + age + gender + edu + smoking  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC6<- mediate(med.PC6, out.PC6, treat = "cumbur6_zscore2", mediator = "PC6", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC6)

#Mediate PC7
out.PC7=lm(BMI_impute~cumbur6_zscore2 + age + gender + edu + smoking + PC7 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC7=lm(PC7~cumbur6_zscore2 + age + gender + edu + smoking   + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC7<- mediate(med.PC7, out.PC7, treat = "cumbur6_zscore2", mediator = "PC7", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC7)

#Mediate PC8
out.PC8=lm(BMI_impute~cumbur6_zscore2 + age + gender + edu + smoking  + PC8 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC8=lm(PC8~cumbur6_zscore2 + age + gender + edu + smoking  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC8<- mediate(med.PC8, out.PC8, treat = "cumbur6_zscore2", mediator = "PC8", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC8)

#Mediate PC9
out.PC9=lm(BMI_impute~cumbur6_zscore2 + age + gender + edu + smoking  + PC9 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC9=lm(PC9~cumbur6_zscore2 + age + gender + edu + smoking  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC9<- mediate(med.PC9, out.PC9, treat = "cumbur6_zscore2", mediator = "PC9", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC9)

#Mediate PC10
out.PC10=lm(BMI_impute~cumbur6_zscore2 + age + gender + edu + smoking + PC10 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC10=lm(PC10~cumbur6_zscore2 + age + gender + edu + smoking  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC10<- mediate(med.PC10, out.PC10, treat = "cumbur6_zscore2", mediator = "PC10", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC10)

####Mediation without Adjusting for smoking#####

#Mediate PC1
out.PC1=lm(BMI_impute~cumbur6_zscore2 + age + gender + edu  + PC1 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC1=lm(PC1~cumbur6_zscore2 + age + gender + edu  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC1_1<- mediate(med.PC1, out.PC1, treat = "cumbur6_zscore2", mediator = "PC1", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
print(summary(med.out_PC1_1))


#Mediate PC2
out.PC2=lm(BMI_impute~cumbur6_zscore2 + age + gender + edu  + PC2 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC2=lm(PC2~cumbur6_zscore2 + age + gender + edu  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC2<- mediate(med.PC2, out.PC2, treat = "cumbur6_zscore2", mediator = "PC2", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC2)

#Mediate PC3
out.PC3=lm(BMI_impute~cumbur6_zscore2 + age + gender + edu  + PC3 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC3=lm(PC3~cumbur6_zscore2 + age + gender + edu  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC3<- mediate(med.PC3, out.PC3, treat = "cumbur6_zscore2", mediator = "PC3", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC3)

#Mediate PC4
out.PC4=lm(BMI_impute~cumbur6_zscore2 + age + gender + edu + PC4 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC4=lm(PC4~cumbur6_zscore2 + age + gender + edu   + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC4<- mediate(med.PC4, out.PC4, treat = "cumbur6_zscore2", mediator = "PC4", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC4)

#Mediate PC5
out.PC5=lm(BMI_impute~cumbur6_zscore2 + age + gender + edu + PC5 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC5=lm(PC5~cumbur6_zscore2 + age + gender + edu   + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC5<- mediate(med.PC5, out.PC5, treat = "cumbur6_zscore2", mediator = "PC5", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC5)

#Mediate PC6
out.PC6=lm(BMI_impute~cumbur6_zscore2 + age + gender + edu + PC6 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC6=lm(PC6~cumbur6_zscore2 + age + gender + edu  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC6<- mediate(med.PC6, out.PC6, treat = "cumbur6_zscore2", mediator = "PC6", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC6)

#Mediate PC7
out.PC7=lm(BMI_impute~cumbur6_zscore2 + age + gender + edu + PC7 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC7=lm(PC7~cumbur6_zscore2 + age + gender + edu  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC7<- mediate(med.PC7, out.PC7, treat = "cumbur6_zscore2", mediator = "PC7", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC7)

#Mediate PC8
out.PC8=lm(BMI_impute~cumbur6_zscore2 + age + gender + edu  + PC8 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC8=lm(PC8~cumbur6_zscore2 + age + gender + edu  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC8<- mediate(med.PC8, out.PC8, treat = "cumbur6_zscore2", mediator = "PC8", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC8)

#Mediate PC9
out.PC9=lm(BMI_impute~cumbur6_zscore2 + age + gender + edu  + PC9 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC9=lm(PC9~cumbur6_zscore2 + age + gender + edu  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC9<- mediate(med.PC9, out.PC9, treat = "cumbur6_zscore2", mediator = "PC9", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC9)

#Mediate PC10
out.PC10=lm(BMI_impute~cumbur6_zscore2 + age + gender + edu + PC10 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC10=lm(PC10~cumbur6_zscore2 + age + gender + edu  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC10<- mediate(med.PC10, out.PC10, treat = "cumbur6_zscore2", mediator = "PC10", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC10)

#Calculating total Variability of Stress on BMI
with=lm(BMI_impute ~ cumbur6_zscore2 + age + gender + edu + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=bmimediators2 )
without=lm(BMI_impute ~ age + gender + edu + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=bmimediators2 )

((summary(with)$r.squared-summary(without)$r.squared)/summary(with)$r.squared)*100
((summary(with)$adj.r.squared-summary(without)$adj.r.squared)/summary(with)$adj.r.squared)*100


print.summary.mediate <- mediation:::print.summary.mediate
fix(print.summary.mediate)
printCoefmat(smat, digits = 3)


########################################WC################################################################

#Reading in WC phenotype, CpG mediators, and merging data
wc = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/HIMA/WCdata.csv', header=T)
methy = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/HIMA/cpgs/WCmediators_new.csv', header=T)

xwalk = read.csv('/net/orion/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Methylation_Apr_2021/transfer2/xwalk.csv', header=T)
plateinfo = read.csv('/net/orion/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Methylation_Apr_2021/transfer04.13.2021/finalPD.csv', header=T)[,c(23, 3, 4, 9, 8)]
cells = read.csv('/net/orion/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Methylation_Apr_2021/transfer2/DNAm_cellpropest_HRSn4018.csv', header=T)
info = merge(merge(xwalk, plateinfo, by='FID'), cells, by.x='id', by.y='sample')
info$row = as.factor(substr(info$BCD_Well, 1, 1))
info$col = as.factor(substr(info$BCD_Well, 2, 3))

ptype2 = merge(wc, info[,-c(4,5)], by='id')
ptype2$Slide = as.factor(ptype2$Slide)
#PCs
PC.comb = read.csv('/net/orion/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Genotype_Jan2015_Processing/PCA/Phase1234/Top_PCs.csv')
colnames(PC.comb) = c('local_id', 'PC1.x', 'PC2.x', 'PC3.x', 'PC4.x', 'PC5.x', 'PC6.x', 'PC7.x', 'PC8.x', 'PC9.x', 'PC10.x')

hrs_wc = merge(ptype2, PC.comb, by='local_id')
hrs_wc=hrs_wc[,-10:-22]
hrs_wc=hrs_wc[,-1]

wcmediators=merge(hrs_wc,methy, by="id")
wcmediators1=wcmediators[,-1]
wcmediators2=na.omit(wcmediators1)
meds=wcmediators2[18:63]

#PC analysis
PC_wc = prcomp(meds, retx = T, center = T, scale. = F, rank.=10)
PC_wcnew=as.data.frame(PC_wc$x)

PC1=cbind(wcmediators2,PC_wcnew)
PC1_1=PC1[,-18:-63]

summary(PC_wc)

####Mediation Adjusting for smoking#####

#Mediate PC1
out.PC1=lm(WC_new~cumbur6_zscore2 + age + gender + edu + smoking + PC1 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC1=lm(PC1~cumbur6_zscore2 + age + gender + edu + smoking  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC1_1<- mediate(med.PC1, out.PC1, treat = "cumbur6_zscore2", mediator = "PC1", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC1_1)

#Mediate PC2
out.PC2=lm(WC_new~cumbur6_zscore2 + age + gender + edu + smoking + PC2 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC2=lm(PC2~cumbur6_zscore2 + age + gender + edu + smoking  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC2<- mediate(med.PC2, out.PC2, treat = "cumbur6_zscore2", mediator = "PC2", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC2)

#Mediate PC3
out.PC3=lm(WC_new~cumbur6_zscore2 + age + gender + edu + smoking + PC3 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC3=lm(PC3~cumbur6_zscore2 + age + gender + edu + smoking  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC3<- mediate(med.PC3, out.PC3, treat = "cumbur6_zscore2", mediator = "PC3", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC3)

#Mediate PC4
out.PC4=lm(WC_new~cumbur6_zscore2 + age + gender + edu + smoking + PC4 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC4=lm(PC4~cumbur6_zscore2 + age + gender + edu + smoking  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC4<- mediate(med.PC4, out.PC4, treat = "cumbur6_zscore2", mediator = "PC4", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC4)

#Mediate PC5
out.PC5=lm(WC_new~cumbur6_zscore2 + age + gender + edu + smoking + PC5 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC5=lm(PC5~cumbur6_zscore2 + age + gender + edu + smoking  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC5<- mediate(med.PC5, out.PC5, treat = "cumbur6_zscore2", mediator = "PC5", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC5)

#Mediate PC6
out.PC6=lm(WC_new~cumbur6_zscore2 + age + gender + edu + smoking + PC6 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC6=lm(PC6~cumbur6_zscore2 + age + gender + edu + smoking  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC6<- mediate(med.PC6, out.PC6, treat = "cumbur6_zscore2", mediator = "PC6", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC6)

#Mediate PC7
out.PC7=lm(WC_new~cumbur6_zscore2 + age + gender + edu + smoking + PC7 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC7=lm(PC7~cumbur6_zscore2 + age + gender + edu + smoking  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC7<- mediate(med.PC7, out.PC7, treat = "cumbur6_zscore2", mediator = "PC7", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC7)

#Mediate PC8
out.PC8=lm(WC_new~cumbur6_zscore2 + age + gender + edu + smoking + PC8 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC8=lm(PC8~cumbur6_zscore2 + age + gender + edu + smoking  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC8<- mediate(med.PC8, out.PC8, treat = "cumbur6_zscore2", mediator = "PC8", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC8)

#Mediate PC9
out.PC9=lm(WC_new~cumbur6_zscore2 + age + gender + edu + smoking + PC9 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC9=lm(PC9~cumbur6_zscore2 + age + gender + edu + smoking  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC9<- mediate(med.PC9, out.PC9, treat = "cumbur6_zscore2", mediator = "PC9", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC9)

#Mediate PC10
out.PC10=lm(WC_new~cumbur6_zscore2 + age + gender + edu + smoking + PC10 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC10=lm(PC10~cumbur6_zscore2 + age + gender + edu + smoking  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC10<- mediate(med.PC10, out.PC10, treat = "cumbur6_zscore2", mediator = "PC10", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC10)



####Mediation without Adjusting for smoking#####


#Mediate PC1
out.PC1=lm(WC_new~cumbur6_zscore2 + age + gender + edu  + PC1 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC1=lm(PC1~cumbur6_zscore2 + age + gender + edu  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC1_1<- mediate(med.PC1, out.PC1, treat = "cumbur6_zscore2", mediator = "PC1", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC1_1)

#Mediate PC2
out.PC2=lm(WC_new~cumbur6_zscore2 + age + gender + edu  + PC2 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC2=lm(PC2~cumbur6_zscore2 + age + gender + edu  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC2<- mediate(med.PC2, out.PC2, treat = "cumbur6_zscore2", mediator = "PC2", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC2)

#Mediate PC3
out.PC3=lm(WC_new~cumbur6_zscore2 + age + gender + edu  + PC3 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC3=lm(PC3~cumbur6_zscore2 + age + gender + edu  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC3<- mediate(med.PC3, out.PC3, treat = "cumbur6_zscore2", mediator = "PC3", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC3)

#Mediate PC4
out.PC4=lm(WC_new~cumbur6_zscore2 + age + gender + edu + PC4 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC4=lm(PC4~cumbur6_zscore2 + age + gender + edu   + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC4<- mediate(med.PC4, out.PC4, treat = "cumbur6_zscore2", mediator = "PC4", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC4)

#Mediate PC5
out.PC5=lm(WC_new~cumbur6_zscore2 + age + gender + edu + PC5 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC5=lm(PC5~cumbur6_zscore2 + age + gender + edu   + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC5<- mediate(med.PC5, out.PC5, treat = "cumbur6_zscore2", mediator = "PC5", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC5)

#Mediate PC6
out.PC6=lm(WC_new~cumbur6_zscore2 + age + gender + edu + PC6 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC6=lm(PC6~cumbur6_zscore2 + age + gender + edu  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC6<- mediate(med.PC6, out.PC6, treat = "cumbur6_zscore2", mediator = "PC6", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC6)

#Mediate PC7
out.PC7=lm(WC_new~cumbur6_zscore2 + age + gender + edu + PC7 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC7=lm(PC7~cumbur6_zscore2 + age + gender + edu  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC7<- mediate(med.PC7, out.PC7, treat = "cumbur6_zscore2", mediator = "PC7", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC7)

#Mediate PC8
out.PC8=lm(WC_new~cumbur6_zscore2 + age + gender + edu  + PC8 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC8=lm(PC8~cumbur6_zscore2 + age + gender + edu  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC8<- mediate(med.PC8, out.PC8, treat = "cumbur6_zscore2", mediator = "PC8", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC8)

#Mediate PC9
out.PC9=lm(WC_new~cumbur6_zscore2 + age + gender + edu  + PC9 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC9=lm(PC9~cumbur6_zscore2 + age + gender + edu  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC9<- mediate(med.PC9, out.PC9, treat = "cumbur6_zscore2", mediator = "PC9", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC9)

#Mediate PC10
out.PC10=lm(WC_new~cumbur6_zscore2 + age + gender + edu + PC10 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC10=lm(PC10~cumbur6_zscore2 + age + gender + edu  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC10<- mediate(med.PC10, out.PC10, treat = "cumbur6_zscore2", mediator = "PC10", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC10)


#Calculating total Variability of Stress on WC
with=lm(WC_new ~  cumbur6_zscore2 + age + gender + edu + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=wcmediators2 )
without=lm(WC_new ~ age + gender + edu + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=wcmediators2 )

((summary(with)$r.squared-summary(without)$r.squared)/summary(with)$r.squared)*100
((summary(with)$adj.r.squared-summary(without)$adj.r.squared)/summary(with)$adj.r.squared)*100

########################################HDL-C################################################################

#Reading in HDL-C phenotype, CpG mediators, and merging data
hdl = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/HIMA/HDLdata.csv', header=T)
mediators = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/HIMA/cpgs/HDLmediators_new.csv', header=F, skip=1)
header = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/HIMA/cpgs/HDLmediators_new.csv', header=F, nrows=1)

methy = data.frame(id=unlist(header)[-1], t(mediators[,-c(1)]), row.names = NULL)
colnames(methy)[-1] = as.vector(mediators$V1)
rm(mediators)


xwalk = read.csv('/net/orion/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Methylation_Apr_2021/transfer2/xwalk.csv', header=T)
plateinfo = read.csv('/net/orion/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Methylation_Apr_2021/transfer04.13.2021/finalPD.csv', header=T)[,c(23, 3, 4, 9, 8)]
cells = read.csv('/net/orion/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Methylation_Apr_2021/transfer2/DNAm_cellpropest_HRSn4018.csv', header=T)
info = merge(merge(xwalk, plateinfo, by='FID'), cells, by.x='id', by.y='sample')
info$row = as.factor(substr(info$BCD_Well, 1, 1))
info$col = as.factor(substr(info$BCD_Well, 2, 3))

ptype2 = merge(hdl, info[,-c(4,5)], by='id')
ptype2$Slide = as.factor(ptype2$Slide)
#PCs
PC.comb = read.csv('/net/orion/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Genotype_Jan2015_Processing/PCA/Phase1234/Top_PCs.csv')
colnames(PC.comb) = c('local_id', 'PC1.x', 'PC2.x', 'PC3.x', 'PC4.x', 'PC5.x', 'PC6.x', 'PC7.x', 'PC8.x', 'PC9.x', 'PC10.x')

hrs_hdl = merge(ptype2, PC.comb, by='local_id')
hrs_hdl=hrs_hdl[,-12:-24]
hrs_hdl=hrs_hdl[,-1]
hrs_hdl=hrs_hdl[,-6]


hdlpmediators=merge(hrs_hdl,methy, by="id")
hdlpmediators1=hdlpmediators[,-1]
hdlpmediators2=na.omit(hdlpmediators1)
meds=hdlpmediators2[19:25]

#PC analysis 
PC_hdl = prcomp(meds, retx = T, center = T, scale. = F, rank.=5)
PC_hdlnew=as.data.frame(PC_hdl$x)

PC1=cbind(hdlpmediators2,PC_hdlnew)
PC1_1=PC1[,-19:-25]

summary(PC_hdl)


####Mediation Adjusting for smoking#####

#Mediate PC1
out.PC1=lm(PHDLD~cumbur6_zscore2 + age + gender + edu + smoking  + PFASTYN + PC1 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC1=lm(PC1~cumbur6_zscore2 + age + gender + edu + smoking   + PFASTYN + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC1_1<- mediate(med.PC1, out.PC1, treat = "cumbur6_zscore2", mediator = "PC1", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC1_1)

#Mediate PC2
out.PC2=lm(PHDLD~cumbur6_zscore2 + age + gender + edu + smoking  + PFASTYN + PC2 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC2=lm(PC2~cumbur6_zscore2 + age + gender + edu + smoking   + PFASTYN + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC2<- mediate(med.PC2, out.PC2, treat = "cumbur6_zscore2", mediator = "PC2", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC2)

#Mediate PC3
out.PC3=lm(PHDLD~cumbur6_zscore2 + age + gender + edu + smoking  + PFASTYN + PC3 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC3=lm(PC3~cumbur6_zscore2 + age + gender + edu + smoking   + PFASTYN + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC3<- mediate(med.PC3, out.PC3, treat = "cumbur6_zscore2", mediator = "PC3", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC3)

#Mediate PC4
out.PC4=lm(PHDLD~cumbur6_zscore2 + age + gender + edu + smoking  + PFASTYN+ PC4 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC4=lm(PC4~cumbur6_zscore2 + age + gender + edu + smoking  + PFASTYN  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC4<- mediate(med.PC4, out.PC4, treat = "cumbur6_zscore2", mediator = "PC4", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC4)

#Mediate PC5
out.PC5=lm(PHDLD~cumbur6_zscore2 + age + gender + edu + smoking  + PFASTYN + PC5 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC5=lm(PC5~cumbur6_zscore2 + age + gender + edu + smoking  + PFASTYN + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC5<- mediate(med.PC5, out.PC5, treat = "cumbur6_zscore2", mediator = "PC5", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC5)

####Mediation without Adjusting for smoking#####
#Mediate PC1
out.PC1=lm(PHDLD~cumbur6_zscore2 + age + gender + edu  + PC1 + indicator2010 + PFASTYN  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC1=lm(PC1~cumbur6_zscore2 + age + gender + edu  + indicator2010  + PFASTYN + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC1_1<- mediate(med.PC1, out.PC1, treat = "cumbur6_zscore2", mediator = "PC1", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC1_1)

#Mediate PC2
out.PC2=lm(PHDLD~cumbur6_zscore2 + age + gender + edu  + PC2 + indicator2010 + PFASTYN  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC2=lm(PC2~cumbur6_zscore2 + age + gender + edu  + indicator2010 + PFASTYN + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC2<- mediate(med.PC2, out.PC2, treat = "cumbur6_zscore2", mediator = "PC2", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC2)

#Mediate PC3
out.PC3=lm(PHDLD~cumbur6_zscore2 + age + gender + edu  + PC3 + indicator2010 + PFASTYN  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC3=lm(PC3~cumbur6_zscore2 + age + gender + edu  + indicator2010 + PFASTYN + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC3<- mediate(med.PC3, out.PC3, treat = "cumbur6_zscore2", mediator = "PC3", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC3)

#Mediate PC4
out.PC4=lm(PHDLD~cumbur6_zscore2 + age + gender + edu + PC4 + indicator2010 + PFASTYN  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC4=lm(PC4~cumbur6_zscore2 + age + gender + edu   + indicator2010 + PFASTYN + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC4<- mediate(med.PC4, out.PC4, treat = "cumbur6_zscore2", mediator = "PC4", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC4)

#Mediate PC5
out.PC5=lm(PHDLD~cumbur6_zscore2 + age + gender + edu + PC5 + indicator2010 + PFASTYN  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC5=lm(PC5~cumbur6_zscore2 + age + gender + edu   + indicator2010 + PFASTYN + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC5<- mediate(med.PC5, out.PC5, treat = "cumbur6_zscore2", mediator = "PC5", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC5)


#Calculating total Variability of Stress on HDL-C
with=lm(PHDLD ~  PFASTYN+ cumbur6_zscore2 + age + gender + edu + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=hdlpmediators )
without=lm(PHDLD ~ PFASTYN+ age + gender + edu + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=hdlpmediators )

((summary(with)$r.squared-summary(without)$r.squared)/summary(with)$r.squared)*100
((summary(with)$adj.r.squared-summary(without)$adj.r.squared)/summary(with)$adj.r.squared)*100


########################################CRP################################################################
#Reading in CRP phenotype, CpG mediators, and merging data
crp = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/HIMA/CRPdata.csv', header=T)
mediators = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/HIMA/cpgs/CRPmediators_new.csv', header=F, skip=1)
header = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/HIMA/cpgs/CRPmediators_new.csv', header=F, nrows=1)

methy = data.frame(id=unlist(header)[-1], t(mediators[,-c(1)]), row.names = NULL)
colnames(methy)[-1] = as.vector(mediators$V1)
rm(mediators)

#log transforming CRP measure
crp$crp2=log(crp$PCRP)
crp=crp[,-8]

xwalk = read.csv('/net/orion/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Methylation_Apr_2021/transfer2/xwalk.csv', header=T)
plateinfo = read.csv('/net/orion/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Methylation_Apr_2021/transfer04.13.2021/finalPD.csv', header=T)[,c(23, 3, 4, 9, 8)]
cells = read.csv('/net/orion/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Methylation_Apr_2021/transfer2/DNAm_cellpropest_HRSn4018.csv', header=T)
info = merge(merge(xwalk, plateinfo, by='FID'), cells, by.x='id', by.y='sample')
info$row = as.factor(substr(info$BCD_Well, 1, 1))
info$col = as.factor(substr(info$BCD_Well, 2, 3))

ptype2 = merge(crp, info[,-c(4,5)], by='id')
ptype2$Slide = as.factor(ptype2$Slide)
#PCs
PC.comb = read.csv('/net/orion/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Genotype_Jan2015_Processing/PCA/Phase1234/Top_PCs.csv')
colnames(PC.comb) = c('local_id', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')


hrs_crp = merge(ptype2, PC.comb, by='local_id')
hrs_crp=hrs_crp[,-10:-22]
hrs_crp=hrs_crp[,-1]

crpmediators=merge(hrs_crp,methy, by="id")
crpmediators1=crpmediators[,-1]
crpmediators2=na.omit(crpmediators1)
mediators=crpmediators2[18:29]

#PC analysis 
PC_crp = prcomp(mediators, retx = T, center = T, scale. = F, rank.=5)
PC_crpnew=as.data.frame(PC_crp$x)

PC1=cbind(crpmediators2,PCPC_crpnew_hdlnew)
PC1_1=PC1[,-18:-29]

summary(PC_crp)

####Mediation Adjusting for smoking#####

#Mediate PC1
#Mediate PC1
out.PC1=lm(crp2~cumbur6_zscore2 + age + gender + edu +smoking  + PC1 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC1=lm(PC1~cumbur6_zscore2 + age + gender + edu +smoking + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC1_1<- mediate(med.PC1, out.PC1, treat = "cumbur6_zscore2", mediator = "PC1", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC1_1)

#Mediate PC2
out.PC2=lm(crp2~cumbur6_zscore2 + age + gender + edu +smoking + PC2 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC2=lm(PC2~cumbur6_zscore2 + age + gender + edu +smoking + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC2<- mediate(med.PC2, out.PC2, treat = "cumbur6_zscore2", mediator = "PC2", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC2)

#Mediate PC3
out.PC3=lm(crp2~cumbur6_zscore2 + age + gender + edu +smoking + PC3 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC3=lm(PC3~cumbur6_zscore2 + age + gender + edu +smoking  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC3<- mediate(med.PC3, out.PC3, treat = "cumbur6_zscore2", mediator = "PC3", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC3)

#Mediate PC4
out.PC4=lm(crp2~cumbur6_zscore2 + age + gender + edu +smoking + PC4 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC4=lm(PC4~cumbur6_zscore2 + age + gender + edu +smoking  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC4<- mediate(med.PC4, out.PC4, treat = "cumbur6_zscore2", mediator = "PC4", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC4)

#Mediate PC5
out.PC5=lm(crp2~cumbur6_zscore2 + age + gender + edu +smoking + PC5 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC5=lm(PC5~cumbur6_zscore2 + age + gender + edu +smoking  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC5<- mediate(med.PC5, out.PC5, treat = "cumbur6_zscore2", mediator = "PC5", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC5)

#Mediate PC6
out.PC6=lm(crp2~cumbur6_zscore2 + age + gender + edu +smoking + PC6 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC6=lm(PC6~cumbur6_zscore2 + age + gender + edu +smoking + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC6<- mediate(med.PC6, out.PC6, treat = "cumbur6_zscore2", mediator = "PC6", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC6)

#Mediate PC7
out.PC7=lm(crp2~cumbur6_zscore2 + age + gender + edu  +smoking+ PC7 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC7=lm(PC7~cumbur6_zscore2 + age + gender + edu +smoking + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC7<- mediate(med.PC7, out.PC7, treat = "cumbur6_zscore2", mediator = "PC7", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC7)

#Mediate PC8
out.PC8=lm(crp2~cumbur6_zscore2 + age + gender + edu +smoking + PC8 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC8=lm(PC8~cumbur6_zscore2 + age + gender + edu +smoking + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC8<- mediate(med.PC8, out.PC8, treat = "cumbur6_zscore2", mediator = "PC8", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC8)

#Mediate PC9
out.PC9=lm(crp2~cumbur6_zscore2 + age + gender + edu +smoking  + PC9 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC9=lm(PC9~cumbur6_zscore2 + age + gender + edu +smoking + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC9<- mediate(med.PC9, out.PC9, treat = "cumbur6_zscore2", mediator = "PC9", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC9)

#Mediate PC10
out.PC10=lm(crp2~cumbur6_zscore2 + age + gender + edu +smoking + PC10 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC10=lm(PC10~cumbur6_zscore2 + age + gender + edu +smoking  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC10<- mediate(med.PC10, out.PC10, treat = "cumbur6_zscore2", mediator = "PC10", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC10)





####Mediation without Adjusting for smoking#####

#Mediate PC1
out.PC1=lm(crp2~cumbur6_zscore2 + age + gender + edu  + PC1 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC1=lm(PC1~cumbur6_zscore2 + age + gender + edu  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC1_1<- mediate(med.PC1, out.PC1, treat = "cumbur6_zscore2", mediator = "PC1", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC1_1)

#Mediate PC2
out.PC2=lm(crp2~cumbur6_zscore2 + age + gender + edu  + PC2 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC2=lm(PC2~cumbur6_zscore2 + age + gender + edu  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC2<- mediate(med.PC2, out.PC2, treat = "cumbur6_zscore2", mediator = "PC2", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC2)

#Mediate PC3
out.PC3=lm(crp2~cumbur6_zscore2 + age + gender + edu  + PC3 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC3=lm(PC3~cumbur6_zscore2 + age + gender + edu  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC3<- mediate(med.PC3, out.PC3, treat = "cumbur6_zscore2", mediator = "PC3", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC3)

#Mediate PC4
out.PC4=lm(crp2~cumbur6_zscore2 + age + gender + edu + PC4 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC4=lm(PC4~cumbur6_zscore2 + age + gender + edu   + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC4<- mediate(med.PC4, out.PC4, treat = "cumbur6_zscore2", mediator = "PC4", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC4)

#Mediate PC5
out.PC5=lm(crp2~cumbur6_zscore2 + age + gender + edu + PC5 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC5=lm(PC5~cumbur6_zscore2 + age + gender + edu   + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC5<- mediate(med.PC5, out.PC5, treat = "cumbur6_zscore2", mediator = "PC5", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC5)

#Mediate PC6
out.PC6=lm(crp2~cumbur6_zscore2 + age + gender + edu + PC6 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC6=lm(PC6~cumbur6_zscore2 + age + gender + edu  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC6<- mediate(med.PC6, out.PC6, treat = "cumbur6_zscore2", mediator = "PC6", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC6)

#Mediate PC7
out.PC7=lm(crp2~cumbur6_zscore2 + age + gender + edu + PC7 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC7=lm(PC7~cumbur6_zscore2 + age + gender + edu  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC7<- mediate(med.PC7, out.PC7, treat = "cumbur6_zscore2", mediator = "PC7", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC7)

#Mediate PC8
out.PC8=lm(crp2~cumbur6_zscore2 + age + gender + edu  + PC8 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC8=lm(PC8~cumbur6_zscore2 + age + gender + edu  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC8<- mediate(med.PC8, out.PC8, treat = "cumbur6_zscore2", mediator = "PC8", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC8)

#Mediate PC9
out.PC9=lm(crp2~cumbur6_zscore2 + age + gender + edu  + PC9 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC9=lm(PC9~cumbur6_zscore2 + age + gender + edu  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC9<- mediate(med.PC9, out.PC9, treat = "cumbur6_zscore2", mediator = "PC9", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC9)

#Mediate PC10
out.PC10=lm(crp2~cumbur6_zscore2 + age + gender + edu + PC10 + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.PC10=lm(PC10~cumbur6_zscore2 + age + gender + edu  + indicator2010 + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=PC1_1)
med.out_PC10<- mediate(med.PC10, out.PC10, treat = "cumbur6_zscore2", mediator = "PC10", robustSE = TRUE, sims = 10000, boot.ci.type = "perc")
summary(med.out_PC10)


#Calculating total Variability of Stress on CRP
with=lm(crp2 ~ cumbur6_zscore2 + age + gender + edu + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=crpmediators )
without=lm(crp2 ~ age + gender + edu + indicator2010  + PC1.x + PC2.x + PC3.x + PC4.x + PC5.x + PC6.x + PC7.x + PC8.x + PC9.x + PC10.x, data=crpmediators )

((summary(with)$r.squared-summary(without)$r.squared)/summary(with)$r.squared)*100
((summary(with)$adj.r.squared-summary(without)$adj.r.squared)/summary(with)$adj.r.squared)*100






