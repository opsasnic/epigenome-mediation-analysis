library(MASS)
library(hdi)
library(HDMT)
library(HIMA)


############################################BMI##############################################
#Reading in BMI phenotype, CpG mediators, and merging
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
colnames(PC.comb) = c('local_id', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')

###Adjusting for smoking 
hrs_bmi = merge(ptype2, PC.comb, by='local_id')
hrs_bmi=hrs_bmi[,-1]
hrs_bmi=hrs_bmi[,-10:-22]
hrs_bmi=hrs_bmi[,-6]

bmimediators=merge(hrs_bmi,methy, by="id")
bmimediators1=bmimediators[,-1]
bmimediators2=na.omit(bmimediators1)
mediators=bmimediators2[18:67]
pheno=bmimediators2[1:17]

x=bmimediators2$cumbur6_zscore2
y=bmimediators2$BMI_impute
x=as.numeric(x)
y=as.numeric(y)


BMI_hima2 <- hima2(BMI_impute ~ cumbur6_zscore2 + gender + age + smoking + edu + indicator2010 + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10,
                   data.pheno = pheno,
                   data.M = mediators,
                   outcome.family = "gaussian",
                   mediator.family = "gaussian",
                   penalty = "DBlasso",
                   scale = TRUE,verbose = TRUE,FDRcut = 0.05) 


###NOT adjusting for smoking 
hrs_bmi = merge(ptype2, PC.comb, by='local_id')
hrs_bmi=hrs_bmi[,-1]
hrs_bmi=hrs_bmi[,-10:-22]
hrs_bmi=hrs_bmi[,-8]
hrs_bmi=hrs_bmi[,-6]


bmimediators=merge(hrs_bmi,methy, by="id")
bmimediators1=bmimediators[,-1]
bmimediators2=na.omit(bmimediators1)
mediators=bmimediators2[17:66]
pheno=bmimediators2[1:16]

x=bmimediators2$cumbur6_zscore
y=bmimediators2$BMI_impute
x=as.numeric(x)
y=as.numeric(y)

bmi_hima2 <- hima2(BMI_impute ~ cumbur6_zscore2 + gender + age + edu + indicator2010 + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10,
                   data.pheno = pheno,
                   data.M = mediators,
                   outcome.family = "gaussian",
                   mediator.family = "gaussian",
                   penalty = "DBlasso",
                   scale = TRUE,verbose = TRUE,FDRcut = 0.05) 

path = '/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/Significant/HIMA/'
select_path = paste(path, 'BMI', '.csv', sep='')
write.csv(bmi_hima2, select_path , row.names=T)



####################################WC#########################################################333
#Reading in WC phenotype, CpG mediators, and merging
wc = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/HIMA/WCdata.csv', header=T)
mediators = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/HIMA/cpgs/WCmediators_new.csv', header=T)

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
colnames(PC.comb) = c('local_id', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')


###Adjusting for smoking
hrs_wc = merge(ptype2, PC.comb, by='local_id')
hrs_wc=hrs_wc[,-10:-22]
hrs_wc=hrs_wc[,-1]

wcmediators=merge(hrs_wc,mediators, by="id")
wcmediators2=na.omit(wcmediators)

mediators=wcmediators2[,19:64]
pheno=wcmediators2[1:18]

x=wcmediators2$cumbur6_zscore2
y=wcmediators2$WC_new
x=as.numeric(x)
y=as.numeric(y)
wc_hima2 <- hima2(WC_new ~ cumbur6_zscore2 + gender + age + smoking + edu + indicator2010 + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10,
                  data.pheno = pheno,
                  data.M = mediators,
                  outcome.family = "gaussian",
                  mediator.family = "gaussian",
                  penalty = "DBlasso",
                  scale = TRUE,verbose = TRUE,FDRcut = 0.05) 


###NOT adjusting for smoking 
hrs_wc = merge(ptype2, PC.comb, by='local_id')
hrs_wc=hrs_wc[,-10:-22]
hrs_wc=hrs_wc[,-1]
hrs_wc=hrs_wc[,-7]

wcmediators=merge(hrs_wc,mediators, by="id")
wcmediators2=na.omit(wcmediators)

mediators=wcmediators2[,18:63]
pheno=wcmediators2[1:17]

x=wcmediators2$cumbur6_zscore
y=wcmediators2$WC_new
x=as.numeric(x)
y=as.numeric(y)

WC_hima2 <- hima2(WC_new ~ cumbur6_zscore2 + gender + age +  edu + indicator2010 + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10,
                  data.pheno = pheno,
                  data.M = mediators,
                  outcome.family = "gaussian",
                  mediator.family = "gaussian",
                  penalty = "DBlasso",
                  scale = TRUE,verbose = TRUE,FDRcut = 0.05) 


####################################HDL-C#########################################################333
#Reading in HDL-C phenotype, CpG mediators, and merging
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
colnames(PC.comb) = c('local_id', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')


###Adjusting for smoking and BMI
hrs_hdl = merge(ptype2, PC.comb, by='local_id')
hrs_hdl=hrs_hdl[,-12:-24]
hrs_hdl=hrs_hdl[,-1]
hrs_hdl=hrs_hdl[,-6]


hdlmediators=merge(hrs_hdl,methy, by="id")
hdlmediators2=na.omit(hdlmediators)
mediators=hdlmediators2[,20:26]
pheno=hdlmediators2[,2:19]

x=hdlmediators2$cumbur6_zscore
y=hdlmediators2$PHDLD
x=as.numeric(x)
y=as.numeric(y)

hdl_hima2 <- hima2(PHDLD ~ cumbur6_zscore2 + gender + age + smoking + edu + PFASTYN + indicator2010 + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10,
                   data.pheno = pheno,
                   data.M = mediators,
                   outcome.family = "gaussian",
                   mediator.family = "gaussian",
                   penalty = "DBlasso",
                   scale = TRUE,verbose = TRUE,FDRcut = 0.05) 


###NOT adjusting for smoking 
hrs_hdl = merge(ptype2, PC.comb, by='local_id')
hrs_hdl=hrs_hdl[,-12:-24]
hrs_hdl=hrs_hdl[,-1]
hrs_hdl=hrs_hdl[,-6]
hrs_hdl=hrs_hdl[,-7]


hdlmediators=merge(hrs_hdl,methy, by="id")
hdlmediators2=na.omit(hdlmediators)
mediators=hdlmediators2[,19:25]
pheno=hdlmediators2[,2:18]

x=hdlmediators2$cumbur6_zscore
y=hdlmediators2$PHDLD
x=as.numeric(x)
y=as.numeric(y)

HDL_hima2 <- hima2(PHDLD ~ cumbur6_zscore2 + gender + age +  edu + PFASTYN  + indicator2010 + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10,
                   data.pheno = pheno,
                   data.M = mediators,
                   outcome.family = "gaussian",
                   mediator.family = "gaussian",
                   penalty = "DBlasso",
                   scale = TRUE,verbose = TRUE,FDRcut = 0.05) 



####################################CRP#########################################################333
#Reading in CRP phenotype, CpG mediators, and merging
crp = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/HIMA/CRPdata.csv', header=T)
mediators = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/HIMA/cpgs/CRPmediators_new.csv', header=F, skip=1)
header = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/HIMA/cpgs/CRPmediators_new.csv', header=F, nrows=1)

methy = data.frame(id=unlist(header)[-1], t(mediators[,-c(1)]), row.names = NULL)
colnames(methy)[-1] = as.vector(mediators$V1)
rm(mediators)


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


###Adjusting for smoking
hrs_crp = merge(ptype2, PC.comb, by='local_id')
hrs_crp=hrs_crp[,-10:-22]
hrs_crp=hrs_crp[,-1]

crpmediators=merge(hrs_crp,methy, by="id")
crpmediators1=crpmediators[,-1]
crpmediators2=na.omit(crpmediators1)
mediators=crpmediators2[18:29]
pheno=crpmediators2[1:17]

x=crpmediators2$cumbur6_zscore2
y=crpmediators2$crp2
x=as.numeric(x)
y=as.numeric(y)

crp_hima2 <- hima2(crp2 ~ cumbur6_zscore2 + gender + age + smoking + edu + indicator2010 + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10,
                   data.pheno = pheno,
                   data.M = mediators,
                   outcome.family = "gaussian",
                   mediator.family = "gaussian",
                   penalty = "DBlasso",
                   scale = TRUE,verbose = TRUE,FDRcut = 0.05) 

###NOT adjusting for smoking 
hrs_crp = merge(ptype2, PC.comb, by='local_id')
hrs_crp=hrs_crp[,-10:-22]
hrs_crp=hrs_crp[,-1]

crpmediators=merge(hrs_crp,methy, by="id")
crpmediators1=crpmediators[,-1]
crpmediators1=crpmediators1[,-6]
crpmediators2=na.omit(crpmediators1)
mediators=crpmediators2[17:28]
pheno=crpmediators2[1:16]

x=crpmediators2$cumbur6_zscore2
y=crpmediators2$crp2
x=as.numeric(x)
y=as.numeric(y)

CRP_hima2 <- hima2(crp2 ~ cumbur6_zscore2 + gender + age +  edu + indicator2010 + PC1+PC2+PC3+PC4+PC5+PC8+PC7+PC8+PC9+PC10,
                   data.pheno = pheno,
                   data.M = mediators,
                   outcome.family = "gaussian",
                   mediator.family = "gaussian",
                   penalty = "DBlasso",
                   scale = TRUE,verbose = TRUE,FDRcut = 0.05) 


