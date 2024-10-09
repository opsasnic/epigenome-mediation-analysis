library(data.table)

#input bash arguments
args = commandArgs(trailingOnly=TRUE)
start = as.numeric(args[1])
filenum = start/50000 + 1

#Reading in methylation data 
header = read.csv("/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Residuals2/residuals_all.csv", header=F, nrows=1)
meth = read.csv("/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Residuals2/residuals_all.csv", header=F,  nrows=50000, skip=(start+1))

methy = data.frame(id=unlist(header)[-1], t(meth[,-c(1)]), row.names = NULL)
colnames(methy)[-1] = as.vector(meth$V1)
rm(meth)

n = dim(methy)[2]-1
cpg = colnames(methy)[-1]

#Reading in phenotype file 
ptype = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/aim2analysis_final2.csv', header=T)

ptype$gender=factor(ptype$gender)
ptype$edu=factor(ptype$edu)
ptype$cumbur6_zscore2=as.numeric(ptype$cumbur6_zscore2)
ptype$indicator2010=factor(ptype$indicator2010)
ptype$age=as.numeric(ptype$age)
ptype$PFASTYN=as.factor(ptype$PFASTYN)


#Reading in cross walk data and technical covariates
xwalk = read.csv('/net/orion/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Methylation_Apr_2021/transfer2/xwalk.csv', header=T)
plateinfo = read.csv('/net/orion/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Methylation_Apr_2021/transfer04.13.2021/finalPD.csv', header=T)[,c(23, 3, 4, 9, 8)]
cells = read.csv('/net/orion/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Methylation_Apr_2021/transfer2/DNAm_cellpropest_HRSn4018.csv', header=T)
info = merge(merge(xwalk, plateinfo, by='FID'), cells, by.x='id', by.y='sample')
info$row = as.factor(substr(info$BCD_Well, 1, 1))
info$col = as.factor(substr(info$BCD_Well, 2, 3))

ptype2 = merge(ptype, info[,-c(4,5)], by='id')
ptype2$Slide = as.factor(ptype2$Slide)


#Reading in genetic ancestry PCs
PC.comb = read.csv('/net/orion/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Genotype_Jan2015_Processing/PCA/Phase1234/Top_PCs.csv')
colnames(PC.comb) = c('local_id', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')
data = merge(ptype2, PC.comb, by='local_id')
  
hrs = merge(data, methy, by='id')

#Running models in loop
snpmod1 = function(out) {
  temp = na.omit(hrs[,c(out, 'WC_new','cumbur6_zscore2', 'smoking', 'age', 'edu', 'gender', 'indicator2010', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')])
  MODEL1 = as.formula(parse(text=paste(out, ' ~ cumbur6_zscore2 + smoking + age  + gender + edu + indicator2010 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10', sep=''))[[1]])
  
  lm_result1 = summary(lm(MODEL1, data=temp))
  
  results1=rep(NA,6)
  results1[1] = out
  results1[2] = dim(temp)[1]
  results1[3:6] = lm_result1$coefficients['cumbur6_zscore2', c("Estimate", "Std. Error", "t value", "Pr(>|t|)")]
  
  return (results1) 
}


path = '/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/WC2/NewAdj/'


model1 = matrix(NA, n, 6)
for (i in 1:n){
  M = cpg[i]
  model1[i,] = snpmod1(M)
  print(i)
}

#Saving off results 
colnames(model1) = c('cpg', 'n', 'coef', 'se', 'tval', 'pval')
output = paste0(path, 'HRS.', 'aim2.', 'expmed.', filenum, '.csv')
write.csv(model1, output , row.names=F)