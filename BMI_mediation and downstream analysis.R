library(medScan)
library(HDMT)

###########Reading in exposure-mediator model####################################
all1 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/BMI2/New/HRS.aim2.expmed.1.csv', header=T)
all2 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/BMI2/New/HRS.aim2.expmed.2.csv', header=T)
all3 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/BMI2/New/HRS.aim2.expmed.3.csv', header=T)
all4 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/BMI2/New/HRS.aim2.expmed.4.csv', header=T)
all5 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/BMI2/New/HRS.aim2.expmed.5.csv', header=T)
all6 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/BMI2/New/HRS.aim2.expmed.6.csv', header=T)
all7 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/BMI2/New/HRS.aim2.expmed.7.csv', header=T)
all8 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/BMI2/New/HRS.aim2.expmed.8.csv', header=T)
all9 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/BMI2/New/HRS.aim2.expmed.9.csv', header=T)
all10 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/BMI2/New/HRS.aim2.expmed.10.csv', header=T)
all11 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/BMI2/New/HRS.aim2.expmed.11.csv', header=T)
all12 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/BMI2/New/HRS.aim2.expmed.12.csv', header=T)
all13 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/BMI2/New/HRS.aim2.expmed.13.csv', header=T)
all14 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/BMI2/New/HRS.aim2.expmed.14.csv', header=T)
all15 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/BMI2/New/HRS.aim2.expmed.15.csv', header=T)
all16 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/BMI2/New/HRS.aim2.expmed.16.csv', header=T)

ewas_stress3_fixed=rbind(all1, all2, all3, all4, all5, all6, all7, all8, all9, all10, all11, all12, all13, all14, all15, all16)
z.alpha=ewas_stress3_fixed[,5]
p.alpha = ewas_stress3_fixed$pval
names(ewas_stress3_fixed)[names(ewas_stress3_fixed) == "coef"] <- "alpha_coef"

###########Reading in mediator-outcome model####################################
BMI1 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/BMI/New/HRS.aim2.BMI.1.csv', header=T)
BMI2 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/BMI/New/HRS.aim2.BMI.2.csv', header=T)
BMI3 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/BMI/New/HRS.aim2.BMI.3.csv', header=T)
BMI4 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/BMI/New/HRS.aim2.BMI.4.csv', header=T)
BMI5 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/BMI/New/HRS.aim2.BMI.5.csv', header=T)
BMI6 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/BMI/New/HRS.aim2.BMI.6.csv', header=T)
BMI7 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/BMI/New/HRS.aim2.BMI.7.csv', header=T)
BMI8 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/BMI/New/HRS.aim2.BMI.8.csv', header=T)
BMI9 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/BMI/New/HRS.aim2.BMI.9.csv', header=T)
BMI10 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/BMI/New/HRS.aim2.BMI.10.csv', header=T)
BMI11 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/BMI/New/HRS.aim2.BMI.11.csv', header=T)
BMI12 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/BMI/New/HRS.aim2.BMI.12.csv', header=T)
BMI13 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/BMI/New/HRS.aim2.BMI.13.csv', header=T)
BMI14 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/BMI/New/HRS.aim2.BMI.14.csv', header=T)
BMI15 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/BMI/New/HRS.aim2.BMI.15.csv', header=T)
BMI16 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/BMI/New/HRS.aim2.BMI.16.csv', header=T)

ewas_bmi_fixed=rbind(BMI1, BMI2, BMI3, BMI4, BMI5, BMI6, BMI7, BMI8, BMI9, BMI10,BMI11, BMI12, BMI13, BMI14, BMI15, BMI16)
z.beta=ewas_bmi_fixed[,5]
p.beta = ewas_bmi_fixed$pval
names(ewas_bmi_fixed)[names(ewas_bmi_fixed) == "coef"] <- "beta_coef"


###Doing the high dimensional mediation using medScan######
obj = medScan(z.alpha = z.alpha, z.beta = z.beta, method = "HDMT")
objpvalues=obj$pvalues
pvalues_bmi=as.data.frame(objpvalues)
pvalues_bmi$q=p.adjust(pvalues_bmi$objpvalues, method='fdr')
median(qchisq(pvalues_bmi$objpvalues, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)

# For p-values, calculate chi-squared statistic
chisq <- qchisq(1-obj$pvalues,1)
median(chisq)/qchisq(0.5,1)

######Merging pvalue file with CPG file################
fulldata_bmi=cbind(ewas_stress3_fixed,ewas_bmi_fixed,pvalues_bmi)
fulldata_bmi=fulldata_bmi[,-7]
fulldata_bmi=fulldata_bmi[order(fulldata_bmi$q),]

#Getting significant cpg sites from mediation analysis 
significant_bmi=subset(fulldata_bmi, q<0.05)
sigcpg_bmi=significant_bmi[,1]

#Merging with annotation file
annotation=read.csv('/net/orion/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Methylation_Apr_2021/QC_processing_Smith/infinium-methylationepic-v-1-0-b5-manifest-file NOEXTRAS.csv', header=T)
head(annotation)
annotation=as.data.frame(annotation)
names(annotation)[names(annotation) == "Name"] <- "cpg"

cpgdata_bmi=merge(fulldata_bmi, annotation, by="cpg")
cpgdata1_bmi=cpgdata_bmi[order(cpgdata_bmi$q),]
names(cpgdata1_bmi)[names(cpgdata1_bmi) == "objpvalues"] <- "P"
names(cpgdata1_bmi)[names(cpgdata1_bmi) == "MAPINFO"] <- "BP"
names(cpgdata1_bmi)[names(cpgdata1_bmi) == "IlmnID"] <- "SNP"
cpgdata1_bmi$CHR=as.numeric(cpgdata1_bmi$CHR)
cpgdata1_bmi$BP=as.numeric(cpgdata1_bmi$BP)
cpgdata1_bmi[is.na(cpgdata1_bmi)] <- 23
cpgdata1_bmi$P=as.numeric(cpgdata1_bmi$P)

sigcpgdata_bmi=merge(significant_bmi,cpgdata1_bmi, by="cpg")

#Saving off results 
path = '/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/Significant/'
select_path = paste(path, 'BMIannotated', '.csv', sep='')
write.csv(sigcpgdata_bmi, select_path , row.names=F)

#MANHATTAN PLOT abd QQ
par(mfrow=c(1,2))
manhattan(cpgdata1_bmi, ylim = c(0, 10), cex = 0.6, cex.axis = 0.9, suggestiveline = F, genomewideline = F, highlight=sigcpg_bmi, annotateTop = FALSE)
qqman::qq(obj$pvalues)

#Looking up location of geneomic features
genomicfeature=read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/genomeenrich.csv', header=T)
genomicfeature$enhancer=as.factor(genomicfeature$enhancer)
genomicfeature$promoter=as.factor(genomicfeature$promoter)
genomicfeature$shelfshore=as.factor(genomicfeature$shelfshore)
genomicfeature$cpgisland=as.factor(genomicfeature$cpgisland)
genomicfeature$DNAse_site=as.factor(genomicfeature$DNAse_site)

#Breaking out into two tables - significant and nonsignificant CpGs
significant=fulldata_bmi[1:50,]
nonsignificant=fulldata_bmi[51:789656,]
sig_anno=merge(significant,genomicfeature,by="cpg", all.x=T)
nonsig_anno=merge(nonsignificant,genomicfeature,by="cpg", all.x=T)


#Getting counts of genomic features 
table(sig_anno$enhancer)
table(nonsig_anno$enhancer)

table(sig_anno$promoter)
table(nonsig_anno$promoter)

table(sig_anno$shelfshore)
table(nonsig_anno$shelfshore)

table(sig_anno$cpgisland)
table(nonsig_anno$cpgisland)

table(sig_anno$DNAse_site)
table(nonsig_anno$DNAse_site)

#Functional Gene Enrichment - creating 2x2 tables
dat <- data.frame(
  "promoter" = c(10,184399),
  "notpromoter" = c(40, 605021),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat) <- c("Promoter", "Not promoter")

dat2 <- data.frame(
  "island" = c(5,153469),
  "not island" = c(45, 635951),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat2) <- c("island", "Not island")

dat4 <- data.frame(
  "shelf" = c(12,201058),
  "not" = c(38, 588362),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat4) <- c("Shelf", "Not shelf")

dat5 <- data.frame(
  "DNAse" = c(31,470185),
  "not DNAse" = c(19, 319235),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat5) <- c("DNAse", "Not DNAse")

dat6 <- data.frame(
  "enhancer" = c(8,26032),
  "notenhancer" = c(42, 763388),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat6) <- c("Enhancer", "Not enhancer")


#Performing Fisher's exact tests to assess differences 
fisher.test(dat)
fisher.test(dat2)
fisher.test(dat4)
fisher.test(dat5)
fisher.test(dat6)

#EQTM enrichment - bringing in eqtm file from FHS
eqtm=read.csv("/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/eqtm/eqtm_genes.csv", header=T)
mergedsig=merge(eqtm,significant, all.y=T) 
noeqtm_sig=mergedsig[is.na(mergedsig$gene),] #34 have no eqtm
eqtm_sig=mergedsig[!is.na(mergedsig$gene),] 
eqtm_sig2=eqtm_sig[!duplicated(eqtm_sig$cpg), ] # 16 have eqtm 

mergedall=merge(eqtm,nonsignificant, all.y=T) 
noeqtm_all=mergedall[is.na(mergedall$gene),] #752960 have no eqtm 
eqtm_all=mergedall[!is.na(mergedall$gene),] 
eqtm_all2=eqtm_all[!duplicated(eqtm_all$cpg),] #36646

dat7 <- data.frame(
  "eqtm" = c(16, 36646),
  "noteqtm" = c(34, 752960),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat7) <- c("eQTM", "not eQTM")

fisher.test(dat7)


#GO and KEGG ENRICHMENT
ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
allcpgs <- fulldata_bmi[,1]

gst <- gometh(sig.cpg = sigcpg_bmi, all.cpg = allcpgs, 
              collection = "GO", array.type="EPIC", anno = ann)

kegg <- gometh(sig.cpg = sigcpg_bmi, all.cpg = allcpgs, collection = "KEGG", array.type="EPIC", prior.prob=TRUE, anno=ann)

#EQTM bringing in eqtm file from FHS
eqtm=read.csv("/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/eqtm/eqtm_genes.csv", header=T)
mergedcpgs=merge(eqtm,fulldata_bmi)
backgroundgenes=mergedcpgs
mergedsig=merge(backgroundgenes,significant_bmi)

#Deduplicating all datasets#
background_final=backgroundgenes[!duplicated(backgroundgenes$gene), ]
mergedsig_final=mergedsig[!duplicated(mergedsig$gene), ]

universe=background_final$gene
sig_id=mergedsig_final$gene

eg = bitr(universe, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
egsig = bitr(sig_id, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

eg2=eg$ENTREZID
egsig2=egsig$ENTREZID

ego <- enrichGO(gene          = sig_id,
                universe      = universe,
                OrgDb         = org.Hs.eg.db,
                ont           = "all",
                pAdjustMethod = "BH",
                pvalueCutoff = .01,
                qvalueCutoff = .01,
                keyType = "SYMBOL")
summary_sig=data.frame(ego)
s_ego<-clusterProfiler::simplify(ego)
sum_GO=as.data.frame(s_ego)


kk <- enrichKEGG(gene         = egsig2,
                 universe =eg2,
                 organism     = 'hsa',
                 keyType = "kegg",
                 pvalueCutoff = .01,
                 qvalueCutoff=.01)
summary_keg=data.frame(kk)


#Plotting GO and KEGG enrichment results
plot1=dotplot(s_ego, x="count", showCategory=30,font.size = 12, title="(A)       GO_Ontology")
plot2=dotplot(kk, x="count", showCategory=30,font.size = 12, title="(B)      KEGG_Pathway")
plot1+plot2

replicate = read.csv('/net/orion/home/opsasnic/Aim 2 Supplement/replicate.csv', header=T)
cpgdata_izzy=merge(fulldata_bmi, replicate, by="cpg", all.y=T)

