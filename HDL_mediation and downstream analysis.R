library(medScan)

###########Reading in exposure-mediator model####################################
all1 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/HDL2/New/HRS.aim2.expmed.1.csv', header=T)
all2 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/HDL2/New/HRS.aim2.expmed.2.csv', header=T)
all3 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/HDL2/New/HRS.aim2.expmed.3.csv', header=T)
all4 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/HDL2/New/HRS.aim2.expmed.4.csv', header=T)
all5 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/HDL2/New/HRS.aim2.expmed.5.csv', header=T)
all6 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/HDL2/New/HRS.aim2.expmed.6.csv', header=T)
all7 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/HDL2/New/HRS.aim2.expmed.7.csv', header=T)
all8 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/HDL2/New/HRS.aim2.expmed.8.csv', header=T)
all9 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/HDL2/New/HRS.aim2.expmed.9.csv', header=T)
all10 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/HDL2/New/HRS.aim2.expmed.10.csv', header=T)
all11 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/HDL2/New/HRS.aim2.expmed.11.csv', header=T)
all12 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/HDL2/New/HRS.aim2.expmed.12.csv', header=T)
all13 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/HDL2/New/HRS.aim2.expmed.13.csv', header=T)
all14 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/HDL2/New/HRS.aim2.expmed.14.csv', header=T)
all15 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/HDL2/New/HRS.aim2.expmed.15.csv', header=T)
all16 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/HDL2/New/HRS.aim2.expmed.16.csv', header=T)

ewas_stress3_fixed=rbind(all1, all2, all3, all4, all5, all6, all7, all8, all9, all10, all11, all12, all13, all14, all15, all16)
z.alpha=ewas_stress3_fixed[,5]
p.alpha = ewas_stress3_fixed$pval
names(ewas_stress3_fixed)[names(ewas_stress3_fixed) == "coef"] <- "alpha_coef"


###########Reading in mediator-outcome model####################################
HDL1 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/HDL/New/HRS.aim2.HDL.1.csv', header=T)
HDL2 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/HDL/New/HRS.aim2.HDL.2.csv', header=T)
HDL3 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/HDL/New/HRS.aim2.HDL.3.csv', header=T)
HDL4 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/HDL/New/HRS.aim2.HDL.4.csv', header=T)
HDL5 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/HDL/New/HRS.aim2.HDL.5.csv', header=T)
HDL6 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/HDL/New/HRS.aim2.HDL.6.csv', header=T)
HDL7 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/HDL/New/HRS.aim2.HDL.7.csv', header=T)
HDL8 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/HDL/New/HRS.aim2.HDL.8.csv', header=T)
HDL9 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/HDL/New/HRS.aim2.HDL.9.csv', header=T)
HDL10 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/HDL/New/HRS.aim2.HDL.10.csv', header=T)
HDL11 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/HDL/New/HRS.aim2.HDL.11.csv', header=T)
HDL12 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/HDL/New/HRS.aim2.HDL.12.csv', header=T)
HDL13 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/HDL/New/HRS.aim2.HDL.13.csv', header=T)
HDL14 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/HDL/New/HRS.aim2.HDL.14.csv', header=T)
HDL15 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/HDL/New/HRS.aim2.HDL.15.csv', header=T)
HDL16 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/HDL/New/HRS.aim2.HDL.16.csv', header=T)

ewas_HDL_fixed=rbind(HDL1, HDL2, HDL3, HDL4, HDL5, HDL6, HDL7, HDL8, HDL9, HDL10,HDL11, HDL12, HDL13, HDL14, HDL15, HDL16)
z.beta=ewas_HDL_fixed[,5]
p.beta = ewas_HDL_fixed$pval
ewas_HDL_fixed$q=p.adjust(ewas_HDL_fixed$pval, method='bonferroni')

names(ewas_HDL_fixed)[names(ewas_HDL_fixed) == "coef"] <- "beta_coef"



###Doing the high dimensional mediation using medScan######
obj = medScan(z.alpha = z.alpha, z.beta = z.beta, method = "HDMT")
objpvalues=obj$pvalues
pvalues_HDL=as.data.frame(objpvalues)
pvalues_HDL$q=p.adjust(pvalues_HDL$objpvalues, method='fdr')
median(qchisq(pvalues_HDL$objpvalues, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)


######Merging pvalue file with CPG file################
fulldata_HDL=cbind(ewas_stress3_fixed,ewas_HDL_fixed,pvalues_HDL)
fulldata_HDL=fulldata_HDL[,-7]
fulldata_HDL=fulldata_HDL[order(fulldata_HDL$q.1),]


#Getting significant cpg sites from mediation analysis 
significant_HDL=subset(fulldata_HDL, q.1<0.05)
sigcpg_HDL=significant_HDL[,1]


#Merging with annotation file
annotation=read.csv('/net/orion/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Methylation_Apr_2021/QC_processing_Smith/infinium-methylationepic-v-1-0-b5-manifest-file NOEXTRAS.csv', header=T)
head(annotation)
annotation=as.data.frame(annotation)
names(annotation)[names(annotation) == "Name"] <- "cpg"

cpgdata_HDL=merge(fulldata_HDL, annotation, by="cpg")
cpgdata1_HDL=cpgdata_HDL[order(cpgdata_HDL$q),]
names(cpgdata1_HDL)[names(cpgdata1_HDL) == "objpvalues"] <- "P"
names(cpgdata1_HDL)[names(cpgdata1_HDL) == "MAPINFO"] <- "BP"
names(cpgdata1_HDL)[names(cpgdata1_HDL) == "IlmnID"] <- "SNP"
cpgdata1_HDL$CHR=as.numeric(cpgdata1_HDL$CHR)
cpgdata1_HDL$BP=as.numeric(cpgdata1_HDL$BP)
cpgdata1_HDL[is.na(cpgdata1_HDL)] <- 23
cpgdata1_HDL$P=as.numeric(cpgdata1_HDL$P)

sigcpgdata_hdl=merge(significant_HDL,cpgdata1_HDL, by="cpg")

#Saving off results 
path = '/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/Significant/'
select_path = paste(path, 'HDLannotated_new', '.csv', sep='')
write.csv(sigcpgdata_hdl, select_path , row.names=F)


#MANHATTAN PLOT abd QQ
par(mfrow=c(1,2))
manhattan(cpgdata1_HDL, ylim = c(0, 10), cex = 0.6, cex.axis = 0.9, suggestiveline = F, genomewideline = F, highlight=sigcpg_HDL, annotateTop = FALSE)
qqman::qq(obj$pvalues)

#Looking up location of genomic features
genomicfeature=read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/genomeenrich.csv', header=T)
genomicfeature$enhancer=as.factor(genomicfeature$enhancer)
genomicfeature$promoter=as.factor(genomicfeature$promoter)
genomicfeature$shelfshore=as.factor(genomicfeature$shelfshore)
genomicfeature$cpgisland=as.factor(genomicfeature$cpgisland)
genomicfeature$DNAse_site=as.factor(genomicfeature$DNAse_site)


#Breaking out into two tables - significant and nonsignificant CpGs
significant=fulldata_HDL[1:7,]
nonsignificant=fulldata_HDL[8:789656,]
sig_anno=merge(significant,genomicfeature,by="cpg", all.x=T)
nonsig_anno=merge(nonsignificant,genomicfeature,by="cpg", all.x=T)

#Getting counts of genomic features, creating 2x2 tables, and Fisher's exact test 
table(sig_anno$enhancer)
table(nonsig_anno$enhancer)

dat6 <- data.frame(
  "enhancer" = c(1,26041),
  "notenhancer" = c(8, 763424),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat6) <- c("Enhancer", "Not enhancer")

fisher.test(dat6)

table(sig_anno$promoter)
table(nonsig_anno$promoter)

dat <- data.frame(
  "promoter" = c(3,184406),
  "notpromoter" = c(4, 605057),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat) <- c("Promoter", "Not promoter")

fisher.test(dat)


table(sig_anno$shelfshore)
table(nonsig_anno$shelfshore)

dat4 <- data.frame(
  "shelf" = c(3,201067),
  "not" = c(4, 588396),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat4) <- c("Shelf", "Not shelf")

fisher.test(dat4)

table(sig_anno$cpgisland)
table(nonsig_anno$cpgisland)

dat2 <- data.frame(
  "island" = c(1,153475),
  "not island" = c(8, 635990),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat2) <- c("island", "Not island")

fisher.test(dat2)


table(sig_anno$DNAse_site)
table(nonsig_anno$DNAse_site)

dat5 <- data.frame(
  "DNAse" = c(8,470210),
  "not DNAse" = c(1, 319255),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat5) <- c("DNAse", "Not DNAse")

fisher.test(dat5)

#EQTM enrichment - bringing in eqtm file from FHS
eqtm=read.csv("/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/eqtm/eqtm_genes.csv", header=T)
mergedsig=merge(eqtm,significant, all.y=T) 
noeqtm_sig=mergedsig[is.na(mergedsig$gene),] #4 have no eqtm
eqtm_sig=mergedsig[!is.na(mergedsig$gene),] 
eqtm_sig2=eqtm_sig[!duplicated(eqtm_sig$cpg), ] # 3 have eqtm 

mergedall=merge(eqtm,nonsignificant, all.y=T) 
noeqtm_all=mergedall[is.na(mergedall$gene),] #752989 have no eqtm 
eqtm_all=mergedall[!is.na(mergedall$gene),] 
eqtm_all2=eqtm_all[!duplicated(eqtm_all$cpg),] #36660


#Functional Gene Enrichment - creating 2x2 tables
dat7 <- data.frame(
  "eqtm" = c(2,36660),
  "noteqtm" = c(5, 752989),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat7) <- c("eQTM", "not eQTM")

fisher.test(dat7)


#GO and KEGG ENRICHMENT
ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
allcpgs <- fulldata_HDL[,1]

gst <- gometh(sig.cpg = sigcpg_HDL, all.cpg = allcpgs, 
              collection = "GO", array.type="EPIC", anno = ann)

kegg <- gometh(sig.cpg = sigcpg_HDL, all.cpg = allcpgs, collection = "KEGG", array.type="EPIC", prior.prob=TRUE, anno=ann)


#EQTM bringing in eqtm file from FHS
eqtm=read.csv("/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/eqtm/eqtm_genes.csv", header=T)
mergedcpgs=merge(eqtm,fulldata_HDL)
backgroundgenes=mergedcpgs
mergedsig=merge(backgroundgenes,significant_HDL)

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
dotplot(s_ego, x="count", showCategory=30,font.size = 12, title="                               GO_Ontology")
