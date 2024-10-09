library(medScan)

###########Reading in exposure-mediator model####################################
all1 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/CRP2/New/HRS.aim2.expmed.1.csv', header=T)
all2 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/CRP2/New/HRS.aim2.expmed.2.csv', header=T)
all3 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/CRP2/New/HRS.aim2.expmed.3.csv', header=T)
all4 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/CRP2/New/HRS.aim2.expmed.4.csv', header=T)
all5 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/CRP2/New/HRS.aim2.expmed.5.csv', header=T)
all6 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/CRP2/New/HRS.aim2.expmed.6.csv', header=T)
all7 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/CRP2/New/HRS.aim2.expmed.7.csv', header=T)
all8 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/CRP2/New/HRS.aim2.expmed.8.csv', header=T)
all9 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/CRP2/New/HRS.aim2.expmed.9.csv', header=T)
all10 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/CRP2/New/HRS.aim2.expmed.10.csv', header=T)
all11 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/CRP2/New/HRS.aim2.expmed.11.csv', header=T)
all12 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/CRP2/New/HRS.aim2.expmed.12.csv', header=T)
all13 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/CRP2/New/HRS.aim2.expmed.13.csv', header=T)
all14 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/CRP2/New/HRS.aim2.expmed.14.csv', header=T)
all15 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/CRP2/New/HRS.aim2.expmed.15.csv', header=T)
all16 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/CRP2/New/HRS.aim2.expmed.16.csv', header=T)

ewas_stress3_fixed=rbind(all1, all2, all3, all4, all5, all6, all7, all8, all9, all10, all11, all12, all13, all14, all15, all16)
z.alpha=ewas_stress3_fixed[,5]
p.alpha = ewas_stress3_fixed$pval
names(ewas_stress3_fixed)[names(ewas_stress3_fixed) == "coef"] <- "alpha_coef"

###########Reading in mediator-outcome model####################################
CRP1 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/logCRP/New/HRS.aim2.CRP.1.csv', header=T)
CRP2 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/logCRP/New/HRS.aim2.CRP.2.csv', header=T)
CRP3 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/logCRP/New/HRS.aim2.CRP.3.csv', header=T)
CRP4 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/logCRP/New/HRS.aim2.CRP.4.csv', header=T)
CRP5 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/logCRP/New/HRS.aim2.CRP.5.csv', header=T)
CRP6 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/logCRP/New/HRS.aim2.CRP.6.csv', header=T)
CRP7 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/logCRP/New/HRS.aim2.CRP.7.csv', header=T)
CRP8 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/logCRP/New/HRS.aim2.CRP.8.csv', header=T)
CRP9 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/logCRP/New/HRS.aim2.CRP.9.csv', header=T)
CRP10 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/logCRP/New/HRS.aim2.CRP.10.csv', header=T)
CRP11 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/logCRP/New/HRS.aim2.CRP.11.csv', header=T)
CRP12 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/logCRP/New/HRS.aim2.CRP.12.csv', header=T)
CRP13 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/logCRP/New/HRS.aim2.CRP.13.csv', header=T)
CRP14 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/logCRP/New/HRS.aim2.CRP.14.csv', header=T)
CRP15 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/logCRP/New/HRS.aim2.CRP.15.csv', header=T)
CRP16 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/logCRP/New/HRS.aim2.CRP.16.csv', header=T)
ewas_CRP_fixed=rbind(CRP1, CRP2, CRP3, CRP4, CRP5, CRP6, CRP7, CRP8, CRP9, CRP10,CRP11, CRP12, CRP13, CRP14, CRP15, CRP16)
z.beta=ewas_CRP_fixed[,5]
p.beta = ewas_CRP_fixed$pval
names(ewas_CRP_fixed)[names(ewas_CRP_fixed) == "coef"] <- "beta_coef"


###Doing the high dimensional mediation using medScan######
obj = medScan(z.alpha = z.alpha, z.beta = z.beta, method = "HDMT")
objpvalues=obj$pvalues
pvalues_crp=as.data.frame(objpvalues)
pvalues_crp$q=p.adjust(pvalues_crp$objpvalues, method='fdr')
median(qchisq(pvalues_crp$objpvalues, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)


######Merging pvalue file with CPG file################
fulldata_crp=cbind(ewas_stress3_fixed,ewas_CRP_fixed,pvalues_crp)
fulldata_crp=fulldata_crp[,-7]
fulldata_crp=fulldata_crp[order(fulldata_crp$q),]

#Getting significant cpg sites from mediation analysis 
significant_crp=subset(fulldata_crp, q<0.05)
sigcpg_crp=significant_crp[,1]

#Merging with annotation file
annotation=read.csv('/net/orion/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Methylation_Apr_2021/QC_processing_Smith/infinium-methylationepic-v-1-0-b5-manifest-file NOEXTRAS.csv', header=T)
head(annotation)
annotation=as.data.frame(annotation)
names(annotation)[names(annotation) == "Name"] <- "cpg"

cpgdata_crp=merge(fulldata_crp, annotation, by="cpg")
cpgdata1_crp=cpgdata_crp[order(cpgdata_crp$q),]
names(cpgdata1_crp)[names(cpgdata1_crp) == "objpvalues"] <- "P"
names(cpgdata1_crp)[names(cpgdata1_crp) == "MAPINFO"] <- "BP"
names(cpgdata1_crp)[names(cpgdata1_crp) == "IlmnID"] <- "SNP"
cpgdata1_crp$CHR=as.numeric(cpgdata1_crp$CHR)
cpgdata1_crp$BP=as.numeric(cpgdata1_crp$BP)
cpgdata1_crp[is.na(cpgdata1_crp)] <- 23
cpgdata1_crp$P=as.numeric(cpgdata1_crp$P)

sigcpgdata_crp=merge(significant_crp,cpgdata1_crp, by="cpg")

#Saving off results 
path = '/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/Significant/'
select_path = paste(path, 'CRPannotated_new', '.csv', sep='')
write.csv(sigcpgdata_crp, select_path , row.names=F)

#MANHATTAN PLOT abd QQ
par(mfrow=c(1,2))
manhattan(cpgdata1_crp, ylim = c(0, 10), cex = 0.6, cex.axis = 0.9, suggestiveline = F, genomewideline = F, highlight=sigcpg_crp, annotateTop = FALSE)
qqman::qq(obj$pvalues)

#Looking up location of genomic features
genomicfeature=read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/genomeenrich.csv', header=T)
genomicfeature$enhancer=as.factor(genomicfeature$enhancer)
genomicfeature$promoter=as.factor(genomicfeature$promoter)
genomicfeature$shelfshore=as.factor(genomicfeature$shelfshore)
genomicfeature$cpgisland=as.factor(genomicfeature$cpgisland)
genomicfeature$DNAse_site=as.factor(genomicfeature$DNAse_site)


#Breaking out into two tables - significant and nonsignificant CpGs
significant=fulldata_crp[1:12,]
nonsignificant=fulldata_crp[13:789656,]
sig_anno=merge(significant,genomicfeature,by="cpg", all.x=T)
nonsig_anno=merge(nonsignificant,genomicfeature,by="cpg", all.x=T)

#Getting counts of genomic features and creating 2x2 tables
table(sig_anno$promoter)
table(nonsig_anno$promoter)

dat <- data.frame(
  "promoter" = c(2,184407),
  "notpromoter" = c(10, 605051),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat) <- c("Promoter", "Not promoter")


table(sig_anno$cpgisland)
table(nonsig_anno$cpgisland)

dat2 <- data.frame(
  "island" = c(1,153475),
  "not island" = c(13, 635985),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat2) <- c("island", "Not island")

table(sig_anno$shelfshore)
table(nonsig_anno$shelfshore)

dat4 <- data.frame(
  "shelf" = c(5,201065),
  "not" = c(7, 588393),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat4) <- c("Shelf", "Not shelf")

table(sig_anno$DNAse_site)
table(nonsig_anno$DNAse_site)

dat5 <- data.frame(
  "DNAse" = c(9,470207),
  "not DNAse" = c(3, 319251),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat5) <- c("DNAse", "Not DNAse")

table(sig_anno$enhancer)
table(nonsig_anno$enhancer)

dat6 <- data.frame(
  "enhancer" = c(3,26037),
  "notenhancer" = c(9, 763421),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat6) <- c("Enhancer", "Not enhancer")

#Performing Fisher's exact tests to assess differences 

fisher.test(dat2)
fisher.test(dat)
fisher.test(dat4)
fisher.test(dat5)
fisher.test(dat6)


#EQTM enrichment - bringing in eqtm file from FHS
eqtm=read.csv("/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/eqtm/eqtm_genes.csv", header=T)
mergedsig=merge(eqtm,significant, all.y=T) 
noeqtm_sig=mergedsig[is.na(mergedsig$gene),] #7 have no eqtm
eqtm_sig=mergedsig[!is.na(mergedsig$gene),] 
eqtm_sig2=eqtm_sig[!duplicated(eqtm_sig$cpg), ] # 5 have eqtm 

mergedall=merge(eqtm,nonsignificant, all.y=T) 
noeqtm_all=mergedall[is.na(mergedall$gene),] #752987 have no eqtm 
eqtm_all=mergedall[!is.na(mergedall$gene),]
eqtm_all2=eqtm_all[!duplicated(eqtm_all$cpg),] #36657


mergedall=merge(eqtm,fulldata_crp, all.y=T) 
noeqtm_all=mergedall[is.na(mergedall$gene),] #752994 have no eqtm 
eqtm_all=mergedall[!is.na(mergedall$gene),]
eqtm_all2=eqtm_all[!duplicated(eqtm_all$cpg),] #36662

dat7 <- data.frame(
  "eqtm" = c(5,36657),
  "noteqtm" = c(7, 752987),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat7) <- c("eQTM", "not eQTM")

fisher.test(dat7)

#GO and KEGG ENRICHMENT
ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
allcpgs <- fulldata_crp[,1]

gst <- gometh(sig.cpg = sigcpg_crp, all.cpg = allcpgs, 
              collection = "GO", array.type="EPIC", anno = ann)

kegg <- gometh(sig.cpg = sigcpg_crp, all.cpg = allcpgs, collection = "KEGG", array.type="EPIC", prior.prob=TRUE, anno=ann)


#EQTM bringing in eqtm file from FHS
eqtm=read.csv("/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/eqtm/eqtm_genes.csv", header=T)
mergedcpgs=merge(eqtm,fulldata_crp)
backgroundgenes=mergedcpgs
mergedsig=merge(backgroundgenes,significant_crp)

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
                 pvalueCutoff = .05,
                 qvalueCutoff=.05)
summary_keg=data.frame(kk)


#Plotting GO and KEGG enrichment results
plot1=dotplot(s_ego, x="count", showCategory=30,font.size = 12, title="(A)        GO_Ontology")
plot2=dotplot(kk, x="count", showCategory=30,font.size = 12, title="(B)      KEGG_Pathway")
plot1+plot2


