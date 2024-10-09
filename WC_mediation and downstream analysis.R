library(medScan)

###########Reading in exposure-mediator model####################################
all1 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/WC2/New/HRS.aim2.expmed.1.csv', header=T)
all2 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/WC2/New/HRS.aim2.expmed.2.csv', header=T)
all3 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/WC2/New/HRS.aim2.expmed.3.csv', header=T)
all4 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/WC2/New/HRS.aim2.expmed.4.csv', header=T)
all5 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/WC2/New/HRS.aim2.expmed.5.csv', header=T)
all6 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/WC2/New/HRS.aim2.expmed.6.csv', header=T)
all7 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/WC2/New/HRS.aim2.expmed.7.csv', header=T)
all8 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/WC2/New/HRS.aim2.expmed.8.csv', header=T)
all9 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/WC2/New/HRS.aim2.expmed.9.csv', header=T)
all10 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/WC2/New/HRS.aim2.expmed.10.csv', header=T)
all11 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/WC2/New/HRS.aim2.expmed.11.csv', header=T)
all12 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/WC2/New/HRS.aim2.expmed.12.csv', header=T)
all13 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/WC2/New/HRS.aim2.expmed.13.csv', header=T)
all14 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/WC2/New/HRS.aim2.expmed.14.csv', header=T)
all15 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/WC2/New/HRS.aim2.expmed.15.csv', header=T)
all16 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/WC2/New/HRS.aim2.expmed.16.csv', header=T)

ewas_stress3_fixed=rbind(all1, all2, all3, all4, all5, all6, all7, all8, all9, all10, all11, all12, all13, all14, all15,all16)
z.alpha=ewas_stress3_fixed[,5]
p.alpha = ewas_stress3_fixed$pval
names(ewas_stress3_fixed)[names(ewas_stress3_fixed) == "coef"] <- "alpha_coef"


###########Reading in mediator-outcome model####################################
WC1 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/WC/New/HRS.aim2.WC.1.csv', header=T)
WC2 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/WC/New/HRS.aim2.WC.2.csv', header=T)
WC3 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/WC/New/HRS.aim2.WC.3.csv', header=T)
WC4 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/WC/New/HRS.aim2.WC.4.csv', header=T)
WC5 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/WC/New/HRS.aim2.WC.5.csv', header=T)
WC6 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/WC/New/HRS.aim2.WC.6.csv', header=T)
WC7 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/WC/New/HRS.aim2.WC.7.csv', header=T)
WC8 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/WC/New/HRS.aim2.WC.8.csv', header=T)
WC9 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/WC/New/HRS.aim2.WC.9.csv', header=T)
WC10 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/WC/New/HRS.aim2.WC.10.csv', header=T)
WC11 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/WC/New/HRS.aim2.WC.11.csv', header=T)
WC12 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/WC/New/HRS.aim2.WC.12.csv', header=T)
WC13 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/WC/New/HRS.aim2.WC.13.csv', header=T)
WC14 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/WC/New/HRS.aim2.WC.14.csv', header=T)
WC15 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/WC/New/HRS.aim2.WC.15.csv', header=T)
WC16 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/WC/New/HRS.aim2.WC.16.csv', header=T)
ewas_WC_fixed=rbind(WC1, WC2, WC3, WC4, WC5, WC6, WC7, WC8, WC9, WC10,WC11, WC12, WC13, WC14, WC15, WC16)
z.beta=ewas_WC_fixed[,5]
p.beta = ewas_WC_fixed$pval
names(ewas_WC_fixed)[names(ewas_WC_fixed) == "coef"] <- "beta_coef"


###Doing the high dimensional mediation using medScan######
obj = medScan(z.alpha = z.alpha, z.beta = z.beta, method = "HDMT")
objpvalues=obj$pvalues
pvalues_wc=as.data.frame(objpvalues)
pvalues_wc$q=p.adjust(pvalues_wc$objpvalues, method='fdr')
median(qchisq(pvalues_wc$objpvalues, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)

######Merging pvalue file with CPG file################
fulldata_wc=cbind(ewas_stress3_fixed,ewas_WC_fixed,pvalues_wc)
fulldata_wc=fulldata_wc[,-7]
fulldata_wc=fulldata_wc[order(fulldata_wc$q),]


#Getting significant cpg sites from mediation analysis 
significant_wc=subset(fulldata_wc, q<0.05)
sigcpg_wc=significant_wc[,1]


#Merging with annotation file
annotation=read.csv('/net/orion/skardia_lab/clubhouse/research/projects/Health_Retirement_Study/Methylation_Apr_2021/QC_processing_Smith/infinium-methylationepic-v-1-0-b5-manifest-file NOEXTRAS.csv', header=T)
head(annotation)
annotation=as.data.frame(annotation)
names(annotation)[names(annotation) == "Name"] <- "cpg"

cpgdata_wc=merge(fulldata_wc, annotation, by="cpg")
cpgdata1_wc=cpgdata_wc[order(cpgdata_wc$q),]
names(cpgdata1_wc)[names(cpgdata1_wc) == "objpvalues"] <- "P"
names(cpgdata1_wc)[names(cpgdata1_wc) == "MAPINFO"] <- "BP"
names(cpgdata1_wc)[names(cpgdata1_wc) == "IlmnID"] <- "SNP"
cpgdata1_wc$CHR=as.numeric(cpgdata1_wc$CHR)
cpgdata1_wc$BP=as.numeric(cpgdata1_wc$BP)
cpgdata1_wc[is.na(cpgdata1_wc)] <- 23
cpgdata1_wc$P=as.numeric(cpgdata1_wc$P)

sigcpgdata_wc=merge(significant_wc,cpgdata1_wc, by="cpg")

#Saving off results 
path = '/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/Significant/'
select_path = paste(path, 'WCannotated_new', '.csv', sep='')
write.csv(sigcpgdata_wc, select_path , row.names=F)

#MANHATTAN PLOT abd QQ
par(mfrow=c(1,2))
manhattan(cpgdata1_wc, ylim = c(0, 10), cex = 0.6, cex.axis = 0.9, suggestiveline = F, genomewideline = F, highlight=sigcpg_wc, annotateTop = FALSE)
qqman::qq(obj$pvalues)

#Looking up location of genomic features
genomicfeature=read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/genomeenrich.csv', header=T)
genomicfeature$enhancer=as.factor(genomicfeature$enhancer)
genomicfeature$promoter=as.factor(genomicfeature$promoter)
genomicfeature$shelfshore=as.factor(genomicfeature$shelfshore)
genomicfeature$cpgisland=as.factor(genomicfeature$cpgisland)
genomicfeature$DNAse_site=as.factor(genomicfeature$DNAse_site)


#Breaking out into two tables - significant and nonsignificant CpGs
significant=fulldata_wc[1:46,]
nonsignificant=fulldata_wc[47:789656,]
sig_anno=merge(significant,genomicfeature,by="cpg", all.x=T)
nonsig_anno=merge(nonsignificant,genomicfeature,by="cpg", all.x=T)

#Getting counts of genomic features, creating 2x2 tables, and Fisher's exact test 
table(sig_anno$enhancer)
table(nonsig_anno$enhancer)

dat6 <- data.frame(
  "enhancer" = c(8,26032),
  "notenhancer" = c(38, 763392),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat6) <- c("Enhancer", "Not enhancer")

fisher.test(dat6)

table(sig_anno$promoter)
table(nonsig_anno$promoter)

dat <- data.frame(
  "promoter" = c(6,184403),
  "notpromoter" = c(40, 605021),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat) <- c("Promoter", "Not promoter")

fisher.test(dat)


table(sig_anno$shelfshore)
table(nonsig_anno$shelfshore)

dat4 <- data.frame(
  "shelf" = c(12,201058),
  "not" = c(34, 588366),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat4) <- c("Shelf", "Not shelf")

fisher.test(dat4)

table(sig_anno$cpgisland)
table(nonsig_anno$cpgisland)

dat2 <- data.frame(
  "island" = c(4,153470),
  "not island" = c(42, 635954),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat2) <- c("island", "Not island")

fisher.test(dat2)


table(sig_anno$DNAse_site)
table(nonsig_anno$DNAse_site)

dat5 <- data.frame(
  "DNAse" = c(34,470182),
  "not DNAse" = c(12, 319242),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat5) <- c("DNAse", "Not DNAse")

fisher.test(dat5)

#EQTM enrichment - bringing in eqtm file from FHS
eqtm=read.csv("/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/eqtm/eqtm_genes.csv", header=T)
mergedsig=merge(eqtm,significant, all.y=T) 
noeqtm_sig=mergedsig[is.na(mergedsig$gene),] #34 have no eqtm
eqtm_sig=mergedsig[!is.na(mergedsig$gene),] 
eqtm_sig2=eqtm_sig[!duplicated(eqtm_sig$cpg), ] #12  have eqtm 

mergedall=merge(eqtm,nonsignificant, all.y=T) 
noeqtm_all=mergedall[is.na(mergedall$gene),] #752986 have no eqtm 
eqtm_all=mergedall[!is.na(mergedall$gene),] 
eqtm_all2=eqtm_all[!duplicated(eqtm_all$cpg),] #36657

dat7 <- data.frame(
  "eqtm" = c(12,36650),
  "noteqtm" = c(34, 752960),
  row.names = c("Sig", "Not Sig"),
  stringsAsFactors = FALSE
)
colnames(dat7) <- c("eQTM", "not eQTM")

fisher.test(dat7)

#GO and KEGG ENRICHMENT
ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
allcpgs <- fulldata_wc[,1]

gst <- gometh(sig.cpg = sigcpg_wc, all.cpg = allcpgs, 
              collection = "GO", array.type="EPIC", anno = ann)

kegg <- gometh(sig.cpg = sigcpg_wc, all.cpg = allcpgs, collection = "KEGG", array.type="EPIC", prior.prob=TRUE, anno=ann)


#EQTM bringing in eqtm file from FHS
eqtm=read.csv("/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/data/eqtm/eqtm_genes.csv", header=T)
mergedcpgs=merge(eqtm,fulldata_wc)
backgroundgenes=mergedcpgs
mergedsig=merge(backgroundgenes,significant_wc)

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


