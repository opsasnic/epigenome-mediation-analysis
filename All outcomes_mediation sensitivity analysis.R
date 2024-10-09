library(medScan)
library(HDMT)



####################BMI - adjusted for smoking ###################################

###########Reading in exposure-mediator model####################################
all1 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/BMI2/NewAdj/HRS.aim2.expmed.1.csv', header=T)
all2 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/BMI2/NewAdj/HRS.aim2.expmed.2.csv', header=T)
all3 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/BMI2/NewAdj/HRS.aim2.expmed.3.csv', header=T)
all4 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/BMI2/NewAdj/HRS.aim2.expmed.4.csv', header=T)
all5 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/BMI2/NewAdj/HRS.aim2.expmed.5.csv', header=T)
all6 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/BMI2/NewAdj/HRS.aim2.expmed.6.csv', header=T)
all7 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/BMI2/NewAdj/HRS.aim2.expmed.7.csv', header=T)
all8 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/BMI2/NewAdj/HRS.aim2.expmed.8.csv', header=T)
all9 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/BMI2/NewAdj/HRS.aim2.expmed.9.csv', header=T)
all10 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/BMI2/NewAdj/HRS.aim2.expmed.10.csv', header=T)
all11 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/BMI2/NewAdj/HRS.aim2.expmed.11.csv', header=T)
all12 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/BMI2/NewAdj/HRS.aim2.expmed.12.csv', header=T)
all13 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/BMI2/NewAdj/HRS.aim2.expmed.13.csv', header=T)
all14 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/BMI2/NewAdj/HRS.aim2.expmed.14.csv', header=T)
all15 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/BMI2/NewAdj/HRS.aim2.expmed.15.csv', header=T)
all16 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/BMI2/NewAdj/HRS.aim2.expmed.16.csv', header=T)

ewas_stress3_fixed=rbind(all1, all2, all3, all4, all5, all6, all7, all8, all9, all10, all11, all12, all13, all14, all15, all16)
z.alpha=ewas_stress3_fixed[,5]
p.alpha = ewas_stress3_fixed$pval
names(ewas_stress3_fixed)[names(ewas_stress3_fixed) == "coef"] <- "alpha_coef"

###########Reading in mediator-outcome model####################################
BMI1 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/BMI/NewAdj/HRS.aim2.BMIadj.1.csv', header=T)
BMI2 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/BMI/NewAdj/HRS.aim2.BMIadj.2.csv', header=T)
BMI3 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/BMI/NewAdj/HRS.aim2.BMIadj.3.csv', header=T)
BMI4 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/BMI/NewAdj/HRS.aim2.BMIadj.4.csv', header=T)
BMI5 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/BMI/NewAdj/HRS.aim2.BMIadj.5.csv', header=T)
BMI6 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/BMI/NewAdj/HRS.aim2.BMIadj.6.csv', header=T)
BMI7 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/BMI/NewAdj/HRS.aim2.BMIadj.7.csv', header=T)
BMI8 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/BMI/NewAdj/HRS.aim2.BMIadj.8.csv', header=T)
BMI9 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/BMI/NewAdj/HRS.aim2.BMIadj.9.csv', header=T)
BMI10 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/BMI/NewAdj/HRS.aim2.BMIadj.10.csv', header=T)
BMI11 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/BMI/NewAdj/HRS.aim2.BMIadj.11.csv', header=T)
BMI12 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/BMI/NewAdj/HRS.aim2.BMIadj.12.csv', header=T)
BMI13 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/BMI/NewAdj/HRS.aim2.BMIadj.13.csv', header=T)
BMI14 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/BMI/NewAdj/HRS.aim2.BMIadj.14.csv', header=T)
BMI15 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/BMI/NewAdj/HRS.aim2.BMIadj.15.csv', header=T)
BMI16 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/BMI/NewAdj/HRS.aim2.BMIadj.16.csv', header=T)

ewas_bmi_fixed=rbind(BMI1, BMI2, BMI3, BMI4, BMI5, BMI6, BMI7, BMI8, BMI9, BMI10,BMI11, BMI12, BMI13, BMI14, BMI15, BMI16)
z.beta=ewas_bmi_fixed[,5]
p.beta = ewas_bmi_fixed$pval
names(ewas_bmi_fixed)[names(ewas_bmi_fixed) == "coef"] <- "beta_coef"

###Doing the high dimensional mediation using medScan######
obj = medScan(z.alpha = z.alpha, z.beta = z.beta, method = "HDMT")
objpvalues=obj$pvalues
pvalues_bmi=as.data.frame(objpvalues)
pvalues_bmi$q=p.adjust(pvalues_bmi$objpvalues, method='fdr')

######Merging pvalue file with CPG file################
fulldata_bmiadj=cbind(ewas_stress3_fixed,ewas_bmi_fixed,pvalues_bmi)
fulldata_bmiadj=fulldata_bmiadj[,-7]

#Getting significant cpg sites from mediation analysis 
significant_bmiadj=subset(fulldata_bmiadj, q<0.05)
significant_bmiadj=significant_bmiadj[,1]

#Extracting CpGs that were FDR significant in Model 1
bmi_adj=subset(fulldata_bmiadj, cpg =="cg00420390" |
                 cpg =="cg00851028" |
                 cpg =="cg01406381" |
                 cpg =="cg02097120" |
                 cpg =="cg02120866" |
                 cpg =="cg02370100" |
                 cpg =="cg02508743" |
                 cpg =="cg03483314" |
                 cpg =="cg03699074" |
                 cpg =="cg04257841" |
                 cpg =="cg04803208" |
                 cpg =="cg04881642" |
                 cpg =="cg07321742" |
                 cpg =="cg07699300" |
                 cpg =="cg08274633" |
                 cpg =="cg08366476" |
                 cpg =="cg09315878" |
                 cpg =="cg09885325" |
                 cpg =="cg09971499" |
                 cpg =="cg11454468" |
                 cpg =="cg11961918" |
                 cpg =="cg12582959" |
                 cpg =="cg12693436" |
                 cpg =="cg12937391" |
                 cpg =="cg13509283" |
                 cpg =="cg13770461" |
                 cpg =="cg13808803" |
                 cpg =="cg13843547" |
                 cpg =="cg15095917" |
                 cpg =="cg15527308" |
                 cpg =="cg15781610" |
                 cpg =="cg16593554" |
                 cpg =="cg16760163" |
                 cpg =="cg17489748" |
                 cpg =="cg18129971" |
                 cpg =="cg18944984" |
                 cpg =="cg19748455" |
                 cpg =="cg20097985" |
                 cpg =="cg20138077" |
                 cpg =="cg20445774" |
                 cpg =="cg20603452" |
                 cpg =="cg21327712" |
                 cpg =="cg22108563" |
                 cpg =="cg22540135" |
                 cpg =="cg23111342" |
                 cpg =="cg23919111" |
                 cpg =="cg24990400" |
                 cpg =="cg25607249" |
                 cpg =="cg27004870" |
                 cpg =="cg27244120")



################################ WC - adjusted for smoking ########################3

###########Reading in exposure-mediator model####################################
all1 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/WC2/NewAdj/HRS.aim2.expmed.1.csv', header=T)
all2 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/WC2/NewAdj/HRS.aim2.expmed.2.csv', header=T)
all3 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/WC2/NewAdj/HRS.aim2.expmed.3.csv', header=T)
all4 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/WC2/NewAdj/HRS.aim2.expmed.4.csv', header=T)
all5 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/WC2/NewAdj/HRS.aim2.expmed.5.csv', header=T)
all6 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/WC2/NewAdj/HRS.aim2.expmed.6.csv', header=T)
all7 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/WC2/NewAdj/HRS.aim2.expmed.7.csv', header=T)
all8 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/WC2/NewAdj/HRS.aim2.expmed.8.csv', header=T)
all9 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/WC2/NewAdj/HRS.aim2.expmed.9.csv', header=T)
all10 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/WC2/NewAdj/HRS.aim2.expmed.10.csv', header=T)
all11 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/WC2/NewAdj/HRS.aim2.expmed.11.csv', header=T)
all12 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/WC2/NewAdj/HRS.aim2.expmed.12.csv', header=T)
all13 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/WC2/NewAdj/HRS.aim2.expmed.13.csv', header=T)
all14 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/WC2/NewAdj/HRS.aim2.expmed.14.csv', header=T)
all15 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/WC2/NewAdj/HRS.aim2.expmed.15.csv', header=T)
all16 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/WC2/NewAdj/HRS.aim2.expmed.16.csv', header=T)

ewas_stress3_fixed=rbind(all1, all2, all3, all4, all5, all6, all7, all8, all9, all10, all11, all12, all13, all14, all15,all16)
z.alpha=ewas_stress3_fixed[,5]
p.alpha = ewas_stress3_fixed$pval
names(ewas_stress3_fixed)[names(ewas_stress3_fixed) == "coef"] <- "alpha_coef"

###########Reading in mediator-outcome model####################################
WC1 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/WC/NewAdj/HRS.aim2.WCadj.1.csv', header=T)
WC2 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/WC/NewAdj/HRS.aim2.WCadj.2.csv', header=T)
WC3 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/WC/NewAdj/HRS.aim2.WCadj.3.csv', header=T)
WC4 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/WC/NewAdj/HRS.aim2.WCadj.4.csv', header=T)
WC5 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/WC/NewAdj/HRS.aim2.WCadj.5.csv', header=T)
WC6 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/WC/NewAdj/HRS.aim2.WCadj.6.csv', header=T)
WC7 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/WC/NewAdj/HRS.aim2.WCadj.7.csv', header=T)
WC8 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/WC/NewAdj/HRS.aim2.WCadj.8.csv', header=T)
WC9 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/WC/NewAdj/HRS.aim2.WCadj.9.csv', header=T)
WC10 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/WC/NewAdj/HRS.aim2.WCadj.10.csv', header=T)
WC11 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/WC/NewAdj/HRS.aim2.WCadj.11.csv', header=T)
WC12 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/WC/NewAdj/HRS.aim2.WCadj.12.csv', header=T)
WC13 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/WC/NewAdj/HRS.aim2.WCadj.13.csv', header=T)
WC14 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/WC/NewAdj/HRS.aim2.WCadj.14.csv', header=T)
WC15 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/WC/NewAdj/HRS.aim2.WCadj.15.csv', header=T)
WC16 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/WC/NewAdj/HRS.aim2.WCadj.16.csv', header=T)
ewas_WC_fixed=rbind(WC1, WC2, WC3, WC4, WC5, WC6, WC7, WC8, WC9, WC10,WC11, WC12, WC13, WC14, WC15, WC16)
z.beta=ewas_WC_fixed[,5]
p.beta = ewas_WC_fixed$pval
names(ewas_WC_fixed)[names(ewas_WC_fixed) == "coef"] <- "beta_coef"

###Doing the high dimensional mediation using medScan######
obj = medScan(z.alpha = z.alpha, z.beta = z.beta, method = "HDMT")
objpvalues=obj$pvalues
pvalues_wcadjfix=as.data.frame(objpvalues)
pvalues_wcadjfix$q=p.adjust(pvalues_wcadjfix$objpvalues, method='fdr')

######Merging pvalue file with CPG file################
fulldata_wcadj=cbind(ewas_stress3_fixed,ewas_WC_fixed,pvalues_wcadjfix)
fulldata_wcadj=fulldata_wcadj[,-7]


#Getting significant cpg sites from mediation analysis 
significant_wcadj=subset(fulldata_wcadj, q<0.05)
sigcpg_wcadj=significant_wcadj[,1]


#Extracting CpGs that were FDR significant in Model 1
wc_adj=subset(fulldata_wcadj, 
              cpg =="cg00420390" |
                cpg =="cg00950718" |
                cpg =="cg01406381" |
                cpg =="cg01695954" | 
                cpg =="cg02097120" | 
                cpg =="cg02370100" | 
                cpg =="cg03037271" | 
                cpg =="cg04051206" | 
                cpg =="cg04803208" | 
                cpg =="cg05241923" | 
                cpg =="cg05308643" | 
                cpg =="cg06088445" | 
                cpg =="cg06809522" | 
                cpg =="cg07321742" | 
                cpg =="cg07699300" | 
                cpg =="cg07787851" | 
                cpg =="cg08848140" | 
                cpg =="cg09315878" | 
                cpg =="cg09971499" | 
                cpg =="cg10529845" | 
                cpg =="cg11454468" | 
                cpg =="cg12581512" | 
                cpg =="cg12582959" | 
                cpg =="cg12693436" | 
                cpg =="cg13509283" | 
                cpg =="cg13843547" | 
                cpg =="cg15527308" | 
                cpg =="cg15781610" | 
                cpg =="cg16593554" | 
                cpg =="cg19748455" | 
                cpg =="cg20097985" | 
                cpg =="cg20138077" | 
                cpg =="cg20603452" | 
                cpg =="cg20745693" | 
                cpg =="cg21646392" | 
                cpg =="cg22540135" | 
                cpg =="cg23086613" | 
                cpg =="cg23691006" | 
                cpg =="cg23733945" | 
                cpg =="cg24507742" | 
                cpg =="cg24704287" | 
                cpg =="cg24837149" | 
                cpg =="cg25607249" | 
                cpg =="cg25922180" | 
                cpg =="cg27004870" |
                cpg == "cg27244120") 


############################HDL-C adjusted for smoking ###################################
###########Reading in exposure-mediator model####################################
all1 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/HDL2/NewAdj/HRS.aim2.expmed.1.csv', header=T)
all2 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/HDL2/NewAdj/HRS.aim2.expmed.2.csv', header=T)
all3 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/HDL2/NewAdj/HRS.aim2.expmed.3.csv', header=T)
all4 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/HDL2/NewAdj/HRS.aim2.expmed.4.csv', header=T)
all5 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/HDL2/NewAdj/HRS.aim2.expmed.5.csv', header=T)
all6 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/HDL2/NewAdj/HRS.aim2.expmed.6.csv', header=T)
all7 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/HDL2/NewAdj/HRS.aim2.expmed.7.csv', header=T)
all8 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/HDL2/NewAdj/HRS.aim2.expmed.8.csv', header=T)
all9 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/HDL2/NewAdj/HRS.aim2.expmed.9.csv', header=T)
all10 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/HDL2/NewAdj/HRS.aim2.expmed.10.csv', header=T)
all11 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/HDL2/NewAdj/HRS.aim2.expmed.11.csv', header=T)
all12 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/HDL2/NewAdj/HRS.aim2.expmed.12.csv', header=T)
all13 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/HDL2/NewAdj/HRS.aim2.expmed.13.csv', header=T)
all14 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/HDL2/NewAdj/HRS.aim2.expmed.14.csv', header=T)
all15 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/HDL2/NewAdj/HRS.aim2.expmed.15.csv', header=T)
all16 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/HDL2/NewAdj/HRS.aim2.expmed.16.csv', header=T)

ewas_stress3_fixed=rbind(all1, all2, all3, all4, all5, all6, all7, all8, all9, all10, all11, all12, all13, all14, all15, all16)
z.alpha=ewas_stress3_fixed[,5]
p.alpha = ewas_stress3_fixed$pval
names(ewas_stress3_fixed)[names(ewas_stress3_fixed) == "coef"] <- "alpha_coef"

###########Reading in mediator-outcome model####################################
HDL1 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/HDL/NewAdj/HRS.aim2.HDLadj.1.csv', header=T)
HDL2 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/HDL/NewAdj/HRS.aim2.HDLadj.2.csv', header=T)
HDL3 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/HDL/NewAdj/HRS.aim2.HDLadj.3.csv', header=T)
HDL4 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/HDL/NewAdj/HRS.aim2.HDLadj.4.csv', header=T)
HDL5 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/HDL/NewAdj/HRS.aim2.HDLadj.5.csv', header=T)
HDL6 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/HDL/NewAdj/HRS.aim2.HDLadj.6.csv', header=T)
HDL7 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/HDL/NewAdj/HRS.aim2.HDLadj.7.csv', header=T)
HDL8 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/HDL/NewAdj/HRS.aim2.HDLadj.8.csv', header=T)
HDL9 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/HDL/NewAdj/HRS.aim2.HDLadj.9.csv', header=T)
HDL10 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/HDL/NewAdj/HRS.aim2.HDLadj.10.csv', header=T)
HDL11 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/HDL/NewAdj/HRS.aim2.HDLadj.11.csv', header=T)
HDL12 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/HDL/NewAdj/HRS.aim2.HDLadj.12.csv', header=T)
HDL13 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/HDL/NewAdj/HRS.aim2.HDLadj.13.csv', header=T)
HDL14 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/HDL/NewAdj/HRS.aim2.HDLadj.14.csv', header=T)
HDL15 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/HDL/NewAdj/HRS.aim2.HDLadj.15.csv', header=T)
HDL16 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/HDL/NewAdj/HRS.aim2.HDLadj.16.csv', header=T)

ewas_HDL_fixed=rbind(HDL1, HDL2, HDL3, HDL4, HDL5, HDL6, HDL7, HDL8, HDL9, HDL10,HDL11, HDL12, HDL13, HDL14, HDL15, HDL16)
z.beta=ewas_HDL_fixed[,5]
p.beta = ewas_HDL_fixed$pval
names(ewas_HDL_fixed)[names(ewas_HDL_fixed) == "coef"] <- "beta_coef"

###Doing the high dimensional mediation using medScan######
obj = medScan(z.alpha = z.alpha, z.beta = z.beta, method = "HDMT")
objpvalues=obj$pvalues
pvalues_HDLadj=as.data.frame(objpvalues)
pvalues_HDLadj$q=p.adjust(pvalues_HDLadj$objpvalues, method='fdr')

######Merging pvalue file with CPG file################
fulldata_HDLadj=cbind(ewas_stress3_fixed,ewas_HDL_fixed,pvalues_HDLadj)
fulldata_HDLadj=fulldata_HDLadj[,-7]


#Getting significant cpg sites from mediation analysis 
significant_HDLadj=subset(fulldata_HDLadj, q<0.05)
sigcpg_HDLadj=significant_HDLadj[,1]


#Extracting CpGs that were FDR significant in Model 1
hdl_adj=subset(fulldata_HDLadj, cpg == "cg01695954" | cpg == "cg04051649" | cpg == "cg24507742" | cpg == "cg04803208" | cpg == "cg06809522" | cpg == "cg08274633" | cpg == "cg13843547" | cpg == "cg25607249")



################################# CRP - adjusted for smoking ##############################
###########Reading in exposure-mediator model####################################
all1 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/CRP2/NewAdj/HRS.aim2.expmed.1.csv', header=T)
all2 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/CRP2/NewAdj/HRS.aim2.expmed.2.csv', header=T)
all3 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/CRP2/NewAdj/HRS.aim2.expmed.3.csv', header=T)
all4 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/CRP2/NewAdj/HRS.aim2.expmed.4.csv', header=T)
all5 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/CRP2/NewAdj/HRS.aim2.expmed.5.csv', header=T)
all6 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/CRP2/NewAdj/HRS.aim2.expmed.6.csv', header=T)
all7 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/CRP2/NewAdj/HRS.aim2.expmed.7.csv', header=T)
all8 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/CRP2/NewAdj/HRS.aim2.expmed.8.csv', header=T)
all9 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/CRP2/NewAdj/HRS.aim2.expmed.9.csv', header=T)
all10 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/CRP2/NewAdj/HRS.aim2.expmed.10.csv', header=T)
all11 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/CRP2/NewAdj/HRS.aim2.expmed.11.csv', header=T)
all12 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/CRP2/NewAdj/HRS.aim2.expmed.12.csv', header=T)
all13 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/CRP2/NewAdj/HRS.aim2.expmed.13.csv', header=T)
all14 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/CRP2/NewAdj/HRS.aim2.expmed.14.csv', header=T)
all15 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/CRP2/NewAdj/HRS.aim2.expmed.15.csv', header=T)
all16 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/ExposureMediator/CRP2/NewAdj/HRS.aim2.expmed.16.csv', header=T)

ewas_stress3_fixed=rbind(all1, all2, all3, all4, all5, all6, all7, all8, all9, all10, all11, all12, all13, all14, all15, all16)
z.alpha=ewas_stress3_fixed[,5]
p.alpha = ewas_stress3_fixed$pval
names(ewas_stress3_fixed)[names(ewas_stress3_fixed) == "coef"] <- "alpha_coef"


###########Reading in mediator-outcome model####################################
CRP1 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/logCRP/NewAdj/HRS.aim2.CRPadj.1.csv', header=T)
CRP2 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/logCRP/NewAdj/HRS.aim2.CRPadj.2.csv', header=T)
CRP3 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/logCRP/NewAdj/HRS.aim2.CRPadj.3.csv', header=T)
CRP4 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/logCRP/NewAdj/HRS.aim2.CRPadj.4.csv', header=T)
CRP5 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/logCRP/NewAdj/HRS.aim2.CRPadj.5.csv', header=T)
CRP6 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/logCRP/NewAdj/HRS.aim2.CRPadj.6.csv', header=T)
CRP7 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/logCRP/NewAdj/HRS.aim2.CRPadj.7.csv', header=T)
CRP8 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/logCRP/NewAdj/HRS.aim2.CRPadj.8.csv', header=T)
CRP9 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/logCRP/NewAdj/HRS.aim2.CRPadj.9.csv', header=T)
CRP10 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/logCRP/NewAdj/HRS.aim2.CRPadj.10.csv', header=T)
CRP11 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/logCRP/NewAdj/HRS.aim2.CRPadj.11.csv', header=T)
CRP12 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/logCRP/NewAdj/HRS.aim2.CRPadj.12.csv', header=T)
CRP13 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/logCRP/NewAdj/HRS.aim2.CRPadj.13.csv', header=T)
CRP14 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/logCRP/NewAdj/HRS.aim2.CRPadj.14.csv', header=T)
CRP15 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/logCRP/NewAdj/HRS.aim2.CRPadj.15.csv', header=T)
CRP16 = read.csv('/net/orion/skardia_lab/treehouse/science/projects/Health_Retirement_Study/Epigenetics/Stress/results/Aim2/logCRP/NewAdj/HRS.aim2.CRPadj.16.csv', header=T)
ewas_CRP_fixed=rbind(CRP1, CRP2, CRP3, CRP4, CRP5, CRP6, CRP7, CRP8, CRP9, CRP10,CRP11, CRP12, CRP13, CRP14, CRP15, CRP16)
z.beta=ewas_CRP_fixed[,5]
p.beta = ewas_CRP_fixed$pval
names(ewas_CRP_fixed)[names(ewas_CRP_fixed) == "coef"] <- "beta_coef"


###Doing the high dimensional mediation using medScan######
obj = medScan(z.alpha = z.alpha, z.beta = z.beta, method = "HDMT")
objpvalues=obj$pvalues
pvalues_crpadj=as.data.frame(objpvalues)
pvalues_crpadj$q=p.adjust(pvalues_crpadj$objpvalues, method='fdr')

######Merging pvalue file with CPG file################
fulldata_crpadj=cbind(ewas_stress3_fixed,ewas_CRP_fixed,pvalues_crpadj)
fulldata_crpadj=fulldata_crpadj[,-7]

#Getting significant cpg sites from mediation analysis 
significant_crpadj=subset(fulldata_crpadj, q<0.05)
sigcpg_crpadj=significant_crpadj[,1]


#Extracting CpGs that were FDR significant in Model 1
crp_cpgadj=subset(fulldata_crpadj, cpg == "cg04803208" | cpg == "cg00420390" | cpg == "cg25607249" | cpg == "cg04051649" | cpg == "cg15781610" | cpg == "cg09276842" | cpg == "cg02508743" | cpg == "cg24837149"
                  | cpg == "cg03699074" | cpg == "cg26010590" | cpg == "cg11956636" | cpg == "cg21327712" | cpg == "cg01406381")


