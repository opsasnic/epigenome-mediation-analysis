The following code corresponds to our manuscript entitled: "Epigenome-Wide Mediation Analysis of the Relationship between Psychosocial Stress and Cardiometabolic Risk Factors in the Health and Retirement Study (HRS)". 
This document will detail each of the R files that have been included: 
1) Demographic Tables - This code calcluates the sample demographic characteristics (Table 1). It also characterizes the difference in demographics between those included vs. excluded in the final analysis from the entire DNA methylation sample (n=4018) (Table S1). 
2) Total effects_Stress CVD risk factors - This code calculates the total effects of psychosocial stress on each of the 10 cardiometabolic risk factors in Model 1 and Model 2
3) "Outcome"_exposure mediator - There are 4 separate files for BMI, WC, HDL-C and CRP. This code runs the first model needed for mediation analysis, performing epigenome wide associations between stress and DNA methylation. We have four different sets of code because each outcome has a slightly different sample size, which needs to be accounted for in mediation. 
4) "Outcome"_exposure mediator_smoking - There are 4 separate files for BMI, WC, HDL-C and CRP. This code runs the first model needed for the sensitivity mediation analysis (adjusting for smoking), performing epigenome wide associations between stress and DNA methylation.
5) "Outcome"_mediator outcome -
6) "Outcome"_mediator outcome_smoking -
