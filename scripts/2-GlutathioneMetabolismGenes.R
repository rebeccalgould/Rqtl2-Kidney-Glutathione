# Script created by Rebecca Gould, UGA
# Using Rdata environment created by 1-Rqtl2-run.R script


####################################################
## Load in necessary packages
####################################################

library(qtl2)
library (tidyverse)
library (readxl)
library(yaml)
library(devtools)
library(jsonlite)
library (data.table)
library (RcppEigen)
library (pander)
library (writexl)
library (RSQLite)


##################################################################
## Evaluate gluthatione synthesis and recycling genes
##################################################################

## Kidney GSH ##
  
  #Glutathione Peroxidase (Gpx1) - Chr 9 59.24 cM
  chr = 9
  coef_blup_KidneyGSH_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_KidneyGSH_chr9, map = control$gmap, scan1_output = qtlscan_KidneyGSH, main = "Kidney GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(58.5,60)
  plot_coefCC(x = coef_blup_KidneyGSH_chr9, map = control$gmap, scan1_output = qtlscan_KidneyGSH, main = "Gpx1 Position -- Kidney GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutathione cysteine ligase - catalytic subunit (Gclc) - Chr 9 43.36 cM
  chr = 9
  #coef_blup_KidneyGSH_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  #plot_coefCC(x = coef_blup_KidneyGSH_chr9, map = control$gmap, scan1_output = qtlscan_KidneyGSH, main = "Kidney GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(42.5,44)
  plot_coefCC(x = coef_blup_KidneyGSH_chr9, map = control$gmap, scan1_output = qtlscan_KidneyGSH, main = "Gclc Position -- Kidney GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutamate-cysteine ligase – modifier subunit (Gclm) - Chr 3 52.94 cM
  chr = 3
  coef_blup_KidneyGSH_chr3 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_KidneyGSH_chr3, map = control$gmap, scan1_output = qtlscan_KidneyGSH, main = "Kidney GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(51.5,53.5)
  plot_coefCC(x = coef_blup_KidneyGSH_chr3, map = control$gmap, scan1_output = qtlscan_KidneyGSH, main = "Gclm Position -- Kidney GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutathione synthetase (Gss) - Chr 2 77.26 cM
  chr = 2
  coef_blup_KidneyGSH_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_KidneyGSH_chr2, map = control$gmap, scan1_output = qtlscan_KidneyGSH, main = "Kidney GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(76.5,78)
  plot_coefCC(x = coef_blup_KidneyGSH_chr2, map = control$gmap, scan1_output = qtlscan_KidneyGSH, main = "Gss Position -- Kidney GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutathione reductase  (Gsr) - Chr 8 20.69 cM
  chr = 8
  coef_blup_KidneyGSH_chr8 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_KidneyGSH_chr8, map = control$gmap, scan1_output = qtlscan_KidneyGSH, main = "Kidney GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(19,22.5)
  plot_coefCC(x = coef_blup_KidneyGSH_chr8, map = control$gmap, scan1_output = qtlscan_KidneyGSH, main = "Gsr Position -- Kidney GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)


## Kidney GSSG ##
  
  #Glutathione Peroxidase (Gpx1) - Chr 9 59.24 cM
  chr = 9
  coef_blup_KidneyGSSG_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_KidneyGSSG_chr9, map = control$gmap, scan1_output = qtlscan_KidneyGSSG, main = "Kidney GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(58.5,60)
  plot_coefCC(x = coef_blup_KidneyGSSG_chr9, map = control$gmap, scan1_output = qtlscan_KidneyGSSG, main = "Gpx1 Position -- Kidney GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutathione cysteine ligase - catalytic subunit (Gclc) - Chr 9 43.36 cM
  chr = 9
  #coef_blup_KidneyGSSG_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  #plot_coefCC(x = coef_blup_KidneyGSSG_chr9, map = control$gmap, scan1_output = qtlscan_KidneyGSSG, main = "Kidney GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(42.5,44)
  plot_coefCC(x = coef_blup_KidneyGSSG_chr9, map = control$gmap, scan1_output = qtlscan_KidneyGSSG, main = "Gclc Position -- Kidney GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutamate-cysteine ligase – modifier subunit (Gclm) - Chr 3 52.94 cM
  chr = 3
  coef_blup_KidneyGSSG_chr3 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_KidneyGSSG_chr3, map = control$gmap, scan1_output = qtlscan_KidneyGSSG, main = "Kidney GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(51.5,53.5)
  plot_coefCC(x = coef_blup_KidneyGSSG_chr3, map = control$gmap, scan1_output = qtlscan_KidneyGSSG, main = "Gclm Position -- Kidney GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutathione synthetase (Gss) - Chr 2 77.26 cM
  chr = 2
  coef_blup_KidneyGSSG_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_KidneyGSSG_chr2, map = control$gmap, scan1_output = qtlscan_KidneyGSSG, main = "Kidney GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(76.5,78)
  plot_coefCC(x = coef_blup_KidneyGSSG_chr2, map = control$gmap, scan1_output = qtlscan_KidneyGSSG, main = "Gss Position -- Kidney GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutathione reductase  (Gsr) - Chr 8 20.69 cM
  chr = 8
  coef_blup_KidneyGSSG_chr8 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyGSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_KidneyGSSG_chr8, map = control$gmap, scan1_output = qtlscan_KidneyGSSG, main = "Kidney GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(19,22.5)
  plot_coefCC(x = coef_blup_KidneyGSSG_chr8, map = control$gmap, scan1_output = qtlscan_KidneyGSSG, main = "Gsr Position -- Kidney GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  

## Kidney Total Glutathione ##
  
  #Glutathione Peroxidase (Gpx1) - Chr 9 59.24 cM
  chr = 9
  coef_blup_KidneyTotalGSH_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyTotalGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_KidneyTotalGSH_chr9, map = control$gmap, scan1_output = qtlscan_KidneyTotalGSH, main = "Kidney Total GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(58.5,60)
  plot_coefCC(x = coef_blup_KidneyTotalGSH_chr9, map = control$gmap, scan1_output = qtlscan_KidneyTotalGSH, main = "Gpx1 Position -- Kidney Total GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutathione cysteine ligase - catalytic subunit (Gclc) - Chr 9 43.36 cM
  chr = 9
  #coef_blup_KidneyTotalGSH_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyTotalGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  #plot_coefCC(x = coef_blup_KidneyTotalGSH_chr9, map = control$gmap, scan1_output = qtlscan_KidneyTotalGSH, main = "Kidney Total GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(42.5,44)
  plot_coefCC(x = coef_blup_KidneyTotalGSH_chr9, map = control$gmap, scan1_output = qtlscan_KidneyTotalGSH, main = "Gclc Position -- Kidney Total GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutamate-cysteine ligase – modifier subunit (Gclm) - Chr 3 52.94 cM
  chr = 3
  coef_blup_KidneyTotalGSH_chr3 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyTotalGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_KidneyTotalGSH_chr3, map = control$gmap, scan1_output = qtlscan_KidneyTotalGSH, main = "Kidney Total GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(51.5,53.5)
  plot_coefCC(x = coef_blup_KidneyTotalGSH_chr3, map = control$gmap, scan1_output = qtlscan_KidneyTotalGSH, main = "Gclm Position -- Kidney Total GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutathione synthetase (Gss) - Chr 2 77.26 cM
  chr = 2
  coef_blup_KidneyTotalGSH_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyTotalGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_KidneyTotalGSH_chr2, map = control$gmap, scan1_output = qtlscan_KidneyTotalGSH, main = "Kidney Total GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(76.5,78)
  plot_coefCC(x = coef_blup_KidneyTotalGSH_chr2, map = control$gmap, scan1_output = qtlscan_KidneyTotalGSH, main = "Gss Position -- Kidney Total GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutathione reductase  (Gsr) - Chr 8 20.69 cM
  chr = 8
  coef_blup_KidneyTotalGSH_chr8 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyTotalGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_KidneyTotalGSH_chr8, map = control$gmap, scan1_output = qtlscan_KidneyTotalGSH, main = "Kidney Total GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(19,22.5)
  plot_coefCC(x = coef_blup_KidneyTotalGSH_chr8, map = control$gmap, scan1_output = qtlscan_KidneyTotalGSH, main = "Gsr Position -- Kidney Total GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  

## Kidney GSH/GSSG Ratio ##
  
  #Glutathione Peroxidase (Gpx1) - Chr 9 59.24 cM
  chr = 9
  coef_blup_KidneyGSH_GSSGRatio_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyGSH_GSSGRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_KidneyGSH_GSSGRatio_chr9, map = control$gmap, scan1_output = qtlscan_KidneyGSH_GSSGRatio, main = "Kidney GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(58.5,60)
  plot_coefCC(x = coef_blup_KidneyGSH_GSSGRatio_chr9, map = control$gmap, scan1_output = qtlscan_KidneyGSH_GSSGRatio, main = "Gpx1 Position -- Kidney GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutathione cysteine ligase - catalytic subunit (Gclc) - Chr 9 43.36 cM
  chr = 9
  #coef_blup_KidneyGSH_GSSGRatio_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyGSH_GSSGRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  #plot_coefCC(x = coef_blup_KidneyGSH_GSSGRatio_chr9, map = control$gmap, scan1_output = qtlscan_KidneyGSH_GSSGRatio, main = "Kidney GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(42.5,44)
  plot_coefCC(x = coef_blup_KidneyGSH_GSSGRatio_chr9, map = control$gmap, scan1_output = qtlscan_KidneyGSH_GSSGRatio, main = "Gclc Position -- Kidney GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutamate-cysteine ligase – modifier subunit (Gclm) - Chr 3 52.94 cM
  chr = 3
  coef_blup_KidneyGSH_GSSGRatio_chr3 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyGSH_GSSGRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_KidneyGSH_GSSGRatio_chr3, map = control$gmap, scan1_output = qtlscan_KidneyGSH_GSSGRatio, main = "Kidney GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(51.5,53.5)
  plot_coefCC(x = coef_blup_KidneyGSH_GSSGRatio_chr3, map = control$gmap, scan1_output = qtlscan_KidneyGSH_GSSGRatio, main = "Gclm Position -- Kidney GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutathione synthetase (Gss) - Chr 2 77.26 cM
  chr = 2
  coef_blup_KidneyGSH_GSSGRatio_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyGSH_GSSGRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_KidneyGSH_GSSGRatio_chr2, map = control$gmap, scan1_output = qtlscan_KidneyGSH_GSSGRatio, main = "Kidney GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(76.5,78)
  plot_coefCC(x = coef_blup_KidneyGSH_GSSGRatio_chr2, map = control$gmap, scan1_output = qtlscan_KidneyGSH_GSSGRatio, main = "Gss Position -- Kidney GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutathione reductase  (Gsr) - Chr 8 20.69 cM
  chr = 8
  coef_blup_KidneyGSH_GSSGRatio_chr8 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyGSH_GSSGRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_KidneyGSH_GSSGRatio_chr8, map = control$gmap, scan1_output = qtlscan_KidneyGSH_GSSGRatio, main = "Kidney GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(19,22.5)
  plot_coefCC(x = coef_blup_KidneyGSH_GSSGRatio_chr8, map = control$gmap, scan1_output = qtlscan_KidneyGSH_GSSGRatio, main = "Gsr Position -- Kidney GSH/GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  
## Kidney Glutathione Redox Potential ##
  
  #Glutathione Peroxidase (Gpx1) - Chr 9 59.24 cM
  chr = 9
  coef_blup_KidneyRedoxPotentialGSSG2GSH_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyRedoxPotentialGSSG2GSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_KidneyRedoxPotentialGSSG2GSH_chr9, map = control$gmap, scan1_output = qtlscan_KidneyRedoxPotentialGSSG2GSH, main = "Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(58.5,60)
  plot_coefCC(x = coef_blup_KidneyRedoxPotentialGSSG2GSH_chr9, map = control$gmap, scan1_output = qtlscan_KidneyRedoxPotentialGSSG2GSH, main = "Gpx1 Position -- Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutathione cysteine ligase - catalytic subunit (Gclc) - Chr 9 43.36 cM
  chr = 9
  #coef_blup_KidneyRedoxPotentialGSSG2GSH_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyRedoxPotentialGSSG2GSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  #plot_coefCC(x = coef_blup_KidneyRedoxPotentialGSSG2GSH_chr9, map = control$gmap, scan1_output = qtlscan_KidneyRedoxPotentialGSSG2GSH, main = "Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(42.5,44)
  plot_coefCC(x = coef_blup_KidneyRedoxPotentialGSSG2GSH_chr9, map = control$gmap, scan1_output = qtlscan_KidneyRedoxPotentialGSSG2GSH, main = "Gclc Position -- Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutamate-cysteine ligase – modifier subunit (Gclm) - Chr 3 52.94 cM
  chr = 3
  coef_blup_KidneyRedoxPotentialGSSG2GSH_chr3 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyRedoxPotentialGSSG2GSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_KidneyRedoxPotentialGSSG2GSH_chr3, map = control$gmap, scan1_output = qtlscan_KidneyRedoxPotentialGSSG2GSH, main = "Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(51.5,53.5)
  plot_coefCC(x = coef_blup_KidneyRedoxPotentialGSSG2GSH_chr3, map = control$gmap, scan1_output = qtlscan_KidneyRedoxPotentialGSSG2GSH, main = "Gclm Position -- Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutathione synthetase (Gss) - Chr 2 77.26 cM
  chr = 2
  coef_blup_KidneyRedoxPotentialGSSG2GSH_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyRedoxPotentialGSSG2GSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_KidneyRedoxPotentialGSSG2GSH_chr2, map = control$gmap, scan1_output = qtlscan_KidneyRedoxPotentialGSSG2GSH, main = "Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(76.5,78)
  plot_coefCC(x = coef_blup_KidneyRedoxPotentialGSSG2GSH_chr2, map = control$gmap, scan1_output = qtlscan_KidneyRedoxPotentialGSSG2GSH, main = "Gss Position -- Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #Glutathione reductase  (Gsr) - Chr 8 20.69 cM
  chr = 8
  coef_blup_KidneyRedoxPotentialGSSG2GSH_chr8 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyRedoxPotentialGSSG2GSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_KidneyRedoxPotentialGSSG2GSH_chr8, map = control$gmap, scan1_output = qtlscan_KidneyRedoxPotentialGSSG2GSH, main = "Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(19,22.5)
  plot_coefCC(x = coef_blup_KidneyRedoxPotentialGSSG2GSH_chr8, map = control$gmap, scan1_output = qtlscan_KidneyRedoxPotentialGSSG2GSH, main = "Gsr Position -- Kidney Redox Potential GSSG/2GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  
  
  
  