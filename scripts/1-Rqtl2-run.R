# Script created by Rebecca Gould, UGA


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



####################################################
## Read in the control file (gm.json)
####################################################

  control <- read_cross2(file = "~/Rqtl2-Kidney-Glutathione/data/control.json")

  

####################################################
## Genotype probabilities and allele probabilities - calculated by the Jackson Laboratory
####################################################
  
  probs <- readRDS("~/Rqtl2-Kidney-Glutathione/data/Pazdro_GigaMUGA_genoprobs_qced_8state_sorted.rds")


  
####################################################
## Variant files
####################################################

  # set the download timeout to 900 seconds from the default 60 seconds to help with large file downloads
  options(timeout=900) 
  # download the data files needed for SNP association mapping and obtaining genes in QTL intervals
  download.file(url="https://ndownloader.figshare.com/files/18533342", destfile="./data/cc_variants.sqlite") 
  download.file(url="https://ndownloader.figshare.com/files/24607961", destfile="./data/mouse_genes.sqlite")
  download.file(url="https://ndownloader.figshare.com/files/24607970", destfile="./data/mouse_genes_mgi.sqlite")
  
  #Will need these for the final lesson episodes on SNP association mapping and QTL analysis in Diversity Outbred mice. Make sure they are the most updated versions!
  query_variants <- create_variant_query_func("./data/cc_variants.sqlite")
  query_genes_mgi <- create_gene_query_func("./data/mouse_genes_mgi.sqlite")
  query_genes <- create_gene_query_func("./data/mouse_genes.sqlite")

  

####################################################
## Calculating kinship
####################################################

  kinship_loco <- calc_kinship(probs = probs, "loco", use_allele_probs = TRUE, cores = 10)

  

####################################################
## Read in pheno file and adjust for R/qtl2
####################################################

  #to have the phenotype file for reference - can be used when plotting the data to see if it needs to be transformed
  pheno <- read.csv(file = "./data/pheno_covar.csv", header = TRUE)
  
  #make row names the ID of each sample
  rownames(pheno) <- pheno$id
  
  #change sex to numeric variables
  pheno$sex[pheno$sex == "M"] <- 1
  pheno$sex[pheno$sex == "F"] <- 0
  
  #both added covariates must be numeric, not characters 
  pheno$sex <- as.numeric(pheno$sex)
  pheno$generation <- as.numeric(pheno$generation)  


  
####################################################
## Transform data
####################################################
  
  # RankZ function
  rankZ <- function(x) {x <- rank(x, na.last = "keep") / (length(x) - sum(is.na(x)) + 1)
  return(qnorm(x))}
  
  # Rank Z transformations
  pheno$zKidneyGSH = rankZ(pheno$Kidney_GSH)
  pheno$zKidneyGSSG = rankZ(pheno$Kidney_GSSG)
  pheno$zKidneyTotalGSH = rankZ(pheno$Kidney_Total_GSH)
  pheno$zKidneyGSH_GSSGRatio = rankZ(pheno$Kidney_GSH_GSSG_Ratio)
  pheno$zKidneyRedoxPotentialGSSG2GSH = rankZ(pheno$Kidney_Redox_Potential_GSSG_2GSH)
  pheno$zBUN = rankZ(pheno$BUN)


  
####################################################
## Add covariates
####################################################

  #adding sex and generation as covariates
  sexgen = model.matrix(~ sex + generation, data = pheno)[,-1]

  

####################################################
## Kidney GSH
## Plot Genome Scans with Permutation Tests
####################################################

  qtlscan_KidneyGSH <- scan1(genoprobs = probs, pheno = pheno["zKidneyGSH"], kinship = kinship_loco, addcovar = sexgen, cores=2)
  perm_KidneyGSH <- scan1perm(genoprobs = probs, pheno = pheno["zKidneyGSH"], addcovar = sexgen, n_perm = 1000, cores=10)
  
  Xcovar = get_x_covar(control)
  perm_strata = mat2strata(Xcovar)
  perm_X_KidneyGSH <- scan1perm(genoprobs = probs, pheno = pheno["zKidneyGSH"], addcovar = sexgen,n_perm = 1000, perm_Xsp = TRUE, perm_strata = perm_strata, chr_lengths = chr_lengths(control$gmap), cores=10)
  
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_KidneyGSH = summary(perm_KidneyGSH, alpha = c(0.2, 0.1, 0.05))
  threshold_X_KidneyGSH = summary(perm_X_KidneyGSH, alpha = c(0.2, 0.1, 0.05))
  
  plot_scan1(x = qtlscan_KidneyGSH, map = control$gmap,  main = "Genome Scan for Kidney GSH (Autosome vs X)", ylim = c(0,11))
  segments(x0 = 0, y0 = threshold_X_KidneyGSH$A, x1 = 1695, y1 =   threshold_X_KidneyGSH$A, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
  segments(x0 = 1695, y0 = threshold_X_KidneyGSH$X, x1 = 2000, y1 = threshold_X_KidneyGSH$X, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
  
  #using gmap (cM)
  gmap_peaksGSH <- find_peaks(scan1_output = qtlscan_KidneyGSH, map = control$gmap, threshold = summary(perm_KidneyGSH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  
  #using pmap (Mbp)
  peaksGSH <- find_peaks(scan1_output = qtlscan_KidneyGSH, map = control$pmap, threshold = summary(perm_KidneyGSH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  
# Kidney GSH --- Chromosome X
  par(mar=c(4.1, 4.1, 2.6, 2.6))

  #using gmap (cM)
  chr = "X"
  coef_blup_KidneyGSH_chrX <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_KidneyGSH_chrX, map = control$gmap, scan1_output = qtlscan_KidneyGSH, main = "Kidney GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(15,40)
  plot_coefCC(x = coef_blup_KidneyGSH_chrX, map = control$gmap, scan1_output = qtlscan_KidneyGSH, main = "Kidney GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = "X"
  variants_KidneyGSH_chrX <- query_variants(chr, 48.2, 52.9)
  out_snps_KidneyGSH_chrX <- scan1snps(genoprobs = probs, map = control$pmap, pheno = pheno["zKidneyGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                       chr = chr, start = 48.2, end = 52.9, keep_all_snps = TRUE)
  plot_snpasso(out_snps_KidneyGSH_chrX$lod, out_snps_KidneyGSH_chrX$snpinfo, main = "Kidney GSH SNPs")
  
  KidneyGSH_Genes_MGI_chrX <- query_genes_mgi(chr = chr, start = 48.2, end = 52.9)
  plot(out_snps_KidneyGSH_chrX$lod, out_snps_KidneyGSH_chrX$snpinfo, drop_hilit=1.5, genes = KidneyGSH_Genes_MGI_chrX, main = "Kidney GSH Genes MGI")

# Kidney GSH --- Chromosome 11
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = 11
  coef_blup_KidneyGSH_chr11 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_KidneyGSH_chr11, map = control$gmap, scan1_output = qtlscan_KidneyGSH, main = "Kidney GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(45,75)
  plot_coefCC(x = coef_blup_KidneyGSH_chr11, map = control$gmap, scan1_output = qtlscan_KidneyGSH, main = "Kidney GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = 11
  variants_KidneyGSH_chr11 <- query_variants(chr, 99, 102)
  out_snps_KidneyGSH_chr11 <- scan1snps(genoprobs = probs, map = control$pmap, pheno = pheno["zKidneyGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                        chr = chr, start = 99, end = 102, keep_all_snps = TRUE)
  plot_snpasso(out_snps_KidneyGSH_chr11$lod, out_snps_KidneyGSH_chr11$snpinfo, main = "Kidney GSH SNPs")
  
  KidneyGSH_Genes_MGI_chr11 <- query_genes_mgi(chr = chr, start = 99, end = 102)
  plot(out_snps_KidneyGSH_chr11$lod, out_snps_KidneyGSH_chr11$snpinfo, drop_hilit=1.5, genes = KidneyGSH_Genes_MGI_chr11, main = "Kidney GSH Genes MGI")
  
       
       
####################################################
## Kidney GSSG
## Plot Genome Scans with Permutation Tests
####################################################
  
  qtlscan_KidneyGSSG <- scan1(genoprobs = probs, pheno = pheno["zKidneyGSSG"], kinship = kinship_loco, addcovar = sexgen, cores=2)
  perm_KidneyGSSG <- scan1perm(genoprobs = probs, pheno = pheno["zKidneyGSSG"], addcovar = sexgen, n_perm = 1000, cores=10)
  
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_KidneyGSSG = summary(perm_KidneyGSSG, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_KidneyGSSG, map = control$gmap,  main = "Genome Scan for Kidney GSSG", ylim = c(0,11))
  abline(h = threshold_KidneyGSSG, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
  
  #using gmap (cM)
  gmap_peaksGSSG <- find_peaks(scan1_output = qtlscan_KidneyGSSG, map = control$gmap, threshold = summary(perm_KidneyGSSG, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  
  #using pmap (Mbp)
  peaksGSSG <- find_peaks(scan1_output = qtlscan_KidneyGSSG, map = control$pmap, threshold = summary(perm_KidneyGSSG, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  
  
  
####################################################
## Kidney Total Glutathione
## Plot Genome Scans with Permutation Tests
####################################################
  
  qtlscan_KidneyTotalGSH<- scan1(genoprobs = probs, pheno = pheno["zKidneyTotalGSH"], kinship = kinship_loco, addcovar = sexgen, cores=2)
  perm_KidneyTotalGSH <- scan1perm(genoprobs = probs, pheno = pheno["zKidneyTotalGSH"], addcovar = sexgen, n_perm = 1000, cores=10)
  
  Xcovar = get_x_covar(control)
  perm_strata = mat2strata(Xcovar)
  perm_X_KidneyTotalGSH <- scan1perm(genoprobs = probs, pheno = pheno["zKidneyTotalGSH"], addcovar = sexgen, n_perm = 1000, perm_Xsp = TRUE, perm_strata = perm_strata, chr_lengths = chr_lengths(control$gmap), cores=10)
  
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_KidneyTotalGSH = summary(perm_KidneyTotalGSH, alpha = c(0.2, 0.1, 0.05))
  threshold_X_KidneyTotalGSH = summary(perm_X_KidneyTotalGSH, alpha = c(0.2, 0.1, 0.05))

  plot_scan1(x = qtlscan_KidneyTotalGSH, map = control$gmap,  main = "Genome Scan for Kidney Total GSH (Autosome vs X)", ylim = c(0,11))
  segments(x0 = 0, y0 = threshold_X_KidneyTotalGSH$A, x1 = 1695, y1 =   threshold_X_KidneyTotalGSH$A, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
  segments(x0 = 1695, y0 = threshold_X_KidneyTotalGSH$X, x1 = 2000, y1 = threshold_X_KidneyTotalGSH$X, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
  
  #using gmap (cM)
  gmap_peaksTotalGSH <- find_peaks(scan1_output = qtlscan_KidneyTotalGSH, map = control$gmap, threshold = summary(perm_KidneyTotalGSH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  
  #using pmap (Mbp)
  peaksTotalGSH <- find_peaks(scan1_output = qtlscan_KidneyTotalGSH, map = control$pmap, threshold = summary(perm_KidneyTotalGSH, alpha = 0.2), peakdrop = 1.0, prob = 0.95, expand2markers = FALSE)
  
  
# Kidney Total Glutathione --- Chromosome X 
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = "X"
  coef_blup_KidneyTotalGSH_chrX <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyTotalGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_KidneyTotalGSH_chrX, map = control$gmap, scan1_output = qtlscan_KidneyTotalGSH, main = "Kidney Total GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(15,40)
  plot_coefCC(x = coef_blup_KidneyTotalGSH_chrX, map = control$gmap, scan1_output = qtlscan_KidneyTotalGSH, main = "Kidney Total GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = "X"
  variants_KidneyTotalGSH_chrX <- query_variants(chr, 48.2, 52.9)
  out_snps_KidneyTotalGSH_chrX <- scan1snps(genoprobs = probs, map = control$pmap, pheno = pheno["zKidneyTotalGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                            chr = chr, start = 48.2, end = 52.9, keep_all_snps = TRUE)
  plot_snpasso(out_snps_KidneyTotalGSH_chrX$lod, out_snps_KidneyTotalGSH_chrX$snpinfo, main = "Kidney Total GSH SNPs")
  
  KidneyTotalGSH_Genes_MGI_chrX <- query_genes_mgi(chr = chr, start = 48.2, end = 52.9)
  plot(out_snps_KidneyTotalGSH_chrX$lod, out_snps_KidneyTotalGSH_chrX$snpinfo, drop_hilit=1.5, genes = KidneyTotalGSH_Genes_MGI_chrX, main = "Kidney Total GSH Genes MGI")

# Kidney GSH --- Chromosome 11
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = 11
  coef_blup_KidneyTotalGSH_chr11 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyTotalGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_KidneyTotalGSH_chr11, map = control$gmap, scan1_output = qtlscan_KidneyTotalGSH, main = "Kidney Total GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(45,75)
  plot_coefCC(x = coef_blup_KidneyTotalGSH_chr11, map = control$gmap, scan1_output = qtlscan_KidneyTotalGSH, main = "Kidney Total GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = 11
  variants_KidneyTotalGSH_chr11 <- query_variants(chr, 99, 102)
  out_snps_KidneyTotalGSH_chr11 <- scan1snps(genoprobs = probs, map = control$pmap, pheno = pheno["zKidneyTotalGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                             chr = chr, start = 99, end = 102, keep_all_snps = TRUE)
  plot_snpasso(out_snps_KidneyTotalGSH_chr11$lod, out_snps_KidneyTotalGSH_chr11$snpinfo, main = "Kidney Total GSH SNPs")
  
  KidneyTotalGSH_Genes_MGI_chr11 <- query_genes_mgi(chr = chr, start = 99, end = 102)
  plot(out_snps_KidneyTotalGSH_chr11$lod, out_snps_KidneyTotalGSH_chr11$snpinfo, drop_hilit=1.5, genes = KidneyGSH_Genes_MGI_chr11, main = "Kidney Total GSH Genes MGI")
  
  
  
####################################################
## Kidney GSH/GSSG Ratio
## Plot Genome Scans with Permutation Tests
####################################################
  
  qtlscan_KidneyGSH_GSSGRatio<- scan1(genoprobs = probs, pheno = pheno["zKidneyGSH_GSSGRatio"], kinship = kinship_loco, addcovar = sexgen, cores=2)
  perm_KidneyGSH_GSSGRatio <- scan1perm(genoprobs = probs, pheno = pheno["zKidneyGSH_GSSGRatio"], addcovar = sexgen, n_perm = 1000, cores=10)
  
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_KidneyGSH_GSSGRatio = summary(perm_KidneyGSH_GSSGRatio, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_KidneyGSH_GSSGRatio, map = control$gmap,  main = "Genome Scan for Kidney GSH/GSSG Ratio", ylim = c(0,11))
  abline(h = threshold_KidneyGSH_GSSGRatio, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
  
  #using gmap (cM)
  gmap_peaksGSH_GSSGRatio <- find_peaks(scan1_output = qtlscan_KidneyGSH_GSSGRatio, map = control$gmap, threshold = summary(perm_KidneyGSH_GSSGRatio, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  
  #using pmap (Mbp)
  peaksGSH_GSSGRatio <- find_peaks(scan1_output = qtlscan_KidneyGSH_GSSGRatio, map = control$pmap, threshold = summary(perm_KidneyGSH_GSSGRatio, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  

  
  
####################################################
## Kidney Glutathione Redox Potential
## Plot Genome Scans with Permutation Tests
####################################################
  
  qtlscan_KidneyRedoxPotentialGSSG2GSH<- scan1(genoprobs = probs, pheno = pheno["zKidneyRedoxPotentialGSSG2GSH"], kinship = kinship_loco, addcovar = sexgen, cores=2)
  perm_KidneyRedoxPotentialGSSG2GSH <- scan1perm(genoprobs = probs, pheno = pheno["zKidneyRedoxPotentialGSSG2GSH"], addcovar = sexgen, n_perm = 1000, cores=10)
  
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_KidneyRedoxPotentialGSSG2GSH = summary(perm_KidneyRedoxPotentialGSSG2GSH, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_KidneyRedoxPotentialGSSG2GSH, map = control$gmap,  main = "Genome Scan for Kidney Redox Potential GSSG/2GSH", ylim = c(0,11))
  abline(h = threshold_KidneyRedoxPotentialGSSG2GSH, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
  
  #using gmap (cM)
  gmap_peaksKidneyRedoxPotentialGSSG2GSH <- find_peaks(scan1_output = qtlscan_KidneyRedoxPotentialGSSG2GSH, map = control$gmap, threshold = summary(perm_KidneyRedoxPotentialGSSG2GSH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  
  #using pmap (Mbp)
  peaksKidneyRedoxPotentialGSSG2GSH <- find_peaks(scan1_output = qtlscan_KidneyRedoxPotentialGSSG2GSH, map = control$pmap, threshold = summary(perm_KidneyRedoxPotentialGSSG2GSH, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  
# Kidney Glutathione Redox Potential --- Chromosome 14
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  #using gmap (cM)
  chr = 14
  coef_blup_KidneyRedoxPotential_chr14 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zKidneyRedoxPotentialGSSG2GSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
  plot_coefCC(x = coef_blup_KidneyRedoxPotential_chr14, map = control$gmap, scan1_output = qtlscan_KidneyRedoxPotentialGSSG2GSH, main = "Kidney Eh BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(0,25)
  plot_coefCC(x = coef_blup_KidneyRedoxPotential_chr14, map = control$gmap, scan1_output = qtlscan_KidneyRedoxPotentialGSSG2GSH, main = "Kidney Eh BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
  #using pmap (Mbp)
  chr = 14
  variants_KidneyRedoxPotential_chr14 <- query_variants(chr, 21, 25)
  out_snps_KidneyRedoxPotential_chr14 <- scan1snps(genoprobs = probs, map = control$pmap, pheno = pheno["zKidneyRedoxPotentialGSSG2GSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                                   chr = chr, start = 21, end = 25, keep_all_snps = TRUE)
  plot_snpasso(out_snps_KidneyRedoxPotential_chr14$lod, out_snps_KidneyRedoxPotential_chr14$snpinfo, main = "Kidney Eh SNPs")
  KidneyEh_Genes_MGI_chr14 <- query_genes_mgi(chr = chr, start = 21, end = 25)
  plot(out_snps_KidneyRedoxPotential_chr14$lod, out_snps_KidneyRedoxPotential_chr14$snpinfo, drop_hilit=1.5, genes = KidneyEh_Genes_MGI_chr14, main = "Kidney Eh Genes MGI")
  


####################################################
## BUN
## Plot Genome Scans with Permutation Tests
####################################################
  
  qtlscan_BUN <- scan1(genoprobs = probs, pheno = pheno["zBUN"], kinship = kinship_loco, addcovar = sexgen, cores=2)
  perm_BUN <- scan1perm(genoprobs = probs, pheno = pheno["zBUN"], addcovar = sexgen, n_perm = 1000, cores=10)
  
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_BUN = summary(perm_BUN, alpha = c(0.2, 0.1, 0.05))
  plot_scan1(x = qtlscan_BUN, map = control$gmap,  main = "Genome Scan for BUN", ylim = c(0,11))
  abline(h = threshold_BUN, col = c("purple", "red", "blue"), lwd = 2, lty = "dashed")
  
  #using gmap (cM)
  gmap_peaksBUN <- find_peaks(scan1_output = qtlscan_BUN, map = control$gmap, threshold = summary(perm_BUN, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  
  #using pmap (Mbp)
  peaksBUN <- find_peaks(scan1_output = qtlscan_BUN, map = control$pmap, threshold = summary(perm_BUN, alpha = 0.2), peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  
  
  
####################################################
## Export all QTL with LOD scores > 6 and all genes in QTL intervals
####################################################

qtl_gmap <- find_peaks(scans, map = control$gmap, threshold = 6, peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
qtl_gmap

write_xlsx(list("GSH chrX" = KidneyGSH_Genes_MGI_chrX, 
                "TotalGSH chrX" = KidneyTotalGSH_Genes_MGI_chrX,
                "GSH chr11" = KidneyGSH_Genes_MGI_chr11,
                "Total GSH chr11" = KidneyTotalGSH_Genes_MGI_chr11,
                "Eh chr14" = KidneyEh_Genes_MGI_chr14),
           "GlutathioneGenesMGI-RankZ-sexgen.xlsx")



####################################################
## Export all QTL with LOD scores > 6 and all genes in QTL intervals
####################################################

qtl_gmap <- find_peaks(scans, map = control$gmap, threshold = 6, peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
qtl_pmap <- find_peaks(scans, map = control$pmap, threshold = 6, peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)

qtl_gmap$marker.id <- find_marker(map = control$gmap, chr = qtl_gmap$chr, pos = qtl_gmap$pos)
qtl_pmap$marker.id <- find_marker(map = control$pmap, chr = qtl_pmap$chr, pos = qtl_pmap$pos)


write_xlsx(list("QTL List RankZ SexGen - cM" = qtl_gmap,
                "QTL List RankZ SexGen - Mbp" = qtl_pmap),
           "QTL List.xlsx")




