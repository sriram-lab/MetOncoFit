#### Lyssiotis Imports ####

library(rio)
# import BxPC3 data as bx
bx <- import("~jennadiegel/Documents/SURE/DES/Lyssiotis_data/CCLs_BxPC3.csv")
# import Capan1 data as cp
cp <- import("~jennadiegel/Documents/SURE/DES/Lyssiotis_data/CCLs_Capan1.csv")
# import PANC1 data as pn
pn <- import("~jennadiegel/Documents/SURE/DES/Lyssiotis_data/CCLs_PANC1.csv")
# import TU8902 data as tu
tu <- import("~jennadiegel/Documents/SURE/DES/Lyssiotis_data/CCLs_TU8902.csv")
# import TU8988T data as tut
tut <- import("~jennadiegel/Documents/SURE/DES/Lyssiotis_data/CCLs_TU8988T.csv")
# import UM2 data as um2
um2 <- import("~jennadiegel/Documents/SURE/DES/Lyssiotis_data/CCLs_UM2.csv")
# import UM90 data as um90
um90 <- import("~jennadiegel/Documents/SURE/DES/Lyssiotis_data/CCLs_UM90.csv")

#### Lyssiotis analysis ####

cells <- list(bx, cp, pn, tu, tut, um2, um90)
titles <- list("Log2 distribution of BxPC3","Log2 distribution of Capan1",
               "Log2 distribution of PANC1","Log2 distribution of TU8902",
               "Log2 distribution of TU8988T","Log2 distribution of UM2",
               "Log2 distribution of UM90")
cell_names <- list("BxPC3","Capan1","PANC1","TU8902","TU8988T","UM2","UM90")

#initiate lists to store data 
all_data <- list()
all_fold_changes_log2 <- list()
gene_names_list <- list()
cell_fold_changes <- list()
all_p_values <- list()
z_scores <- list()
p_values <- list()
genes_sig_list <- list()

count <- 1

pdf("~jennadiegel/Documents/SURE/DES/figures/Lyssiotis_DEA_histograms.pdf")

for (cell in cells) {
  #separate out gene names
  genes <- cell[ ,1]
  
  #separate out KO and WT samples
  cell_KO <- cell[ ,2:4]
  cell_WT <- cell[ ,5:7]
  
  #get row-wise means for KO and WT samples
  KO_mean <- rowMeans(cell_KO)  
  WT_mean <- rowMeans(cell_WT)  
  
  #get fold change KO/WT
  fold_change <- KO_mean/WT_mean #getting NaN with WT = 0 - set to max log 2 fold change - P value is arbitrarily large number - add after pvalue calculation
  
  #transform fold change to get normal distribution using log2
  fold_change_log2 <- log2(fold_change)
  
  #separate out inf values 
  bool_inf_pos <- (is.infinite(fold_change_log2) & sign(fold_change_log2) == 1)
  bool_inf_neg <- (is.infinite(fold_change_log2) & sign(fold_change_log2) == -1)
  genes_pos_inf <- genes[bool_inf_pos]
  genes_neg_inf <- genes[bool_inf_neg]
  
  #find NaN and inf values and keep only the non-NaN and non-inf values
  bool <- (is.nan(fold_change_log2)|is.infinite(fold_change_log2))
  fold_change_real <- fold_change[!bool]
  fold_change_log2_real <- fold_change_log2[!bool]
  genes_real <- c(genes[!bool], genes_pos_inf, genes_neg_inf)
  
  #set fold change of inf genes to the max fold change and combine with real values
  fold_change_comb <- c(fold_change_real, rep(max(fold_change_real),length(genes_pos_inf)), rep(min(fold_change_real),length(genes_neg_inf)))
  fold_change_log2_comb <- c(fold_change_log2_real, rep(max(fold_change_log2_real),length(genes_pos_inf)), rep(min(fold_change_log2_real),length(genes_neg_inf)))
  
  ################# need to specify mean and standard deviation
  mean_fc <- mean(fold_change_log2_real)
  stddev <- sd(fold_change_log2_real)
  p_value <- 2*pnorm(fold_change_log2_real, mean_fc, stddev)
  p_value <- c(p_value, rep(10^(-50), (length(genes_pos_inf)+length(genes_neg_inf))))
  
  #find which values are significant and sort out necessary parameters
  bool_p <- (p_value <= 0.05)
  genes_sig <- genes_real[bool_p]
  p_values_keep <- p_value[bool_p]
  fold_change_sig <- fold_change_comb[bool_p]

  #put data frames and vectors into lists for storage
  cell_fold_changes[[count]] <- fold_change_sig
  genes_sig_list[[count]] <- genes_sig
  all_p_values[[count]] <- p_values_keep
  
  #visualize distribution of values as histograms
  hist(fold_change_log2_real, probability = T, breaks = 50,
       main = titles[[count]],
       xlab = 'Log2 of Fold Change')
  
  count <- count + 1}

dev.off()

#exports to use data in matlab or python
export(genes_sig_list, '~jennadiegel/Documents/SURE/DES/Lyssiotis_data/genes_sig.xlsx')
export(all_p_values, '~jennadiegel/Documents/SURE/DES/Lyssiotis_data/p_values.xlsx')
export(cell_fold_changes, '~jennadiegel/Documents/SURE/DES/Lyssiotis_data/fold_change.xlsx')





#### CCLE Imports ####
library(rio)
ccle = import("~jennadiegel/Documents/SURE/DES/data/CCLE_data.csv")

gene_names <- ccle[,2]
wt_vals <- ccle[ ,3:4]
mut_vals <- ccle[ ,-(1:4)]

#### CCLE plot titles ####
titles <- list("Log2 distribution of DAN-G","Log2 distribution of HPAC",
               "Log2 distribution of HUP-T3","Log2 distribution of HUP-T4",
               "Log2 distribution of Hs 766T","Log2 distribution of PA-TU-8902",
               "Log2 distribution of PA-TU-8988S","Log2 distribution of PA-TU-8988T",
               "Log2 distribution of PSN1","Log2 distribution of Panc 02.03",
               "Log2 distribution of Panc 02.13","Log2 distribution of Panc 03.27",
               "Log2 distribution of SNU-213","Log2 distribution of SNU-324",
               "Log2 distribution of SNU-410","Log2 distribution of SU.86.86",
               "Log2 distribution of SW 1990","Log2 distribution of TCC-PAN2",
               "Log2 distribution of AsPC-1","Log2 distribution of BxPC-3",
               "Log2 distribution of CFPAC-1","Log2 distribution of Capan-1",
               "Log2 distribution of Capan-2","Log2 distribution of HPAF-II",
               "Log2 distribution of MIA PaCa-2","Log2 distribution of PANC-1",
               "Log2 distribution of PK-1","Log2 distribution of Panc 04",
               "Log2 distribution of Panc 05", "Log2 distribution of Panc 08",
               "Log2 distribution of Panc 10","Log2 distribution of SUIT-2",
               "Log2 distribution of T3M-4")

pdf("~jennadiegel/Documents/SURE/DES/figures/CCLE_DEA_histograms.pdf")

#### CCLE analysis ####

#fill in NA values with 0 - database excluded 0 values
wt_vals[is.na(wt_vals)] <- 0
mut_vals[is.na(mut_vals)] <- 0

#calculate WT mean from 2 WT cell lines
WT_mean <- (mut_vals[,1] + mut_vals[,2])/2

#determine fold change for all cell lines
fold_changes <- list()
fold_changes_log2 <- list()
noninf_gene_names <- list()
upreg_inf_genes <- list()
downreg_inf_genes <- list()
sig_genes <- list()
p_values <- list()
all_data <- list()
# starting at 3 to exclude first two cell lines that are not normally distributed
for (i in 3:length(mut_vals)) {
  
  
  fold_change <- mut_vals[ ,i]/WT_mean
  fold_change_log2 <- log2(fold_change)
  
  #find NaN and inf values and keep only the non-NaN and non-inf values
  bool <- (is.nan(fold_change_log2)|is.infinite(fold_change_log2))
  fold_change_keep <- fold_change[!bool]
  fold_change_log2_keep <- fold_change_log2[!bool]
  gene_names_keep <- gene_names[!bool]
  
  #separate out inf values, +inf = upregulation, -inf = downregulation
  bool_inf_pos <- (is.infinite(fold_change_log2) & (sign(fold_change_log2) == 1))
  bool_inf_neg <- (is.infinite(fold_change_log2) & (sign(fold_change_log2) == -1))
  gene_names_pos_inf <- gene_names[bool_inf_pos]
  gene_names_neg_inf <- gene_names[bool_inf_neg]
  
  #determine p values for each fold change
  mean_fc <- mean(fold_change_log2_keep)
  stddev <- sd(fold_change_log2_keep)
  p_value <- 2*pnorm(fold_change_log2_keep,mean_fc,stddev)
  p_values[[i]] <- p_value
  
  data_storage <- data.frame(
    gene_name <- c(gene_names_keep, gene_names_pos_inf, gene_names_neg_inf),
    fold_change <- c(fold_change_keep, rep(max(fold_change_keep),length(gene_names_pos_inf)), rep(min(fold_change_keep),length(gene_names_neg_inf))),
    log2_fold_change <- c(fold_change_log2_keep, rep(max(fold_change_log2_keep),length(gene_names_pos_inf)),rep(min(fold_change_log2_keep),length(gene_names_neg_inf))),
    p_value <- c(p_value, rep(10^(-50),(length(gene_names_pos_inf) + length(gene_names_neg_inf)))))
  
  all_data[[i]] <- data_storage
  fold_changes_log2[[i]] <- fold_change_log2_keep
  fold_changes[[i]] <- fold_change_keep
  upreg_inf_genes[[i]] <- gene_names_pos_inf
  downreg_inf_genes[[i]] <- gene_names_neg_inf
  
  sig_genes[[i]] <- c(gene_names_keep, gene_names_pos_inf, gene_names_neg_inf)
  
  hist(fold_change_log2_keep, probability = T, breaks = 50,
       main = titles[[i]],
       xlab = 'Log2 of Fold Change')
  
  i = i + 1}

export(all_data, '~jennadiegel/Documents/SURE/DES/CCLE_data/CCLE.xlsx')
export(fold_changes, '~jennadiegel/Documents/SURE/DES/CCLE_data/fold_changes.xlsx')
export(sig_genes, '~jennadiegel/Documents/SURE/DES/CCLE_data/sig_genes.xlsx')
dev.off()

##figure out why cell line 1 and 2 have weird distributions
##could be modeled with an exponential distribution?
##ignore first two MUT samples



#### TCGA Imports ####
library("readxl")

#pdf to store histograms
pdf("~jennadiegel/Documents/SURE/DES/figures/TCGA_DEA_histograms.pdf")

# make list of sheet names that correspond with sheet names in excel file with organized TCGA data
sheet_names <- list() #fill in file names here - preferably with a loop
wt_count = 1
mut_count = 1
i = 1
for (sheet in 1:38) {
  if (sheet <= 6) {
    name <- paste("wt_TCGA",wt_count, ".txt",sep = "")
    sheet_names[[i]] <- name
    wt_count = wt_count + 1
  }
  else {
    name <- paste("mut_TCGA",mut_count, ".txt", sep = "")
    sheet_names[[i]] <- name
    mut_count = mut_count + 1
  }
  i = i + 1
}


TCGA_wt_data <- list()
TCGA_mut_data <- list()
wt_TPM = matrix()
mut_TPM = matrix()

# put data into lists, separate out TPM values into matrices for wt and mut
wt_count = 1
mut_count = 1
count = 1
for (sheet in sheet_names) {
  if (count <= 6) {
    TCGA_wt_data[[wt_count]] <- read_excel("~jennadiegel/Documents/SURE/DES/TCGA_TPM_data.xlsx", sheet = sheet_names[[count]])
    if (count == 1) {
      wt_TPM = matrix(NA, nrow = dim(TCGA_wt_data[[1]])[1], ncol = 6)
    }
    wt_TPM[ ,wt_count] <- TCGA_wt_data[[wt_count]]$TPMs
    wt_count = wt_count + 1
  }
  else {
    TCGA_mut_data[[mut_count]] <- read_excel("~jennadiegel/Documents/SURE/DES/TCGA_TPM_data.xlsx", sheet = sheet_names[[count]])
    if (count == 7) {
      mut_TPM = matrix(NA, nrow = dim(TCGA_mut_data[[1]])[1], ncol = 32)
    }
    mut_TPM[ ,mut_count] <- TCGA_mut_data[[mut_count]]$TPMs
    mut_count = mut_count + 1
    
  }
  count = count + 1
}

#### TCGA analysis ####

FCS <- list()
log2FCS <- list()
p_vals <- list()
sig_FCS <- list()
genes_sig_list <- list()

wt_mean <- rowMeans(wt_TPM)
count = 1


for(CL in TCGA_mut_data) {
  
  genes = TCGA_mut_data[[count]]$symbol
  FC = mut_TPM[ ,count] / wt_mean
  
  log2FC = log2(FC)
  bool_inf_pos <- (is.infinite(log2FC) & (sign(log2FC) == 1))
  bool_inf_neg <- (is.infinite(log2FC) & (sign(log2FC) == -1))
  
  genes_pos_inf <- genes[bool_inf_pos]
  genes_neg_inf <- genes[bool_inf_neg]
  
  #find NaN and inf values and keep only the non-NaN and non-inf values
  bool <- (is.nan(log2FC)|is.infinite(log2FC))
  FC_real <- FC[!bool]
  log2FC_real <- log2FC[!bool]
  genes_real <- c(genes[!bool], genes_pos_inf, genes_neg_inf)
  
  #set fold change of inf genes to the max fold change and combine with real values
  FC_comb <- c(FC_real, rep(max(FC_real),length(genes_pos_inf)), rep(min(FC_real),length(genes_neg_inf)))
  log2FC_comb <- c(log2FC_real, rep(max(log2FC_real),length(genes_pos_inf)), rep(min(log2FC_real),length(genes_neg_inf)))
  FCS[[count]] <- FC_comb
  log2FCS[[count]] <- log2FC_comb
  
  #need to specify mean and standard deviation
  mean_fc <- mean(log2FC_real)
  stddev <- sd(log2FC_real)
  p_value <- 2*pnorm(log2FC_real, mean_fc, stddev)
  p_value <- c(p_value, rep(10^(-50), (length(genes_pos_inf)+length(genes_neg_inf))))
  p_vals[[count]] <- p_value
  
  bool_p <- (p_value <= 0.05)
  genes_sig <- genes_real[bool_p]
  p_values_keep <- p_vals[bool_p]
  fold_change_sig <- FC_comb[bool_p]
  
  sig_FCS[[count]] <- fold_change_sig
  genes_sig_list[[count]] <- genes_sig
  
  hist(log2FC_real, probability = T, breaks = 50,
       main = count,
       xlab = 'Log2 of Fold Change')
  
  count = count + 1
}

dev.off()

export(genes_sig_list, '~jennadiegel/Documents/SURE/DES/TCGA_data/genes_sig.xlsx')
export(sig_FCS, '~jennadiegel/Documents/SURE/DES/TCGA_data/fold_change.xlsx')


