#### Imports ####

library(rio)
# import BxPC3 data as bx
bx <- import("~jennadiegel/Documents/SURE/DES/CCLs_BxPC3.csv")
# import Capan1 data as cp
cp <- import("~jennadiegel/Documents/SURE/DES/CCLs_Capan1.csv")
# import PANC1 data as pn
pn <- import("~jennadiegel/Documents/SURE/DES/CCLs_PANC1.csv")
# import TU8902 data as tu
tu <- import("~jennadiegel/Documents/SURE/DES/CCLs_TU8902.csv")
# import TU8988T data as tut
tut <- import("~jennadiegel/Documents/SURE/DES/CCLs_TU8988T.csv")
# import UM2 data as um2
um2 <- import("~jennadiegel/Documents/SURE/DES/CCLs_UM2.csv")
# import UM90 data as um90
um90 <- import("~jennadiegel/Documents/SURE/DES/CCLs_UM90.csv")

#### analysis ####

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
inf_gene_names <- list()
all_p_values <- list()
z_scores <- list()
p_values <- list()
genes_sig_list <- list()

count <- 1

pdf("Lyssiotis_DEA_histograms.pdf")

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
  bool_inf <- is.infinite(fold_change_log2)
  genes_inf <- genes[bool_inf]
  
  #find NaN and inf values and keep only the non-NaN and non-inf values
  bool <- (is.nan(fold_change_log2)|is.infinite(fold_change_log2))
  fold_change_real <- fold_change[!bool]
  fold_change_log2_real <- fold_change_log2[!bool]
  genes_real <- c(genes[!bool], genes_inf)
  
  #set fold change of inf genes to the max fold change and combine with real values
  fold_change_comb <- c(fold_change_real, rep(max(fold_change_real),length(genes_inf)))
  fold_change_log2_comb <- c(fold_change_log2_real, rep(max(fold_change_log2_real),length(genes_inf)))
  
  p_value <- 2*pnorm(fold_change_log2_real)
  p_value <- c(p_value, rep(10^(-50),length(genes_inf)))
  p_values_all <- p_value
  p_values[[count]] <- p_value
  
  #store values for all genes into a data frame
  data_storage <- data.frame(
    gene_name <- genes_real,
    fold_change <- fold_change_comb,
    log2_fold_change <- fold_change_log2_comb,
    p_value <- p_value)
  
  #find which values are significant and sort out necessary parameters
  bool_p <- (p_value <= 0.05)
  genes_sig <- gene_name[bool_p]
  p_values_keep <- p_values_all[bool_p]
  fold_change_sig <- fold_change_comb[bool_p]
  
  gene_data <- data.frame(
    genes_sig <- genes_sig
  )

  #put data frames and vectors into lists for storage
  all_data[[count]] <- data_storage
  cell_fold_changes[[count]] <- fold_change_sig
  gene_names_list[[count]] <- genes_real
  inf_gene_names[[count]] <- genes_inf
  genes_sig_list[[count]] <- genes_sig
  all_p_values[[count]] <- p_values_keep
  
  #visualize distribution of values as histograms
  hist(fold_change_log2_real, probability = T, breaks = 50,
       main = titles[[count]])
  
  count <- count + 1}

dev.off()

#exports to use data in matlab or python
export(all_data, 'lyossiotis.xlsx')
export(genes_sig_list, 'genes_sig.xlsx')
export(all_p_values, 'p_values.xlsx')
export(cell_fold_changes, 'fold_change.xlsx')



