#### Imports ####

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
titles <- list("Log2 distribution of BxPC3","Log2 distribution of Capan1","Log2 distribution of PANC1","Log2 distribution of TU8902","Log2 distribution of TU8988T","Log2 distribution of UM2","Log2 distribution of UM90")

#initiate lists to hold fold changes and gene names outside the loop
all_data <- list()
all_fold_changes_log2 <- list()
gene_names_list <- list()
inf_gene_names <- list()
z_scores <- list()
p_values <- list()
count <- 1

for (cell in cells) {
  #separate out gene names
  gene_names <- cell[ ,1]
  
  #separate out KO and WT samples
  cell_KO <- cell[ ,2:4]
  cell_WT <- cell[ ,5:7]
  
  #get row-wise means for KO and WT samples
  cell_KO_mean <- rowMeans(cell_KO)  
  cell_WT_mean <- rowMeans(cell_WT)  
  
  #get fold change KO/WT
  cell_fold_change <- cell_KO_mean/cell_WT_mean #getting NaN with WT = 0 - set to max log 2 fold change - P value is arbitrarily large number - add after pvalue calculation
  
  #transform fold change to get normal distribution using log2
  cell_fold_change_log2 <- log2(cell_fold_change)
  
  #find NaN and inf values and keep only the non-NaN and non-inf values
  bool <- (is.nan(cell_fold_change_log2)|is.infinite(cell_fold_change_log2))
  cell_fold_change_keep <- cell_fold_change[!bool]
  cell_fold_change_log2_keep <- cell_fold_change_log2[!bool]
  gene_names_keep <- gene_names[!bool]
  
  #separate out inf values 
  bool_inf <- is.infinite(cell_fold_change_log2)
  gene_names_inf <- gene_names[bool_inf]
  
  #calculate z score and p value for non inf values
  z__score <- scale(cell_fold_change_log2_keep)
  z_scores[[count]] <- z__score
  
  p__value <- 2*pnorm(-abs(z__score))
  p_values[[count]] <- p__value
  
  data_storage <- data.frame(
    gene_name <- gene_names_keep,
    fold_change <- cell_fold_change_keep,
    log2_fold_change <- cell_fold_change_log2_keep,
    z_score <- z__score,
    p_value <- p__value)
  
  all_data[[count]] <- data_storage
  
  #all_fold_changes_log2[[count]] <- cell_fold_change_keep #redundant line to work with data outside loop
  gene_names_list[[count]] <- gene_names_keep
  inf_gene_names[[count]] <- gene_names_inf
  hist(cell_fold_change_log2_keep, probability = T, breaks = 50,
    main = titles[[count]])
    
  count <- count + 1}


