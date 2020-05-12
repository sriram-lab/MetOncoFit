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
all_fold_changes_log2 <- list()
count = 1

for (cell in cells) {
  #separate out KO and WT samples
  cell_KO <- cell[,2:4]
  cell_WT <- cell[,5:7]
  
  #get row-wise means for KO and WT samples
  cell_KO_mean <- rowMeans(cell_KO)  
  cell_WT_mean <- rowMeans(cell_WT)  
  
  #get fold change WT/KO
  cell_fold_change <- cell_KO_mean/cell_WT_mean
  
  #transform fold change to get normal distribution using log2
  cell_fold_change_log2 <- log2(cell_fold_change)
  
  all_fold_changes_log2[[count]] <- cell_fold_change_log2
  hist(all_fold_changes_log2[[count]], probability = T, breaks = 50,
    main = titles[[count]])
    
  count <- count + 1
  
}

# #separate fold changes out of list 
# #note: these are the log 2 of fold change
# bx_fc <- all_fold_changes_log2[[1]]
# cp_fc <- all_fold_changes_log2[[2]]
# pn_fc <- all_fold_changes_log2[[3]]
# tu_fc <- all_fold_changes_log2[[4]]
# tut_fc <- all_fold_changes_log2[[5]]
# um2_fc <- all_fold_changes_log2[[6]]
# um90_fc <- all_fold_changes_log2[[7]]
# 
# #histograms to look at distributions
# hist(bx_fc, probability = T, breaks = 30)
# hist(cp_fc, probability = T, breaks = 30)
# hist(pn_fc, probability = T, breaks = 50)
# hist(tu_fc, probability = T, breaks = 30)
# hist(tut_fc, probability = T, breaks = 30)
# hist(um2_fc, probability = T, breaks = 50)
# hist(um90_fc, probability = T, breaks = 50)

bx_fc_t <- t(bx_fc)
sd(bx_fc_t, na.rm = TRUE)
