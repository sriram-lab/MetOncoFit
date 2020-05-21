#### imports ####
library(rio)
ccle = import("~jennadiegel/Documents/SURE/DES/data/CCLE_data.csv")

gene_names <- ccle[,2]
wt_vals <- ccle[ ,3:4]
mut_vals <- ccle[ ,-(1:4)]

#### plot titles ####
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

#### analysis ####

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
p_values <- list()
all_data <- list()
i = 1
for (i in 1:length(mut_vals)) {
  
  
  fold_change <- mut_vals[ ,i]/WT_mean
  fold_change_log2 <- log2(fold_change)
  
  #find NaN and inf values and keep only the non-NaN and non-inf values
  bool <- (is.nan(fold_change_log2)|is.infinite(fold_change_log2))
  fold_change_keep <- fold_change[!bool]
  fold_change_log2_keep <- fold_change_log2[!bool]
  gene_names_keep <- gene_names[!bool]
  
  #separate out inf values, +inf = upregulation, -inf = downregulation
  bool_inf_pos <- (is.infinite(fold_change_log2) & sign(fold_change_log2 == 1))
  bool_inf_neg <- (is.infinite(fold_change_log2) & sign(fold_change_log2 == -1))
  gene_names_pos_inf <- gene_names[bool_inf_pos]
  gene_names_neg_inf <- gene_names[bool_inf_neg]
  
  mean_fc <- mean(fold_change_log2_keep)
  stddev <- sd(fold_change_log2_keep)
  p__value <- 2*pnorm(fold_change_log2_keep,mean_fc,stddev)
  p_values[[i]] <- p__value
  
  data_storage <- data.frame(
    gene_name <- c(gene_names_keep, gene_names_pos_inf, gene_names_neg_inf),
    fold_change <- c(fold_change_keep, rep(max(fold_change_keep),length(gene_names_pos_inf)), rep(min(fold_change_keep),length(gene_names_neg_inf))),
    log2_fold_change <- c(fold_change_log2_keep, rep(max(fold_change_log2_keep),length(gene_names_pos_inf)),rep(min(fold_change_log2_keep),length(gene_names_neg_inf))),
    p_value <- c(p__value, rep(10^(-50),(length(gene_names_pos_inf) + length(gene_names_neg_inf)))))
  
  all_data[[i]] <- data_storage
  fold_changes_log2[[i]] <- fold_change_log2_keep
  fold_changes[[i]] <- fold_change_keep
  upreg_inf_genes[[i]] <- gene_names_pos_inf
  downreg_inf_genes[[i]] <- gene_names_neg_inf
  
  hist(fold_change_log2_keep, probability = T, breaks = 50,
       main = titles[[i]],
       xlab = 'Log2 of Fold Change')
  
  i = i + 1}

export(all_data, '~jennadiegel/Documents/SURE/DES/CCLE_data/CCLE.xlsx')
dev.off()

##figure out why cell line 1 and 2 have weird distributions
##could be modeled with an exponential distribution?
##ignore first two MUT samples

