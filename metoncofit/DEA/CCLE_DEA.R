#### imports ####
ccle = import("~jennadiegel/Documents/SURE/DES/CCLE_data.csv")

gene_names <- ccle[,2]
wt_vals <- ccle[ ,3:4]
mut_vals <- ccle[ ,-(1:4)]

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
inf_gene_names <- list()
z_scores <- list()
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
  
  #separate out inf values 
  bool_inf <- is.infinite(fold_change_log2)
  gene_names_inf <- gene_names[bool_inf]
  
  
  z__score <- scale(fold_change_log2_keep)
  z_scores[[i]] <- z__score
  
  p__value <- 2*pnorm(-abs(z__score))
  p_values[[i]] <- p__value
  
  data_storage <- data.frame(
    gene_name <- c(gene_names_keep, gene_names_inf),
    fold_change <- c(fold_change_keep, rep(max(fold_change_keep),length(gene_names_inf))),
    log2_fold_change <- c(fold_change_log2_keep, rep(max(fold_change_log2_keep),length(gene_names_inf))),
    z_score <- c(z__score, rep(max(z__score),length(gene_names_inf))),
    p_value <- c(p__value, rep(10^(-50),length(gene_names_inf))))
  
  all_data[[i]] <- data_storage
  fold_changes_log2[[i]] <- fold_change_log2_keep
  fold_changes[[i]] <- fold_change_keep
  inf_gene_names[[i]] <- gene_names_inf
  
  hist(fold_change_log2_keep, probability = T, breaks = 50,
       main = i)
  
  i = i + 1}

export(all_data, 'CCLE.xlsx')

##figure out why cell line 1 and 2 have weird distributions
##could be modeled with an exponential distribution?
##ignore first two MUT samples
