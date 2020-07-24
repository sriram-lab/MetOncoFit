#DEA: diffential expression analysis 

# Perform differential expression analysis using RNAseq data. Determines fold change values, 
# log2 of fold change values, and displays data with up to 2 different naming schemes if 
# provided in the data.

# DEA separates out significant and non-significant genes using a p-value of 0.05 after 
# transforming the data using log2 of the fold changes.

# Inputs:
# symbols -> list of gene names/symbols
# WT_mean -> vector of wildtype RNAseq values for each gene
# MUT_mean -> vector of mutant RNAseq values for each gene
# other_IDs -> optional, another ID type mapped to the symbols/names

# Output:
# list containing sig_data and non_sig_data
# sig_data -> dataframe of significant genes and their fold change values, 
# log2 fold change values, p values, and other ID if provided as an input
# non_sig_data -> dataframe of non-significant genes and their fold change values, 
# log2 fold change values, p values, and other ID if provided as an input

# @author: Jenna Diegel

DEA <- function(symbols, WT_mean, MUT_mean, other_IDs) {
  
  # Allow for possibility of mapping to another ID type
  if(missing(other_IDs)) {
    other_IDs <- symbols
    missing_IDS <- TRUE
  } else {
    missing_IDS <- FALSE
  }
  
  # get the fold change values for all genes
  FC = MUT_mean/WT_mean
  log2FC = log2(FC)
  
  # find where there are infinate values
  bool_inf_pos <- (is.infinite(log2FC) & (sign(log2FC) == 1))
  bool_inf_neg <- (is.infinite(log2FC) & (sign(log2FC) == -1))
  
  # get list of symbols and IDs for infinite values
  symbols_pos_inf <- symbols[bool_inf_pos]
  IDS_pos_inf <- other_IDs[bool_inf_pos]
  symbols_neg_inf <- symbols[bool_inf_neg]
  IDS_neg_inf <- other_IDs[bool_inf_neg]

  # find where NaN values occur (0/0 = NaN) and start list of non-significant
  # genes, FCs, etc
  non_sig_symbols <- symbols[is.nan(log2FC)]
  non_sig_FC <- FC[is.nan(log2FC)]
  non_sig_log2FC <- log2FC[is.nan(log2FC)]
  non_sig_IDS <- other_IDs[is.nan(log2FC)]
  
  # find where the non-real FCs are and separate out to leave only real values
  # for analysis
  bool <- (is.nan(log2FC)|is.infinite(log2FC))
  FC_real <- FC[!bool]
  log2FC_real <- log2FC[!bool]
  symbols_real <- symbols[!bool]
  IDS_real <- other_IDs[!bool]
  
  # find the mean and standard deviation to get the 2 sided p values of the log2FCs
  mean_fc <- mean(log2FC_real)
  stddev <- sd(log2FC_real)
  p_value <- 2*pnorm(log2FC_real, mean_fc, stddev)
  
  # make list of p-values, FCs, log2FCs, symbols, and other IDs
  # infinite values -> get arbitrarily small p-value and FC = max/min value seen
  # depending on if it is -inf or +inf
  p_value_comb <- c(p_value, rep(10^(-50), (length(symbols_pos_inf)+length(symbols_neg_inf))))
  FC_comb <- c(FC_real, rep(max(FC_real),length(symbols_pos_inf)), rep(min(FC_real),length(symbols_neg_inf)))
  log2FC_comb <- c(log2FC_real, rep(max(log2FC_real),length(symbols_pos_inf)), rep(min(log2FC_real),length(symbols_neg_inf)))
  symbols_comb <- c(symbols_real, symbols_pos_inf, symbols_neg_inf)
  IDS_comb <- c(IDS_real, IDS_pos_inf, IDS_neg_inf)
  
  # find where there are significant p values
  bool_p <- (p_value_comb <= 0.05)
  
  # separate all non significant data
  # NaN values get an arbitrarily large p-value (the max seen in the set)
  non_sig_p_values <- c(rep(max(p_value_comb), length(non_sig_symbols)), p_value_comb[!bool_p])
  non_sig_symbols <- c(non_sig_symbols, symbols_comb[!bool_p])
  non_sig_FC <- c(non_sig_FC, FC_comb[!bool_p])
  non_sig_log2FC <- c(non_sig_log2FC, log2FC_comb[!bool_p])
  non_sig_IDS <- c(non_sig_IDS, IDS_comb[!bool_p])
  
  # set NaN FC=1 (0/0 -> no change in FC) and log2FC=0 (log1=0)
  non_sig_FC[is.nan(non_sig_FC)] = 1
  non_sig_log2FC[is.nan(non_sig_log2FC)] = 0
  
  # separate out significant data
  p_values_sig <- p_value_comb[bool_p]
  symbols_sig <- symbols_comb[bool_p]
  FC_sig <- FC_comb[bool_p]
  log2FC_sig <- log2FC_comb[bool_p]
  sig_IDS <- IDS_comb[bool_p]
  
  # create dataframes to store data, depending on if otherID is provided
  if(missing_IDS) {
    sig_data = data.frame("Symbols" = symbols_sig, "FC" = FC_sig, "log2FC" = log2FC_sig,
                          "Pvalues" = p_values_sig)
    non_sig_data = data.frame("Symbols" = non_sig_symbols, "FC" = non_sig_FC, "log2FC" = non_sig_log2FC,
                              "Pvalues" = non_sig_p_values)
    
  } else {
    sig_data = data.frame("Symbols" = symbols_sig, "FC" = FC_sig, "log2FC" = log2FC_sig,
                          "Pvalues" = p_values_sig, "otherID" = sig_IDS)
    non_sig_data = data.frame("Symbols" = non_sig_symbols, "FC" = non_sig_FC, "log2FC", non_sig_log2FC,
                              "Pvalues" = non_sig_p_values, "otherID" = non_sig_IDS)
  }
  
  # create list to output both dataframes, R only allows one output
  return(list(sig_data = sig_data, non_sig_data = non_sig_data))
  
}
