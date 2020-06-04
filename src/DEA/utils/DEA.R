#DEA: diffential expression analysis 

# Perform differential expression analysis using RNAseq data. Determines fold change values, 
# log2 of fold change values, and displays data with up to 2 different naming schemes if 
# provided in the data.

# DEA separates out significant and non-significant genes using a p-value of 0.05 after 
# transforming the data using log2 of the fold changes.

# @author: Jenna Diegel

DEA <- function(symbols, WT_mean, KO_mean, other_IDs) {
  
  if(missing(other_IDs)) {
    other_IDs <- symbols
    missing_IDS <- TRUE
  } else {
    missing_IDS <- FALSE
  }
  
  FC = KO_mean/WT_mean
  log2FC = log2(FC)
  bool_inf_pos <- (is.infinite(log2FC) & (sign(log2FC) == 1))
  bool_inf_neg <- (is.infinite(log2FC) & (sign(log2FC) == -1))
  
  symbols_pos_inf <- symbols[bool_inf_pos]
  IDS_pos_inf <- other_IDs[bool_inf_pos]
  symbols_neg_inf <- symbols[bool_inf_neg]
  IDS_neg_inf <- other_IDs[bool_inf_neg]

  non_sig_symbols <- symbols[is.nan(log2FC)]
  non_sig_FC <- FC[is.nan(log2FC)]
  non_sig_log2FC <- log2FC[is.nan(log2FC)]
  non_sig_IDS <- other_IDs[is.nan(log2FC)]
  
  bool <- (is.nan(log2FC)|is.infinite(log2FC))
  FC_real <- FC[!bool]
  log2FC_real <- log2FC[!bool]
  symbols_real <- symbols[!bool]
  IDS_real <- other_IDs[!bool]
    
  mean_fc <- mean(log2FC_real)
  stddev <- sd(log2FC_real)
  p_value <- 2*pnorm(log2FC_real, mean_fc, stddev)
  
  p_value_comb <- c(p_value, rep(10^(-50), (length(symbols_pos_inf)+length(symbols_neg_inf))))
  FC_comb <- c(FC_real, rep(max(FC_real),length(symbols_pos_inf)), rep(min(FC_real),length(symbols_neg_inf)))
  log2FC_comb <- c(log2FC_real, rep(max(log2FC_real),length(symbols_pos_inf)), rep(min(log2FC_real),length(symbols_neg_inf)))
  symbols_comb <- c(symbols_real, symbols_pos_inf, symbols_neg_inf)
  IDS_comb <- c(IDS_real, IDS_pos_inf, IDS_neg_inf)
  
  bool_p <- (p_value_comb <= 0.05)
  
  non_sig_p_values <- c(rep(max(p_value_comb), length(non_sig_symbols)), p_value_comb[!bool_p])
  non_sig_symbols <- c(non_sig_symbols, symbols_comb[!bool_p])
  non_sig_FC <- c(non_sig_FC, FC_comb[!bool_p])
  non_sig_log2FC <- c(non_sig_log2FC, log2FC_comb[!bool_p])
  non_sig_IDS <- c(non_sig_IDS, IDS_comb[!bool_p])
  
  non_sig_FC[is.nan(non_sig_FC)] = 1
  non_sig_log2FC[is.nan(non_sig_log2FC)] = 0
  
  p_values_sig <- p_value_comb[bool_p]
  symbols_sig <- symbols_comb[bool_p]
  FC_sig <- FC_comb[bool_p]
  log2FC_sig <- log2FC_comb[bool_p]
  sig_IDS <- IDS_comb[bool_p]
  
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
  
  return(list(sig_data = sig_data, non_sig_data = non_sig_data))
  
}
