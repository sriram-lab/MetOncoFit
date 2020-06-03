#### DEA function ####
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

#### Lyssiotis Imports ####

library("rio")
library("readxl")
sheets = list('Sheet1','Sheet2','Sheet3','Sheet4','Sheet5','Sheet6','Sheet7')
bx <- as.data.frame(read_excel("~jennadiegel/Documents/SURE/DES/Lyssiotis_data/TPM_lyss.xlsx", sheet = 'Sheet1'))
cp <- as.data.frame(read_excel("~jennadiegel/Documents/SURE/DES/Lyssiotis_data/TPM_lyss.xlsx", sheet = 'Sheet2'))
pn <- as.data.frame(read_excel("~jennadiegel/Documents/SURE/DES/Lyssiotis_data/TPM_lyss.xlsx", sheet = 'Sheet3'))
tu<- as.data.frame(read_excel("~jennadiegel/Documents/SURE/DES/Lyssiotis_data/TPM_lyss.xlsx", sheet = 'Sheet4'))
tut <- as.data.frame(read_excel("~jennadiegel/Documents/SURE/DES/Lyssiotis_data/TPM_lyss.xlsx", sheet = 'Sheet5'))
um2<- as.data.frame(read_excel("~jennadiegel/Documents/SURE/DES/Lyssiotis_data/TPM_lyss.xlsx", sheet = 'Sheet6'))
um90 <- as.data.frame(read_excel("~jennadiegel/Documents/SURE/DES/Lyssiotis_data/TPM_lyss.xlsx", sheet = 'Sheet7'))

#### Lyssiotis analysis ####

cells <- list(bx, cp, pn, tu, tut, um2, um90)
titles <- list("Log2 distribution of BxPC3","Log2 distribution of Capan1",
               "Log2 distribution of PANC1","Log2 distribution of TU8902",
               "Log2 distribution of TU8988T","Log2 distribution of UM2",
               "Log2 distribution of UM90")
cell_names <- list("BxPC3","Capan1","PANC1","TU8902","TU8988T","UM2","UM90")

#initiate lists to store data 
Lyssiotis_data_sig <- list()
Lyssiotis_data_non_sig <- list()

count <- 1

for (cell in cells) {
  
  #separate out gene names
  genes <- cell[ ,1]
  
  #separate out KO and WT samples
  cell_KO <- cell[ ,2:4]
  cell_WT <- cell[ ,5:7]
  
  #get row-wise means for KO and WT samples
  KO_mean <- rowMeans(cell_KO)
  WT_mean <- rowMeans(cell_WT)  
  
  temp_data <- DEA(genes,WT_mean,KO_mean)
  Lyssiotis_data_sig[[count]] <- temp_data[[1]]
  Lyssiotis_data_non_sig[[count]] <- temp_data[[2]]
  
  count <- count + 1}

export(Lyssiotis_data_sig, '~jennadiegel/Documents/SURE/DES/Lyssiotis_data/Lyssiotis_data_sig.xlsx')
export(Lyssiotis_data_non_sig, '~jennadiegel/Documents/SURE/DES/Lyssiotis_data/Lyssiotis_data_non_sig.xlsx')

#### CCLE Imports ####
library(rio)
ccle = import("~jennadiegel/Documents/SURE/DES/CCLE_data/CCLE_data.csv")

gene_names <- ccle[,2]
wt_vals <- ccle[ ,3:4]
mut_vals <- ccle[ ,-(1:4)]
ENSGs <- ccle[ ,1]

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

#### CCLE analysis ####

#fill in NA values with 0 - database excluded 0 values
wt_vals[is.na(wt_vals)] <- 0
mut_vals[is.na(mut_vals)] <- 0

#calculate WT mean from 2 WT cell lines
WT_mean <- (wt_vals[,1] + wt_vals[,2])/2

CCLE_data_sig <- list()
CCLE_data_non_sig <- list()

count = 1
for (i in 1:length(mut_vals)) {
  
  
  mut <- mut_vals[ ,i]
  temp_data <- DEA(gene_names,WT_mean,mut,ENSGs)
  CCLE_data_sig[[count]] <- temp_data[[1]]
  CCLE_data_non_sig[[count]] <- temp_data[[2]]

  count = count + 1
  i = i + 1}

export(CCLE_data_sig, '~jennadiegel/Documents/SURE/DES/CCLE_data/CCLE_data_sig.xlsx')
export(CCLE_data_non_sig, '~jennadiegel/Documents/SURE/DES/CCLE_data/CCLE_data_non_sig.xlsx')

#### TCGA Imports ####
library("readxl")

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
    TCGA_wt_data[[wt_count]] <- read_excel("~jennadiegel/Documents/SURE/DES/TCGA_data/TCGA_TPM_data.xlsx", sheet = sheet_names[[count]])
    if (count == 1) {
      wt_TPM = matrix(NA, nrow = dim(TCGA_wt_data[[1]])[1], ncol = 6)
    }
    wt_TPM[ ,wt_count] <- TCGA_wt_data[[wt_count]]$TPMs
    wt_count = wt_count + 1
  }
  else {
    TCGA_mut_data[[mut_count]] <- read_excel("~jennadiegel/Documents/SURE/DES/TCGA_data/TCGA_TPM_data.xlsx", sheet = sheet_names[[count]])
    if (count == 7) {
      mut_TPM = matrix(NA, nrow = dim(TCGA_mut_data[[1]])[1], ncol = 32)
    }
    mut_TPM[ ,mut_count] <- TCGA_mut_data[[mut_count]]$TPMs
    mut_count = mut_count + 1
    
  }
  count = count + 1
}

#### TCGA analysis ####

TCGA_data_sig <- list()
TCGA_data_non_sig <- list()

wt_mean <- rowMeans(wt_TPM)

count = 1
for(CL in TCGA_mut_data) {
  
  genes = TCGA_mut_data[[count]]$symbol
  mut = mut_TPM[ ,count]
  ENSG = TCGA_mut_data[[count]]$ensID
  
  temp_data <- DEA(genes,wt_mean,mut,ENSG)
  TCGA_data_sig[[count]] <- temp_data[[1]]
  TCGA_data_non_sig[[count]] <- temp_data[[2]]
  
  count = count + 1
}

export(TCGA_data_sig, '~jennadiegel/Documents/SURE/DES/TCGA_data/TCGA_data_sig.xlsx')
export(TCGA_data_non_sig, '~jennadiegel/Documents/SURE/DES/TCGA_data/TCGA_data_non_sig.xlsx')






