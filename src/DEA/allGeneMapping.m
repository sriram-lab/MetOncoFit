%% naming schemes for variables
%{ 
CL = cell line specific data
pan = pan pancreatic data
NS = non significant genes
L = lyssiotis
C = CCLE
T = TCGA
%}

%% RUN FIRST! load models, files, and convert naming scheme if necessary
load(fullfile('..','models','model_human_duarte.mat'));
RECON1_genes = metabolicmodel.genes;

load(fullfile('..','models','Recon2.v05.mat')); 
RECON2_genes = modelR205.genes;

%Need to convert RECON2 naming scheme 
%convert to new style of naming
for a = 1:length(RECON2_genes)
    
    RECON2_genes{a} = replace(RECON2_genes{a},".","_AT");

end

[~, bigg_id, ~] = xlsread(fullfile('gene_name_info','BiGG_ID_to_Name.xlsx'),'A:A');
[~, name, ~] = xlsread(fullfile('gene_name_info','BiGG_ID_to_Name.xlsx'),'B:B');

load(fullfile('..','models','Recon3D.mat'));
RECON3_genes = Recon3D.genes;

load(fullfile('..','models','HumanGEM.mat'));
Human1_genes = ihuman.genes;

%% Lyssiotis %%
% Import gene names and fold change values

sheets = {'Sheet 1', 'Sheet 2', 'Sheet 3', 'Sheet 4', 'Sheet 5', 'Sheet 6', 'Sheet 7'};
%row 1 = sig genes, row 2 = non sig genes
genes = cell(2,length(sheets));
fcs = cell(2,length(sheets));

% unpack data from excel sheets
for sheet = 1:length(sheets)
    
    Lyss_sig = readcell(fullfile('..','Lyssiotis_data','Lyssiotis_data_sig.xlsx'),'Sheet',sheets{sheet},'Range',2);
    Lyss_non_sig = readcell(fullfile('..','Lyssiotis_data','Lyssiotis_data_non_sig.xlsx'),'Sheet',sheets{sheet},'Range',2);
    genes{1, sheet} = Lyss_sig(:,1);
    genes{2, sheet} = Lyss_non_sig(:,1);
    fcs{1, sheet} = cell2mat(Lyss_sig(:,2));
    fcs{2, sheet} = cell2mat(Lyss_non_sig(:,2));
    
end 

%% RECON1

% map significant and non significant genes
[L_RECON1_CL_matches, L_RECON1_pan] = mapping(RECON1_genes, genes(1,:), fcs(1,:), true);
[L_RECON1_CL_matches_NS, L_RECON1_pan_NS] = mapping(RECON1_genes, genes(2,:), fcs(2,:), false);

%% RECON2

all_BIGG = cell(2,length(genes));
all_FC = cell(2,length(genes));

% convert bigg ID to gene symbol 
%b = 1: sig genes, b = 2: non sig genes
for b = 1:2
    for a = 1:length(genes)
        [~, bigg_to_FC] = paired_mapping(name, bigg_id, genes{b,a}, fcs{b,a});
        all_BIGG{b,a} = bigg_to_FC(:,1);
        all_FC{b,a} = cell2mat(bigg_to_FC(:,2));
    end
end

%% 
% map significant and non significant genes
[L_RECON2_CL_matches, L_RECON2_pan] = mapping(RECON2_genes, all_BIGG(1,:), all_FC(1,:), true);
[L_RECON2_CL_matches_NS, L_RECON2_pan_NS] = mapping(RECON2_genes, all_BIGG(2,:), all_FC(2,:), false);

%% RECON3D
% map significant and non significant genes
[L_RECON3_CL_matches, L_RECON3_pan] = mapping(RECON3_genes, all_BIGG(1,:), all_FC(1,:), true);
[L_RECON3_CL_matches_NS, L_RECON3_pan_NS] = mapping(RECON3_genes, all_BIGG(2,:), all_FC(2,:), false);

%% Human1 

% load data to convert to symbol from ENSG ID
ENSG_to_sym = readtable(fullfile('gene_name_info','ENS_to_symbol.txt'));
ENSG = ENSG_to_sym.ENSG_ID;
symbol = ENSG_to_sym.symbol;

all_ENSG = cell(2,length(genes));
all_FC = cell(2,length(genes));

% convert from ENSG ID to gene symbol
for b = 1:2
    for a = 1:length(genes)
        [~, bigg_to_FC] = paired_mapping(symbol, ENSG, genes{b,a}, fcs{b,a});
        all_ENSG{b,a} = bigg_to_FC(:,1);
        all_FC{b,a} = cell2mat(bigg_to_FC(:,2));
    end
end

% map significant and non significant genes
[L_human1_CL_matches, L_human1_pan] = mapping(Human1_genes, all_ENSG(1,:), all_FC(1,:), true);
[L_human1_CL_matches_NS, L_human1_pan_NS] = mapping(Human1_genes, all_ENSG(2,:), all_FC(2,:), false);

%% CCLE %%

% import data
% 33 for 33 cell lines
ENSGs = cell(2,33);
genes = cell(2,33);
fcs = cell(2,33);

% unpack data from excel file
for a = 1:33
    sheet = strcat('Sheet', {' '}, num2str(a));
    
    CCLE_sig = readcell(fullfile('..','CCLE_data','CCLE_data_sig.xlsx'),'Sheet',sheet{1},'Range',2);
    CCLE_non_sig = readcell(fullfile('..','CCLE_data','CCLE_data_non_sig.xlsx'),'Sheet',sheet{1},'Range',2);
    genes{1, a} = CCLE_sig(:,1);
    genes{2, a} = CCLE_non_sig(:,1);
    fcs{1, a} = cell2mat(CCLE_sig(:,2));
    fcs{2, a} = cell2mat(CCLE_non_sig(:,2));
    ENSGs{1, a} = CCLE_sig(:,5);
    ENSGs{2,a} = CCLE_sig(:,5);
    
end

% Recon1
% map significant and non significant genes
[C_RECON1_CL_matches, C_RECON1_pan] = mapping(RECON1_genes, genes(1,:), fcs(1,:), true);
[C_RECON1_CL_matches_NS, C_RECON1_pan_NS] = mapping(RECON1_genes, genes(2,:), fcs(2,:), false);

% Recon2
all_BIGG = cell(2,length(genes));
all_FC = cell(2,length(genes));

for b = 1:2
    for a = 1:length(genes)
        [~, bigg_to_FC] = paired_mapping(name, bigg_id, genes{b,a}, fcs{b,a});
        all_BIGG{b,a} = bigg_to_FC(:,1);
        all_FC{b,a} = cell2mat(bigg_to_FC(:,2));
    end
end

% map significant and non significant genes
[C_RECON2_CL_matches, C_RECON2_pan] = mapping(RECON2_genes, all_BIGG(1,:), all_FC(1,:), true);
[C_RECON2_CL_matches_NS, C_RECON2_pan_NS] = mapping(RECON2_genes, all_BIGG(2,:), all_FC(2,:), false);

% Recon 3
% map significant and non significant genes
[C_RECON3_CL_matches, C_RECON3_pan] = mapping(RECON3_genes, all_BIGG(1,:), all_FC(1,:), true);
[C_RECON3_CL_matches_NS, C_RECON3_pan_NS] = mapping(RECON3_genes, all_BIGG(2,:), all_FC(2,:), false);

% Human1
% map significant and non significant genes
[C_human1_CL_matches, C_human1_pan] = mapping(Human1_genes, ENSGs(1,:), fcs(1,:), true);
[C_human1_CL_matches_NS, C_human1_pan_NS] = mapping(Human1_genes, ENSGs(2,:),fcs(2,:), false);

%% TCGA %%

%lengths are 32 for 32 mutant cell lines
genes = cell(2,32);
fcs = cell(1,32);
ENSGs = cell(2,32);

% unpack data frome excel file
for a = 1:32
    
    sheet = strcat('Sheet', {' '}, num2str(a));
    
    TCGA_sig = readcell(fullfile('..','TCGA_data','TCGA_data_sig.xlsx'),'Sheet',sheet{1},'Range',2);
    TCGA_non_sig = readcell(fullfile('..','TCGA_data','TCGA_data_non_sig.xlsx'),'Sheet',sheet{1},'Range',2);
    genes{1, a} = TCGA_sig(:,1);
    genes{2, a} = TCGA_non_sig(:,1);
    fcs{1, a} = cell2mat(TCGA_sig(:,2));
    fcs{2, a} = cell2mat(TCGA_non_sig(:,2));
    ENSGs{1, a} = TCGA_sig(:,5);
    ENSGs{2,a} = TCGA_sig(:,5);
end 

% RECON1
% map significant and non significant genes
[T_RECON1_CL_matches, T_RECON1_pan] = mapping(RECON1_genes, genes(1,:), fcs(1,:), true);
[T_RECON1_CL_matches_NS, T_RECON1_pan_NS] = mapping(RECON1_genes, genes(1,:), fcs(2,:), false);

% Recon2
all_BIGG = cell(1,length(genes));
all_FC = cell(1,length(genes));

% convert BIGG ID to gene symbol
for b = 1:2
    for a = 1:length(genes)
        [~, bigg_to_FC] = paired_mapping(name, bigg_id, genes{b,a}, fcs{b,a});
        all_BIGG{b,a} = bigg_to_FC(:,1);
        all_FC{b,a} = cell2mat(bigg_to_FC(:,2));
    end
end

% map significant and non significant genes
[T_RECON2_CL_matches, T_RECON2_pan] = mapping(RECON2_genes, all_BIGG(1,:), all_FC(1,:), true);
[T_RECON2_CL_matches_NS, T_RECON2_pan_NS] = mapping(RECON2_genes, all_BIGG(2,:), all_FC(2,:), false);

% Recon3
% map significant and non significant genes
[T_RECON3_CL_matches, T_RECON3_pan] = mapping(RECON3_genes, all_BIGG(1,:), all_FC(1,:), true);
[T_RECON3_CL_matches_NS, T_RECON3_pan_NS] = mapping(RECON3_genes, all_BIGG(2,:), all_FC(2,:), false);

 % Human1
 % map significant and non significant genes
[T_human1_CL_matches, T_human1_pan] = mapping(Human1_genes, ENSGs(1,:), fcs(1,:), true);
[T_human1_CL_matches_NS, T_human1_pan_NS] = mapping(Human1_genes, ENSGs(2,:), fcs(2,:), false);