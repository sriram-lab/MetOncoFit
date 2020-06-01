%% Lyssiotis %%
% Import gene names and fold change values

sheets = {'Sheet 1', 'Sheet 2', 'Sheet 3', 'Sheet 4', 'Sheet 5', 'Sheet 6', 'Sheet 7'};
genes = cell(1,length(sheets));
fc = cell(1,length(sheets));

for sheet = 1:length(sheets)

    [~,genes{1, sheet},~] = xlsread(fullfile('..','Lyssiotis_data','genes_sig.xlsx'),sheets{sheet},'A:A');
    [fc{1, sheet},~ ,~] = xlsread(fullfile('..','Lyssiotis_data','fold_change.xlsx'), sheets{sheet}, 'A:A');
    
end 

%CL = cell line specific data
%pan = pan pancreatic data

% RECON1
load(fullfile('..','models','model_human_duarte.mat'));
RECON1_genes = metabolicmodel.genes;

[L_RECON1_CL_matches, L_RECON1_pan] = mapping(RECON1_genes, genes, fc);

sound(sin(1:3000))

% %% RECON2
% load(fullfile('..','models','Recon2.v05.mat')); 
% RECON2_genes = modelR205.genes;
% 
% %Need to convert RECON2 naming scheme 
% %convert to new style of naming
% for a = 1:length(genes)
%     
%     genes{a} = replace(genes{a},".","_AT");
% 
% end
% 
% %import list of RECON3D genes with bigg ids and names
% [~, bigg_id, ~] = xlsread(fullfile('gene_name_info','BiGG_ID_to_Name.xlsx'),'A:A');
% [~, name, ~] = xlsread(fullfile('gene_name_info','BiGG_ID_to_Name.xlsx'),'B:B');
% 
% gene_str = join(genes); 
% %add space at beginning to use as an indicator of a new word
% gene_str = strcat({' '}, gene_str, {' '}); 
% 
% %initialize an empty cell array to input name matches
% matches = cell(1, length(bigg_id));
% bigg_ids = cell(1, length(bigg_id));
% for a = 1:length(bigg_id)
%         
%     reg = strcat('\s',bigg_id{a},'\s');
%     %find matches and for genes listed more than once, take only one match
%     matches(1, a) = regexp(gene_str,reg,'match','once');
%     bigg_ids(1,a) = bigg_id{a};
% end
% 
% %remove empty cell contents
% bigg_ids = bigg_ids(cellfun('isempty', matches)); 
% unmatched = name(cellfun('isempty', matches)); 
% model_matches = name(~cellfun('isempty', matches)); 
% 
% [RECON2_CL_matches, RECON2_pan] = mapping(model_matches, genes, fc);
% 
% %% RECON3D
% [~,RECON3_genes,~] = xlsread(fullfile('gene_name_info','BiGG_ID_to_Name.xlsx'),'B:B');
% 
% [RECON3_CL_matches, RECON3_pan] = mapping(RECON3_genes, genes, fc);
% 
% %% Human1 
% [~, Human1_genes, ~] = xlsread(fullfile('gene_name_info','gene_name_conversion.xlsx'),'A:A');
% 
% [Human1_CL_matches, Human1_pan] = mapping(Human1_genes, genes, fc);

%% CCLE %%

% import data
%start at 3, do not want first 2 cell lines, end at 33, a total of 33 cell
%lines

sheet = cell(1,31);
count = 1;
for a = 3:33
    
    sheet_name = strcat('Sheet', {' '}, num2str(a));
    sheet(1,count) = sheet_name;
    count = count + 1;
    
end

genes = cell(1,length(sheet));
fcs = cell(1,length(sheet));
for a = 1:length(sheet)
    genes{1,a} = readcell(fullfile('..','CCLE_data','sig_genes.xlsx'),'Sheet',sheet{a});
    fcs{1,a} = readcell(fullfile('..','CCLE_data','fold_changes.xlsx'),'Sheet',sheet{a});
end

load(fullfile('..','models','model_human_duarte.mat'));
RECON1_genes = metabolicmodel.genes;

[C_RECON1_CL_matches, C_RECON1_pan] = mapping(RECON1_genes, genes, fcs);

%% TCGA %%

%lengths are 32 for 32 mutant cell lines
genes = cell(1,32);
fc = cell(1,32);

for a = 1:32
    
    sheet = strcat('Sheet', {' '}, num2str(a));
    
    [~,genes{1, a},~] = xlsread(fullfile('..','TCGA_data','genes_sig.xlsx'),sheet{1} ,'A:A');
    [fc{1, a},~ ,~] = xlsread(fullfile('..','TCGA_data','fold_change.xlsx'), sheet{1}, 'A:A');
    
end 

%CL = cell line specific data
%pan = pan pancreatic data

% RECON1
load(fullfile('..','models','model_human_duarte.mat'));
RECON1_genes = metabolicmodel.genes;

[T_RECON1_CL_matches, T_RECON1_pan] = mapping(RECON1_genes, genes, fc);

sound(sin(1:3000))