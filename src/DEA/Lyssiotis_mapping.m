%% Import gene names and fold change values

sheets = {'Sheet 1', 'Sheet 2', 'Sheet 3', 'Sheet 4', 'Sheet 5', 'Sheet 6', 'Sheet 7'};
genes = cell(1,length(sheets));
fc = cell(1,length(sheets));

for sheet = 1:length(sheets)

    [~,genes{1, sheet},~] = xlsread(fullfile('..','Lyssiotis_data','genes_sig.xlsx'),sheets{sheet},'A:A');
    [fc{1, sheet},~ ,~] = xlsread(fullfile('..','Lyssiotis_data','fold_change.xlsx'), sheets{sheet}, 'A:A');
    
end 

%CL = cell line specific data
%pan = pan pancreatic data

%% RECON1
load(fullfile('..','models','model_human_duarte.mat'));
RECON1_genes = metabolicmodel.genes;

[RECON1_CL_matches, RECON1_pan] = mapping(RECON1_genes, genes, fc);

%% RECON2
load(fullfile('..','models','Recon2.v05.mat')); 
RECON2_genes = modelR205.genes;

%Need to convert RECON2 naming scheme 
%convert to new style of naming
for a = 1:length(genes)
    
    genes{a} = replace(genes{a},".","_AT");

end

%import list of RECON3D genes with bigg ids and names
[~, bigg_id, ~] = xlsread(fullfile('gene_name_info','BiGG_ID_to_Name.xlsx'),'A:A');
[~, name, ~] = xlsread(fullfile('gene_name_info','BiGG_ID_to_Name.xlsx'),'B:B');

gene_str = join(genes); 
%add space at beginning to use as an indicator of a new word
gene_str = strcat({' '}, gene_str, {' '}); 

%initialize an empty cell array to input name matches
matches = cell(1, length(bigg_id));

for a = 1:length(bigg_id)
        
    reg = strcat('\s',bigg_id{a},'\s');
    %find matches and for genes listed more than once, take only one match
    matches(1, a) = regexp(gene_str,reg,'match','once');
    
end

%remove empty cell contents
unmatched = name(cellfun('isempty', matches)); 
model_matches = name(~cellfun('isempty', matches)); 

[RECON2_CL_matches, RECON2_pan] = mapping(RECON2_genes, genes, fc);

%% RECON3D
[~,RECON3_genes,~] = xlsread(fullfile('gene_name_info','BiGG_ID_to_Name.xlsx'),'B:B');

[RECON3_CL_matches, RECON3_pan] = mapping(RECON3_genes, genes, fc);

%% Human1 
[~, Human1_genes, ~] = xlsread(fullfile('gene_name_info','gene_name_conversion.xlsx'),'A:A');

[Human1_CL_matches, Human1_pan] = mapping(Human1_genes, genes, fc);
