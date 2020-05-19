%% Import gene names and p values
%BxPC3
[~,bx_genes,~] = xlsread('map_genes.xlsx','Sheet 1','A:A');
[bx_fc, ~, ~] = xlsread('fold_change.xlsx', 'Sheet 1', 'A:A');
%Capan1
[~,cp_genes,~] = xlsread('map_genes.xlsx','Sheet 2','A:A');
[cp_fc, ~, ~] = xlsread('fold_change.xlsx', 'Sheet 2', 'A:A');
%PANC1
[~,pn_genes,~] = xlsread('map_genes.xlsx','Sheet 3','A:A');
[pn_fc, ~, ~] = xlsread('fold_change.xlsx', 'Sheet 3', 'A:A');
%TU8902
[~,tu_genes,~] = xlsread('map_genes.xlsx','Sheet 4','A:A');
[tu_fc, ~, ~] = xlsread('fold_change.xlsx', 'Sheet 4', 'A:A');
%TU8988T
[~,tut_genes,~] = xlsread('map_genes.xlsx','Sheet 5','A:A');
[tut_fc, ~, ~] = xlsread('fold_change.xlsx', 'Sheet 5', 'A:A');
%UM2
[~,um2_genes,~] = xlsread('map_genes.xlsx','Sheet 6','A:A');
[um2_fc, ~, ~] = xlsread('fold_change.xlsx', 'Sheet 6', 'A:A');
%UM90
[~,um90_genes,~] = xlsread('map_genes.xlsx','Sheet 7','A:A');
[um90_fc, ~, ~] = xlsread('fold_change.xlsx', 'Sheet 7', 'A:A');
%metabolic model
[~, gene_list, ~] = xlsread('gene_name_conversion.xlsx','A:A');

%% BxPC3 matching

gene_names = join(gene_list);
%add space at beginning to use as an indicator of a new word
gene_names = strcat({' '}, gene_names,{' '});

%initialize an empty cell array to input the matches into
bx_matches = cell(1,length(bx_genes));
for a = 1:length(bx_genes) 
    reg = strcat('\s',bx_genes(a),'\s');
    %find matches and for genes listed more than once, take only one match
    bx_matches(1,a) = regexp(gene_names,reg,'match','once');
    
end

%remove empty cell contents
bx_fc_matches = bx_fc(~cellfun('isempty',bx_matches));
bx_matches = bx_matches(~cellfun('isempty',bx_matches))';

%% Capan1 matching

%initialize an empty cell array to input the matches into
cp_matches = cell(1,length(cp_genes));
for a = 1:length(cp_genes) 
    reg = strcat('\s',cp_genes(a),'\s');
    %find matches and for genes listed more than once, take only one match
    cp_matches(1,a) = regexp(gene_names,reg,'match','once');
    
end

%remove empty cell contents
cp_fc_matches = cp_fc(~cellfun('isempty',cp_matches));
cp_matches = cp_matches(~cellfun('isempty',cp_matches))';

%% PANC1 matching

%initialize an empty cell array to input the matches into
pn_matches = cell(1,length(pn_genes));
for a = 1:length(pn_genes) 
    
    reg = strcat('\s',pn_genes(a),'\s');
    %find matches and for genes listed more than once, take only one match
    pn_matches(1,a) = regexp(gene_names,reg,'match','once');
    
end

%remove empty cell contents
pn_fc_matches = pn_fc(~cellfun('isempty',pn_matches));
pn_matches = pn_matches(~cellfun('isempty',pn_matches))';

%% TU8902 matching

%initialize an empty cell array to input the matches into
tu_matches = cell(1,length(tu_genes));
for a = 1:length(tu_genes) 
    
    reg = strcat('\s',tu_genes(a),'\s');
    %find matches and for genes listed more than once, take only one match
    tu_matches(1,a) = regexp(gene_names,reg,'match','once');
    
end

%remove empty cell contents
tu_fc_matches = tu_fc(~cellfun('isempty',tu_matches));
tu_matches = tu_matches(~cellfun('isempty',tu_matches))';

%% TU8988T matching

%initialize an empty cell array to input the matches into
tut_matches = cell(1,length(tut_genes));
for a = 1:length(tut_genes) 
    
    reg = strcat('\s',tut_genes(a),'\s');
    %find matches and for genes listed more than once, take only one match
    tut_matches(1,a) = regexp(gene_names,reg,'match','once');
    
end

%remove empty cell contents
tut_fc_matches = tut_fc(~cellfun('isempty',tut_matches));
tut_matches = tut_matches(~cellfun('isempty',tut_matches))';

%% UM2 matching

%initialize an empty cell array to input the matches into
um2_matches = cell(1,length(um2_genes));
for a = 1:length(um2_genes) 
    
    reg = strcat('\s',um2_genes(a),'\s');
    %find matches and for genes listed more than once, take only one match
    um2_matches(1,a) = regexp(gene_names,reg,'match','once');
    
end

%remove empty cell contents
um2_fc_matches = um2_fc(~cellfun('isempty',um2_matches));
um2_matches = um2_matches(~cellfun('isempty',um2_matches))';

%% UM90 matching

%initialize an empty cell array to input the matches into
um90_matches = cell(1,length(um90_genes));
for a = 1:length(um90_genes) 
    
    reg = strcat('\s',um90_genes(a),'\s');
    %find matches and for genes listed more than once, take only one match
    um90_matches(1,a) = regexp(gene_names,reg,'match','once');
    
end

%remove empty cell contents
um90_fc_matches = um90_fc(~cellfun('isempty',um90_matches));
um90_matches = um90_matches(~cellfun('isempty',um90_matches))';

%% Averaging repeated genes fold change values

all_matches = [bx_matches; cp_matches; pn_matches; tu_matches; tut_matches; um2_matches; um90_matches];
all_fc = [bx_fc_matches; cp_fc_matches; pn_fc_matches; tu_fc_matches; tut_fc_matches; um2_fc_matches; um90_fc_matches];

%alphebetize matches
[all_matches_sorted,i] = sort(all_matches);

%sort fc in same order
all_fc_sorted = all_fc(i);


all_matches_joined = join(all_matches);
all_matches_joined = strcat({' '}, all_matches_joined,{' '});
unique_matches = unique(all_matches);

duplicates = cell(length(unique_matches),1);
bool = zeros(length(unique_matches),1);
num_rep = zeros(length(unique_matches),1);

for b = 1:length(unique_matches)
   
    reg = strcat('\s',unique_matches(b),'\s');
    duplicates(b,1) = regexp(all_matches_joined,reg,'match');
    
    %locate where there are duplicates
    bool(b,1) = (length(duplicates{b,1}) > 1);
    num_rep(b,1) = length(duplicates{b,1});
    
end

index1 = 1;
index2 = 0;
fc_averaged = zeros(length(unique_matches),1);
for c = 1:length(unique_matches)
    
    index2 = index2 + num_rep(c,1);
    fc_average = sum(all_fc_sorted(index1:index2,1))./num_rep(c,1);
    fc_averaged(c,1) = fc_average;
    
    index1 = index2 + 1;
end

genes_and_fcs = table(unique_matches, fc_averaged);

