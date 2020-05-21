%% Import gene names and fold change values

sheets = {'Sheet 1', 'Sheet 2', 'Sheet 3', 'Sheet 4', 'Sheet 5', 'Sheet 6', 'Sheet 7'};

%bx_genes,cp_genes,pn_genes,tu_genes,tut_genes,um2_genes,um90_genes
genes = cell(1,length(sheets));

%bx_fc,cp_fc,pn_fc,tu_fc,tut_fc,um2_fc,um90_fc
fc = cell(1,length(sheets));

for sheet = 1:length(sheets)

    [~,genes{1, sheet},~] = xlsread(fullfile('..','Lyssiotis_data','genes_sig.xlsx'),sheets{sheet},'A:A');
    [fc{1, sheet},~ ,~] = xlsread(fullfile('..','Lyssiotis_data','fold_change.xlsx'), sheets{sheet}, 'A:A');
    
end 

%metabolic model
load(fullfile('..','models','model_human_duarte.mat'));
model_genes = metabolicmodel.genes;


%% matching

gene_str = join(model_genes);

%add space at beginning to use as an indicator of a new word
gene_str = strcat({' '}, gene_str, {' '});

%initialize an empty cell array to input gene names and fold changes
%each cell will hold an array: 1st column = names, 2nd column = fold changes
matches = cell(2, length(genes));
%place to store matches from all cell lines for next step

for a = 1:length(matches)
    
    match = cell(length(genes{a}), 1);
    fcs = fc{a};
    
    for b = 1:length(genes{a}) 
        
        reg = strcat('\s',genes{a}{b},'\s');
        %find matches and for genes listed more than once, take only one match
        match(b, 1) = regexp(gene_str,reg,'match','once');
        
    end
    
    %remove empty cell contents
    fcs = fcs(~cellfun('isempty', match));
    match = match(~cellfun('isempty', match));
    
    matches{1,a} = match;
    matches{2,a} = fcs;
end

%% Averaging repeated genes fold change values

all_matches = [];
all_fc = [];

for c = 1:length(matches)
    
    all_matches = [all_matches; matches{1,c}];
    all_fc = [all_fc; matches{2,c}];
    
end

%alphebetize matches
[all_matches_sorted,i] = sort(all_matches);

%sort fc in same order
all_fc_sorted = all_fc(i);

all_matches_joined = join(all_matches);

%add spaces around all words, ensures regex can match words at beginning
%and end of string 
all_matches_joined = strcat({' '}, all_matches_joined,{' '});

unique_matches = unique(all_matches);

%initializing empty cells and matricies to hold data
duplicates = cell(length(unique_matches),1);
bool = zeros(length(unique_matches),1);
num_rep = zeros(length(unique_matches),1);

for d = 1:length(unique_matches)
   
    reg = strcat('\s',unique_matches(d),'\s');
    duplicates(d,1) = regexp(all_matches_joined,reg,'match');
    
    %locate where there are duplicates
    bool(d,1) = (length(duplicates{d,1}) > 1);
    num_rep(d,1) = length(duplicates{d,1});
    
end

index1 = 1;
index2 = 0;
fc_averaged = zeros(length(unique_matches),1);
for e = 1:length(unique_matches)
    
    index2 = index2 + num_rep(e,1);
    fc_average = sum(all_fc_sorted(index1:index2,1))./num_rep(e,1);
    fc_averaged(e,1) = fc_average;
    
    index1 = index2 + 1;
end

genes_and_fcs = table(unique_matches, fc_averaged);

