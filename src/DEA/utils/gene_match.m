function [ data_table ] = gene_match( model_genes, cell_line_genes, fold_change_data )

%{ 
This function maps genes and fold change values to the genes from a
metabolic model

Inputs: 
model_genes: a cell array containing genes from the desired model
cell_line_genes: a cell array containing genes to map to the model
fold_change_data: a vector containing the fold change values mapped to 
cell_line genes

Outputs: 
data_table: a table that contains the gene names and their fold change
values
%}

%combine the contents of multiple cells into one string
gene = join(model_genes);

%initialize an empty cell array to input the matches into
matching_genes = cell(1,length(cell_line_genes));
for a = 1:length(cell_line_genes) 
    
    %create a regular expression to find if there is a match within the
    %genes string
    %the spaces added to the beginning and end ensure a substring will
    %not be chosen as a match (ex. 'earth' as a match for 'art')
    reg = strcat('\s',cell_line_genes(a),'\s');
    %find matches and for genes listed more than once, take only one match
    matching_genes(1,a) = regexp(gene,reg,'match','once');
    
end

%remove empty cell contents
fold_change_data = fold_change_data(~cellfun('isempty',matching_genes));
matching_genes = matching_genes(~cellfun('isempty',matching_genes));
data_table = table(matching_genes,fold_change_data);