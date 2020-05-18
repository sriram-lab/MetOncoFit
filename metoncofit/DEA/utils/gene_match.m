function [ data_table ] = gene_match( model_genes, cell_line_genes, fold_change_data )

%combine the contents of multiple cells into one string
gene_names = join(model_genes);

%initialize an empty cell array to input the matches into
matching_genes = cell(1,length(cell_line_genes));
for a = 1:length(cell_line_genes) 
    
    reg = strcat('\s',cell_line_genes(a),'\s');
    %find matches and for genes listed more than once, take only one match
    matching_genes(1,a) = regexp(gene_names,reg,'match','once');
    
end

%remove empty cell contents
fold_change_data = fold_change_data(~cellfun('isempty',matching_genes));
matching_genes = matching_genes(~cellfun('isempty',matching_genes));
data_table = table(matching_genes,fold_change_data);
