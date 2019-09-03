% EssentialGenes.m
clear; clc;


% Input the core model
model1=readSBML('./NCI60_SBML_models/DU-145.xml',1000);
model2=readSBML('./NCI60_SBML_models/PC-3.xml',1000);

models = [model1, model2];
names = ['DU-145','PC-3  '];

tol = 1e-6; % Growth rate lower limit

%Create list of gene IDs
genelist=model1.genes;
genelistrounded=(floor(str2double(genelist)));
%writetable(table(unique(genelistrounded),'stable'),'CNS_GeneList.csv');

k = length(models);
FullTable = {};

n=1;
for i = 1:k
    % get growth rate ratio for gene KO
    grRatio = singleGeneDeletion(models(i));
    grRatio(isnan(grRatio))=0;
    gr_list(1,:)=grRatio;
    
    % create list of cell lines equal to length of genes    
    CL    = cell(length(genelistrounded),1);
    % cell line name is from n to n+x
    CL(:) = {names(n:n+5)};
        
    % create gene ko table with list of genes, grRatio and cell line
    ko = table(genelistrounded,grRatio,CL);
    % get rid of duplicate genes
    ko_unq = unique(ko,'stable');
    
    %add data from each cell line to table
    FullTable=[FullTable;ko_unq];
    n=n+6;
end

%write entire table to csv file
writetable(FullTable,'./Gene KO Tables/Prostate Gene Knockout Table.csv')

