% EssentialGenes.m
clear; clc;


% Input the core model
model1=readSBML('./NCI60_SBML_models/SF268.xml',1000);
model2=readSBML('./NCI60_SBML_models/SF539.xml',1000);
model3=readSBML('./NCI60_SBML_models/SNB-19.xml',1000);
model4=readSBML('./NCI60_SBML_models/SNB-75.xml',1000);
model5=readSBML('./NCI60_SBML_models/U251.xml',1000);

% Create list of models and cell lines
% All cell line names must be same length in order to create table

models = [model1, model2, model3, model4, model5];
names = ['SF268 ','SF539 ','SNB-19','SNB-75','U251  '];

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
writetable(FullTable,'./Gene KO Tables/CNS Gene Knockout Table.csv')


