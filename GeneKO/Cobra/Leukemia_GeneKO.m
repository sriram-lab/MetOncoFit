% EssentialGenes.m
clear; clc;


% Input the core model
model1=readSBML('./NCI60_SBML_models/CCRF-CEM.xml',1000);
model2=readSBML('./NCI60_SBML_models/HL-60.xml',1000);
model3=readSBML('./NCI60_SBML_models/K-562.xml',1000);
model4=readSBML('./NCI60_SBML_models/MOLT-4.xml',1000);
model5=readSBML('./NCI60_SBML_models/RPMI-8226.xml',1000);
model6=readSBML('./NCI60_SBML_models/SR.xml',1000);

models = [model1, model2, model3, model4, model5, model6];
names = ['CCRF-CEM ','HL-60    ','K-562    ','MOLT-4   ','RPMI-8226','SR       '];

tol = 1e-6; % Growth rate lower limit

%Create list of gene IDs
genelist=model1.genes;
genelistrounded=(floor(str2double(genelist)));
%writetable(table(genelistrounded),'LeukemiaGeneList.csv');

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
    CL(:) = {names(n:n+8)};
        
    % create gene ko table with list of genes, grRatio and cell line
    ko = table(genelistrounded,grRatio,CL);
    % get rid of duplicate genes
    ko_unq = unique(ko,'stable');
    
    %add data from each cell line to table
    FullTable=[FullTable;ko_unq];
    n=n+9;
end

%write entire table to csv file
writetable(FullTable,'./Gene KO Tables/Leukemia Gene Knockout Table.csv')



