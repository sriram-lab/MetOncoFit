%need to initliaze Cobra and set Gurobi as solver before running

% EssentialGenes.m
clear; clc;


% Input the core models
model1=readSBML('./NCI60_SBML_models/BT-549.xml',1000);
model2=readSBML('./NCI60_SBML_models/T47D.xml',1000);
model3=readSBML('./NCI60_SBML_models/MCF7.xml',1000);
model4=readSBML('./NCI60_SBML_models/MDA-MB-231.xml',1000);
model5=readSBML('./NCI60_SBML_models/HS-578-T.xml',1000);

% Create list of models and cell lines
% All cell line names must be same length in order to create table

models = [model1, model2, model3, model4, model5];
names = ['BT-549    ','T47D      ','MCF7      ','MDA-MB-231','HS-578-T  '];

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
    CL(:) = {names(n:n+9)};
    
    % create gene ko table with list of genes, grRatio and cell line
    ko = table(genelistrounded,grRatio,CL);
    % get rid of duplicate genes
    ko_unq = unique(ko,'stable');
    
    %add data from each cell line to table
    FullTable=[FullTable;ko_unq];
    n=n+10;
end

%write entire table to csv file
writetable(FullTable,'./Gene KO Tables/Breast Gene Knockout Table.csv')


