% EssentialGenes.m
clear; clc;


% Input the core model
model1=readSBML('./NCI60_SBML_models/786-0.xml',1000);
model2=readSBML('./NCI60_SBML_models/A498.xml',1000);
model3=readSBML('./NCI60_SBML_models/ACHN.xml',1000);
model4=readSBML('./NCI60_SBML_models/CAKI-1.xml',1000);
model5=readSBML('./NCI60_SBML_models/SN12C.xml',1000);
model6=readSBML('./NCI60_SBML_models/TK10.xml',1000);
model7=readSBML('./NCI60_SBML_models/UO-31.xml',1000);


models = [model1, model2, model3, model4, model5, model6, model7];
names = ['786-0 ','A498  ','ACHN  ','CAKI-1','SN12C ','TK10  ','UO-31 '];

tol = 1e-6; % Growth rate lower limit

%Create list of gene IDs
genelist=model1.genes;
genelistrounded=(floor(str2double(genelist)));
%writetable(table(unique(genelistrounded),'stable'),'Renal_GeneList.csv');

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
writetable(FullTable,'./Gene KO Tables/Renal Gene Knockout Table.csv')

