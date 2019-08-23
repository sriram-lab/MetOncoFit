% EssentialGenes.m
clear; clc;


% Input the core model
model_LOXIMVI=readSBML('./NCI60_SBML_models/LOXIMVI.xml',1000);
model_MALME_3M=readSBML('./NCI60_SBML_models/MALME-3M.xml',1000);
model_MDA_MB_435=readSBML('./NCI60_SBML_models/MDA-MB-435.xml',1000);
model_SK_MEL_2=readSBML('./NCI60_SBML_models/SK-MEL-2.xml',1000);
model_SK_MEL_28=readSBML('./NCI60_SBML_models/SK-MEL-28.xml',1000);
model_SK_MEL_5=readSBML('./NCI60_SBML_models/SK-MEL-5.xml',1000);
model_UACC_257=readSBML('./NCI60_SBML_models/UACC-257.xml',1000);
model_UACC_62=readSBML('./NCI60_SBML_models/UACC-62.xml',1000);

models = [model_LOXIMVI, model_MALME_3M, model_MDA_MB_435, model_SK_MEL_2, model_SK_MEL_28, model_SK_MEL_5, model_UACC_257, model_UACC_62];
names = ['LOXIMVI   ','MALME-3M  ','MDA-MB-435','SK-MEL-2  ','SK-MEL-28 ','SK-MEL-5  ','UACC-257  ','UACC-62   '];

tol = 1e-6; % Growth rate lower limit

%Create list of gene IDs
genelist=model_LOXIMVI.genes;
genelistrounded=(floor(str2double(genelist)));
%writetable(table(genelistrounded),'MelanomaGeneList.csv');

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
writetable(FullTable,'./Gene KO Tables/Melanoma Gene Knockout Table.csv')

