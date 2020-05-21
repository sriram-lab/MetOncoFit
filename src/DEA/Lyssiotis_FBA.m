initCobraToolbox(false);

%load model
load(fullfile('..','models','model_human_duarte.mat'));

%% Pan Pancreatic model
%load genes and fcs
load('RECON1_matches.mat');

%starting with genes_and_fcs from mapping file
genes = genes_and_fcs.unique_matches;
fcs = genes_and_fcs.fc_averaged;

%find up and down regulated genes
down_reg = genes(fcs < 1);
up_reg = genes(fcs > 1);

%use CFR (constrained flux regulation)
fields = {'eps','kap','rho','mode','eps2','pfba','kap2'};
empty_cell = cell(length(fields),1);
hyperparams = cell2struct(empty_cell,fields);
[constrained_model_pan, solution_pan] = CFR(metabolicmodel, hyperparams, up_reg, down_reg);

%use knockOut
geneKO_pan = knockOut(constrained_model_pan, 'GeneKO');

%% Cell line models

%load genes and fcs
load('all_cell_line_matches.mat');

models_and_sols = cell(2, length(matches));
KO_sols = cell(1, length(matches));

for a = 1:length(matches)
    
    genes = matches{1,a};
    fcs = matches{2,a};

    %find up and down regulated genes
    down_reg = genes(fcs < 1);
    up_reg = genes(fcs > 1);

    [constrained_model, solution] = CFR(metabolicmodel, hyperparams, up_reg, down_reg);

    %use knockOut
    geneKO = knockOut(constrained_model, 'GeneKO');
    
    models_and_sols{1,a} = constrained_model;
    models_and_sols{2,a} = solution;
    KO_sols{1,a} = geneKO;
    
end
