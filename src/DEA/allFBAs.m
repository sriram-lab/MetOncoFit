%% load in RECON1 model
load(fullfile('..','models','model_human_duarte.mat'));

%% Lyssiotis

% load in file containing match data
load(fullfile('RECON1_matches.mat'));

% starting with genes_and_fcs from mapping file
LR1_panGenes = genes_and_fcs.unique_matches;
LR1_panFcs = genes_and_fcs.fc_averaged;

LR1_clGenes = matches(1,:);
LR1_clFcs = matches(2,:);

% use FBA function 
% pan: A struct with fields 'Model', 'Solution', and 'KO_solution'.
% cl: A struct with fields 'Models', 'Solutions', and 'KO_solutions'.
[LR1_pan,LR1_cl] = FBA(metabolicmodel,LR1_panGenes,LR1_panFcs, LR1_clGenes, LR1_clFcs);

%% CCLE 

load(fullfile('mat_files','CCLE_mapping_dat_5_29.mat'));

% use FBA function 
% pan: A struct with fields 'Model', 'Solution', and 'KO_solution'.
% cl: A struct with fields 'Models', 'Solutions', and 'KO_solutions'.
[pan_model, cl_models] = FBA(metabolicmodel, pan_data.unique_names, ...
    pan_data.values_averaged, matches(1,:), matches(2,:));

%% TCGA

load(fullfile('mat_files','TCGA_mapping_6_01.mat'));

%starting with genes_and_fcs from mapping file
R1_panGenes = genes_and_fcs.unique_matches;
R1_panFcs = genes_and_fcs.fc_averaged;

R1_clGenes = matches(1,:);
R1_clFcs = matches(2,:);

% use FBA function 
% pan: A struct with fields 'Model', 'Solution', and 'KO_solution'.
% cl: A struct with fields 'Models', 'Solutions', and 'KO_solutions'.
[R1_pan,R1_cl] = FBA(metabolicmodel,R1_panGenes,R1_panFcs, R1_clGenes, R1_clFcs);