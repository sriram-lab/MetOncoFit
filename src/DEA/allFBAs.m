%% Lyssiotis
load(fullfile('..','models','model_human_duarte.mat'));
load(fullfile('..','models','Recon2.v05.mat'));
load(fullfile('..','models','Recon3D.mat'));
load(fullfile('..','models','HumanGEM.mat'));

% Recon1
load(fullfile('RECON1_matches.mat'));

%starting with genes_and_fcs from mapping file
LR1_panGenes = genes_and_fcs.unique_matches;
LR1_panFcs = genes_and_fcs.fc_averaged;

LR1_clGenes = matches(1,:);
LR1_clFcs = matches(2,:);

[LR1_pan,LR1_cl] = FBA(metabolicmodel,LR1_panGenes,LR1_panFcs, LR1_clGenes, LR1_clFcs);

% Recon2

load(fullfile('RECON2_matches.mat'));

%starting with genes_and_fcs from mapping file
LR2_panGenes = genes_and_fcs.unique_matches;
LR2_panFcs = genes_and_fcs.fc_averaged;

LR2_clGenes = matches(1,:);
LR2_clFcs = matches(2,:);

[LR2_pan,LR2_cl] = FBA(modelR205,LR2_panGenes,LR2_panFcs, LR2_clGenes, LR2_clFcs);

% Recon3

load(fullfile('RECON3D_matches.mat'));

%starting with genes_and_fcs from mapping file
LR3_panGenes = genes_and_fcs.unique_matches;
LR3_panFcs = genes_and_fcs.fc_averaged;

LR3_clGenes = matches(1,:);
LR3_clFcs = matches(2,:);

[LR3_pan,LR3_cl] = FBA(Recon3D,LR3_panGenes,LR3_panFcs, LR3_clGenes, LR3_clFcs);

% Human1

load(fullfile('Human1_matches.mat'));

%starting with genes_and_fcs from mapping file
LH_panGenes = genes_and_fcs.unique_matches;
LH_panFcs = genes_and_fcs.fc_averaged;

LH_clGenes = matches(1,:);
LH_clFcs = matches(2,:);

[LH_pan,LH_cl] = FBA(ihuman,LH_panGenes,LH_panFcs, LH_clGenes, LH_clFcs);

%% CCLE 

load(fullfile('..','models','model_human_duarte.mat'));

load(fullfile('mat_files','CCLE_mapping_dat_5_29.mat'));

[pan_model, cl_models] = FBA(metabolicmodel, pan_data.unique_names, ...
    pan_data.values_averaged, matches(1,:), matches(2,:));

%% TCGA

load(fullfile('..','models','model_human_duarte.mat'));

load(fullfile('mat_files','TCGA_mapping_6_01.mat'));

%starting with genes_and_fcs from mapping file
R1_panGenes = genes_and_fcs.unique_matches;
R1_panFcs = genes_and_fcs.fc_averaged;

R1_clGenes = matches(1,:);
R1_clFcs = matches(2,:);

[R1_pan,R1_cl] = FBA(metabolicmodel,R1_panGenes,R1_panFcs, R1_clGenes, R1_clFcs);