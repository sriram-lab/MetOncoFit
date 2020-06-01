% %each column holds a different file's data as structs
% TCGA_mut_data = cell(1,32);
% TCGA_wt_data = cell(1,6);
% 
% %make a cell array corresponding to all of the file names
% file_mut = cell(32,1);
% for a = 1:length(file_mut)
%     
%     file_name = strcat('mut_TCGA', num2str(a),'.txt');
%     file_mut{a,1} = file_name;
%     
% end
% 
% file_wt = cell(6,1);
% for a = 1:length(file_wt)
%     
%     file_name = strcat('wt_TCGA', num2str(a),'.txt');
%     file_wt{a,1} = file_name;
%     
% end
% 
file = vertcat(file_wt, file_mut);
% 
% %import all files as structs and store in a cell array
% for b = 1:length(file_mut)
%     
%     TCGA_mut_data{1,b} = importdata(fullfile('TCGA_data','mut_samples_TCGA',file_mut{b,1}));
%     
% end
% 
% %import all files as structs and store in a cell array
% for b = 1:length(file_wt)
%     
%     TCGA_wt_data{1,b} = importdata(fullfile('TCGA_data','normal_samples_TCGA',file_wt{b,1}));
%     
% end

%combine wt and mut data - first 6 = wt
TCGA_data = [TCGA_wt_data, TCGA_mut_data];

%import ens ID to symbol data and remove missing values
ens_to_sym = readcell(fullfile('gene_mapping','gene_name_info','ensgID_symbol.txt'));
ens_to_sym = ens_to_sym(all(cellfun(@ischar, ens_to_sym), 2), :);

for c = 1:length(TCGA_data)
    
    curr_data = TCGA_data{c}.data;
    TPM = FPKMtoTPM(curr_data);
    TCGA_data{c}.TPMs = TPM;
     
    %remove .# ending after ENSG ID
    curr_ID = TCGA_data{c}.textdata;
    for d = 1:length(curr_ID)
        reg = '\.\d+$';
        curr_ID{d} = regexprep(curr_ID{d},reg,'');
    end
     
    [IDs, symbols, bool] = gene_match(curr_ID, ens_to_sym(:,2), ens_to_sym(:,1));
    
    %add struct fields for TPMs, ensID, and symbols, remove unwanted struct
    %fields from import
    TCGA_data{c}.TPMs = TCGA_data{c}.TPMs(bool);
    TCGA_data{c}.ensID = IDs';
    TCGA_data{c}.symbol = symbols;
    TCGA_data{c} = rmfield(TCGA_data{c}, {'data','textdata','rowheaders'});
  
    writetable(struct2table(TCGA_data{c}), 'TCGA_TPM_data.xlsx','Sheet',file{c});
    
end

