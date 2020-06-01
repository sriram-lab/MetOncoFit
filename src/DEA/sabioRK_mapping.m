%% loading data and filtering out unwanted values/missing data
%for legibility, data from the sabio data file are labeled with _s
sabio_data = readcell('homosapien.xlsx');
param_type_s = sabio_data(2:end,10);
uniprotIDs_s = sabio_data(2:end,8);
param_val_s = sabio_data(2:end,12);
enzyme_type_s = sabio_data(2:end,7);
param_unit_s = sabio_data(2:end,15);

%boolean vectors to find the data that is wanted (done in multiple lines to
%improve legibility
has_id = cellfun(@ischar,uniprotIDs_s);
is_wt = contains(enzyme_type_s, 'wildtype');
is_kcat = strcmp('kcat',param_type_s);
is_km = strcmp('Km',param_type_s);
is_vmax = strcmp('Vmax',param_type_s);
has_unit = ~strcmp('-',param_unit_s);
has_param_value = cellfun(@isnumeric,param_val_s);

%combining data that is desired, must have id, be wt, have a unit and be either a
%kcat, km, or vmax parameter type
keep = has_param_value & has_id & is_wt & has_unit & (is_kcat | is_km | is_vmax);

%remove parameters that are not wanted
param_type_s = param_type_s(keep);
uniprotIDs_s = uniprotIDs_s(keep);
param_val_s = param_val_s(keep);
enzyme_type_s = enzyme_type_s(keep);
param_unit_s = param_unit_s(keep);

%% analyzing and converting units 
%dif_units = unique(param_unit_s);
%{'M';'M*s^(-1)';'M*s^(-1)*M^(-1)';'M*s^(-1)*g^(-1)';'M^2';'g/(s*g)';'katal*g^(-1)';'katal*l^(-1)';'katal*s^(-1)';'l*g^(-1)*s^(-1)';'mg/ml';'mol';'mol*s^(-1)*g^(-1)';'mol*s^(-1)*mol^(-1)';'mol/mol';'mol/s';'s^(-1)';'s^(-1)*g^(-1)'}

%katal = mol/s

%units to keep:
%for km: 'M' --> convert to mM to match BRENDA
%for kcat: 's^(-1)', 'g/(s*g)', 'mol*s^(-1)*mol^(-1)', 'M*s^(-1)*M^(-1)'
%for vmax: 'katal*g^(-1)', 'mol*s^(-1)*g^(-1)'
units_keep = {'s^(-1)', 'g/(s*g)', 'mol*s^(-1)*mol^(-1)', ...
    'M*s^(-1)*M^(-1)', 'katal*g^(-1)', 'mol*s^(-1)*g^(-1)'};

unit_bool = zeros(length(param_unit_s),1);
%cannot remove param rows as we go must remove at the end
for a = 1:length(param_unit_s)
    if strcmp('M',param_unit_s{a})
        %convert M to mM
        param_val_s{a} = param_val_s{a}*1000;
        unit_bool(a) = 1;
    else
        for b = 1:length(units_keep)
            if strcmp(units_keep{b}, param_unit_s{a})
                unit_bool(a) = 1;
            end
        end
    end        
end

param_unit_s = param_unit_s(logical(unit_bool));
param_type_s = param_type_s(logical(unit_bool));
uniprotIDs_s = uniprotIDs_s(logical(unit_bool));
param_val_s = param_val_s(logical(unit_bool));
enzyme_type_s = enzyme_type_s(logical(unit_bool));

%% 
%import list of uniprotIDs and gene symbols

%IDs from the ID_list and symbols from the symbol_list correspond

%for legibility variables that are from the uniprot to symbol conversion
%file are labeled with _c
uniprotID_conv = readcell('uniprotToSymbol.xlsx');
ID_list_c = uniprotID_conv(:,2);
symbol_list_c = uniprotID_conv(:,1);

%remove symbols that dont have an ID
symbol_list_c = symbol_list_c(cellfun(@ischar,ID_list_c));
ID_list_c = ID_list_c(cellfun(@ischar,ID_list_c));

%find matches and store locations of matches 
%preset locations matrix to arbitrarily large size
%first row is the location within the sabio data 
%second row is the location within the conversion data
locations = zeros(2,50000);
count = 1;
for a = 1:length(uniprotIDs_s)
    for b = 1:length(ID_list_c)
        match = regexp(ID_list_c{b}, uniprotIDs_s{a}, 'match');
        
        %store a list of locations to match symbols to parameter values
        if ~isempty(match)
            locations(1, count) = a;
            locations(2, count) = b;
            count = count + 1;
        end
        
    end
    
end

%first column = location within uniprotID, second column = location within
%ID_list
%remove empty/zero entries from large locations matrix
locs = [nonzeros(locations(1,:)), nonzeros(locations(2,:))];

%% organize data by parameter type
% data = col 1 = symbol col 2 = param type col 3 = param value

%using the locations to match symbol to parameter type and parameter value
km_data = cell(length(locs),3);
km_count = 1;
kcat_data = cell(length(locs),3);
kcat_count = 1;
vmax_data = cell(length(locs),3);
vmax_count = 1;
for c = 1:length(locs)
    
    if isequal(param_type_s{locs(c,1)},'Km')
        km_data{km_count, 1} = symbol_list_c{locs(c,2)};
        km_data{km_count, 2} = param_type_s{locs(c,1)};
        km_data{km_count, 3} = param_val_s{locs(c,1)};
        km_count = km_count + 1;
        
    elseif isequal(param_type_s{locs(c,1)},'kcat')
        kcat_data{kcat_count, 1} = symbol_list_c{locs(c,2)};
        kcat_data{kcat_count, 2} = param_type_s{locs(c,1)};
        kcat_data{kcat_count, 3} = param_val_s{locs(c,1)};
        kcat_count = kcat_count + 1;
        
    else
        vmax_data{vmax_count, 1} = symbol_list_c{locs(c,2)};
        vmax_data{vmax_count, 2} = param_type_s{locs(c,1)};
        vmax_data{vmax_count, 3} = param_val_s{locs(c,1)};
        vmax_count = vmax_count + 1;
        
    end
    
end

%remove empty cell contents and average values with more than one entry
%write data to excel files for visualization in R
%%
km_data = km_data(~all(cellfun('isempty', km_data), 2), :);
km_table = average_repeats(km_data(:,1),vertcat(km_data{:,3}));
writetable(km_table,'km_sabio_data.xlsx','FileType', 'spreadsheet');

kcat_data = kcat_data(~all(cellfun('isempty', kcat_data), 2), :);
kcat_table = average_repeats(kcat_data(:,1),vertcat(kcat_data{:,3}));
writetable(kcat_table,'kcat_sabio_data.xlsx','FileType', 'spreadsheet');

vmax_data = vmax_data(~all(cellfun('isempty', vmax_data), 2), :);
vmax_table = average_repeats(vmax_data(:,1),vertcat(vmax_data{:,3}));
writetable(vmax_table,'vmax_sabio_data.xlsx','FileType', 'spreadsheet');
