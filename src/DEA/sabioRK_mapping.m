%% loading data and filtering out unwanted values/missing data
sabio_data = readcell('homosapien.xlsx');
param_type = sabio_data(2:end,10);
uniprotIDs = sabio_data(2:end,8);
param_val = sabio_data(2:end,12);
enzyme_type = sabio_data(2:end,7);

%boolean vectors to find the data that is wanted (done in multiple lines to
%improve legibility
has_id = cellfun(@ischar,uniprotIDs);
is_wt = contains(enzyme_type, 'wildtype');
is_kcat = strcmp('kcat',param_type);
is_km = strcmp('Km',param_type);
is_vmax = strcmp('Vmax',param_type);

%combining data that is desired
keep = has_id & is_wt & (is_kcat | is_km | is_vmax);

%remove parameters that are not wanted
param_type = param_type(keep);
uniprotIDs = uniprotIDs(keep);
param_val = param_val(keep);
enzyme_type = enzyme_type(keep);


%% 
%import list of uniprotIDs and gene symbols
uniprotID_conv = readcell('uniprotToSymbol.xlsx');
ID_list = uniprotID_conv(:,2);
symbol_list = uniprotID_conv(:,1);

%remove symbols that dont have an ID
symbol_list = symbol_list(cellfun(@ischar,ID_list));
ID_list = ID_list(cellfun(@ischar,ID_list));

%find matches and store locations of matches 
%preset locations matrix to arbitrarily large size
locations = zeros(2,50000);
count = 1;
for a = 1:length(uniprotIDs)
    for b = 1:length(ID_list)
        match = regexp(ID_list{b}, uniprotIDs{a}, 'match');
        if ~isempty(match)
            locations(1, count) = a;
            locations(2, count) = b;
            count = count + 1;
        end
        
    end
    
end

%first column = location within uniprotID, second column = location within
%ID_list
locs = [nonzeros(locations(1,:)), nonzeros(locations(2,:))];

%% organize data by parameter type
% data = col 1 = symbol col 2 = param type col 3 = param value
km_data = cell(length(locs),3,3);
km_count = 1;
kcat_data = cell(length(locs),3);
kcat_count = 1;
vmax_data = cell(length(locs),3);
vmax_count = 1;
for c = 1:length(locs)
    
    if isequal(param_type{locs(c,1)},'Km')
        km_data{km_count, 1} = symbol_list{locs(c,2)};
        km_data{km_count, 2} = param_type{locs(c,1)};
        km_data{km_count, 3} = param_val{locs(c,1)};
        km_count = km_count + 1;
        
    elseif isequal(param_type{locs(c,1)},'kcat')
        kcat_data{kcat_count, 1} = symbol_list{locs(c,2)};
        kcat_data{kcat_count, 2} = param_type{locs(c,1)};
        kcat_data{kcat_count, 3} = param_val{locs(c,1)};
        kcat_count = kcat_count + 1;
        
    else
        vmax_data{vmax_count, 1} = symbol_list{locs(c,2)};
        vmax_data{vmax_count, 2} = param_type{locs(c,1)};
        vmax_data{vmax_count, 3} = param_val{locs(c,1)};
        vmax_count = vmax_count + 1;
        
    end
    
end

%remove empty cell contents
km_data = km_data(~all(cellfun('isempty', km_data), 2), :);
kcat_data = kcat_data(~all(cellfun('isempty', kcat_data), 2), :);
vmax_data = vmax_data(~all(cellfun('isempty', vmax_data), 2), :);