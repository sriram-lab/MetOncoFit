%% averaging repeated km values
[~,km_ec] = xlsread('km.xlsx','A2:A20886');
km = xlsread('km.xlsx','B:B');

%want only positive kinetics parameter values
km_ec = km_ec(km>0);
km = km(km>0);

%to export non negative values for data visualization in R
km_tab = table(km_ec, km);
writetable(km_tab, 'brenda_km.xlsx','FileType', 'spreadsheet');

ecs_kms = average_repeats(km_ec, km);
ecs_kms = table2cell(ecs_kms);

%% averaging repeated kcat values

[~,kcat_ec] = xlsread('kcat.xlsx','G2:G10381');
kcat = xlsread('kcat.xlsx','F:F');

%keep only postive values
kcat_ec = kcat_ec(kcat>0);
kcat = kcat(kcat>0);

%to export non negative values for data visualization in R
kcat_tab = table(kcat_ec, kcat);
writetable(kcat_tab, 'brenda_kcat.xlsx','FileType', 'spreadsheet');

ecs_kcats = average_repeats(kcat_ec, kcat);
ecs_kcats = table2cell(ecs_kcats);

%% ec number to gene symbol

%variables of conversion data are labeled with _c for clarity
names = readcell('geneNamesAndSymbols.xlsx');
ecs_c = names(:,3);
symbols_c = names(:,1);

%remove <missing> from ecs 
symbols_c = symbols_c(cellfun(@ischar,ecs_c));
ecs_c = ecs_c(cellfun(@ischar,ecs_c));

%remove commas from string
reg = ',';

for a = 1:length(ecs_c)
    
    index = regexp(ecs_c{a}, reg);
    ecs_c{a}(index) = ' ';
    
end

%load ecs from BRENDA data
[~,km_ec] = xlsread('km.xlsx','A2:A20886');
km_ec = join(km_ec);
[~,kcat_ec] = xlsread('kcat.xlsx','G2:G10381');
kcat_ec = join(kcat_ec);

% example reg = '(1\.1\.1\.1 |1\.1\.1\.1$| 1\.1\.1\.1$)' 
%to account for being the first entry, last entry, or only entry
%because some gene symbols map to multiple ec numbers

%find all . and add a \ before them
%so the regular expression will not use . as any character
reg = '\.';
for c = 1:length(ecs_c)
    
    ecs_c{c} = regexprep(ecs_c{c},reg,'\\.');
    
end

km_ecs = cell(1, length(ecs_c));
kcat_ecs = cell(1, length(ecs_c));
for b = 1:length(ecs_c)
    
    %make expression as explained in comment above
    reg = strcat(ecs_c{b}, {' |'},ecs_c{b}, {'$| '}, ecs_c{b}, '$');
    km_ecs(1,b) = regexp(km_ec, reg, 'match','once');
    kcat_ecs(1,b) = regexp(kcat_ec, reg, 'match','once');
    
end

% some ec numbers map to multiple symbols, I kept all possible symbols to
% increase the chance of mapping to a gene symbol in the model
km_symbols = symbols_c(~cellfun('isempty',km_ecs));
km_ecs = km_ecs(~cellfun('isempty',km_ecs))';
km_sym_ec = [km_symbols, km_ecs, cell(length(km_ecs),1)];

kcat_symbols = symbols_c(~cellfun('isempty',kcat_ecs));
kcat_ecs = kcat_ecs(~cellfun('isempty',kcat_ecs))';
kcat_sym_ec = [kcat_symbols, kcat_ecs, cell(length(kcat_ecs),1)];

%% matching kinetics parameters to symbols
%% km

%find all . and add a \ before them
reg = '\.';
for c = 1:length(ecs_kms)
    
    ecs_kms{c,1} = regexprep(ecs_kms{c,1},reg,'\\.');
    
end

for d = 1:length(ecs_kms)
    for e = 1:length(km_sym_ec)
        reg = ecs_kms{d,1};
        match = regexp(km_sym_ec{e,2}, reg);
        
        if ~isempty(match)
            km_sym_ec{e,3} = ecs_kms{d,2};
        end
        
    end
end


%% kcat
%find all . and add a \ before them

reg = '\.';
for c = 1:length(ecs_kcats)
    
    ecs_kcats{c,1} = regexprep(ecs_kcats{c,1},reg,'\\.');
    
end

for d = 1:length(ecs_kcats)
    for e = 1:length(kcat_sym_ec)
        reg = ecs_kcats{d,1};
        match = regexp(kcat_sym_ec{e,2}, reg);
        
        if ~isempty(match)
            kcat_sym_ec{e,3} = ecs_kcats{d,2};
        end
        
    end
end
