%% averaging repeated km values
[~,km_ec] = xlsread('km.xlsx','A2:A20886');
km = xlsread('km.xlsx','B:B');

%want only positive kinetics parameter values
km_ec = km_ec(km>0);
km = km(km>0);

%sort km in the order unique will sort ec
[km_ec,i] = sort(km_ec);
km = km(i);

% join into one string to search through
ec_joined = join(km_ec);

%get list of unique ec numbers
unique_ec = unique(km_ec);

%initializing empty cells and matricies to hold data
duplicates = cell(length(unique_ec),1);
num_rep = zeros(length(unique_ec),1);
    
for a = 1:length(unique_ec)
    
    reg = strcat('\s',unique_ec(a),'\s');
    duplicates(a,1) = regexp(ec_joined,reg,'match');
    
    %locate where there are duplicates
    num_rep(a,1) = length(duplicates{a,1});
    
end

index1 = 1;
index2 = 0;
    
%average km value if there are multiple entries
km_averaged = zeros(length(unique_ec),1);
for b = 1:length(unique_ec)

    index2 = index2 + num_rep(b,1);
    km_average = sum(km(index1:index2,1))./num_rep(b,1);
    km_averaged(b,1) = km_average;

    index1 = index2 + 1;
end

%store data in a cell array
km_averaged = num2cell(km_averaged);
ecs_kms = [unique_ec, km_averaged];

%% averaging repeated kcat values

[~,kcat_ec] = xlsread('kcat.xlsx','G2:G10381');
kcat = xlsread('kcat.xlsx','F:F');

%keep only postive values
kcat_ec = kcat_ec(kcat>0);
kcat = kcat(kcat>0);

%sort km in the order unique will sort ec
[kcat_ec,i] = sort(kcat_ec);
kcat = kcat(i);

% join into one string to search through
ec_joined = join(kcat_ec);

% get list of unique ec numbers
unique_ec = unique(kcat_ec);

%initializing empty cells and matricies to hold data
duplicates = cell(length(unique_ec),1);
num_rep = zeros(length(unique_ec),1);
    
for a = 1:length(unique_ec)
    
    reg = strcat('\s',unique_ec(a),'\s');
    duplicates(a,1) = regexp(ec_joined,reg,'match');
    
    %locate where there are duplicates
    num_rep(a,1) = length(duplicates{a,1});
    
end

index1 = 1;
index2 = 0;
    
kcat_averaged = zeros(length(unique_ec),1);
for b = 1:length(unique_ec)

    index2 = index2 + num_rep(b,1);
    kcat_average = sum(kcat(index1:index2,1))./num_rep(b,1);
    kcat_averaged(b,1) = kcat_average;

    index1 = index2 + 1;
end

kcat_averaged = num2cell(kcat_averaged);
unique_ec = strtrim(unique_ec);
ecs_kcats = [unique_ec, kcat_averaged];

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
