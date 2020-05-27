%% averaging repeated km values
[~,km_ec] = xlsread('km.xlsx','A2:A20886');
km = xlsread('km.xlsx','B:B');

km_ec = km_ec(km>0);
km = km(km>0);

%sort km in the order unique will sort ec
[km_ec,i] = sort(km_ec);
km = km(i);

ec_joined = join(km_ec);

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
    
km_averaged = zeros(length(unique_ec),1);
for b = 1:length(unique_ec)

    index2 = index2 + num_rep(b,1);
    km_average = sum(km(index1:index2,1))./num_rep(b,1);
    km_averaged(b,1) = km_average;

    index1 = index2 + 1;
end

km_averaged = num2cell(km_averaged);
ecs_kms = [unique_ec, km_averaged];

%% averaging repeated kcat values

[~,kcat_ec] = xlsread('kcat.xlsx','G2:G10381');
kcat = xlsread('kcat.xlsx','F:F');

kcat_ec = kcat_ec(kcat>0);
kcat = kcat(kcat>0);

%sort km in the order unique will sort ec
[kcat_ec,i] = sort(kcat_ec);
kcat = kcat(i);

ec_joined = join(kcat_ec);

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

names = readcell('geneNamesAndSymbols.xlsx');
ecs = names(:,3);
symbols = names(:,1);

%remove <missing> from ecs 
symbols = symbols(cellfun(@ischar,ecs));
ecs = ecs(cellfun(@ischar,ecs));

%remove commas 
reg = ',';

for a = 1:length(ecs)
    
    index = regexp(ecs{a}, reg);
    ecs{a}(index) = ' ';
    
end

%load ecs from BRENDA data
[~,km_ec] = xlsread('km.xlsx','A2:A20886');
km_ec = join(km_ec);
[~,kcat_ec] = xlsread('kcat.xlsx','G2:G10381');
kcat_ec = join(kcat_ec);

% example reg = '(1\.1\.1\.1 |1\.1\.1\.1$| 1\.1\.1\.1$)';

%find all . and add a \ before them
reg = '\.';
for c = 1:length(ecs)
    
    ecs{c} = regexprep(ecs{c},reg,'\\.');
    
end

km_ecs = cell(1, length(ecs));
kcat_ecs = cell(1, length(ecs));
for b = 1:length(ecs)
    
    reg = strcat(ecs{b}, {' |'},ecs{b}, {'$| '}, ecs{b}, '$');
    km_ecs(1,b) = regexp(km_ec, reg, 'match','once');
    kcat_ecs(1,b) = regexp(kcat_ec, reg, 'match','once');
    
end

% some ec numbers map to multiple symbols, I kept all possible symbols to
% increase the chance of mapping to a gene symbol in the model
km_symbols = symbols(~cellfun('isempty',km_ecs));
km_ecs = km_ecs(~cellfun('isempty',km_ecs))';
km_sym_ec = [km_symbols, km_ecs, cell(length(km_ecs),1)];

kcat_symbols = symbols(~cellfun('isempty',kcat_ecs));
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
