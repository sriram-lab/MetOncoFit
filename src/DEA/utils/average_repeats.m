function [ averaged_values ] = average_repeats( names, values )

%alphebetize matches
[~,i] = sort(names); 

%sort fc in same order
values_sorted = values(i); 


names_joined = join(names); 
names_joined = strcat({' '}, names_joined,{' '}); 
unique_names = unique(names);

duplicates = cell(length(unique_names),1);
bool = zeros(length(unique_names),1);
num_rep = zeros(length(unique_names),1);

for b = 1:length(unique_names)
   
    reg = strcat('\s',unique_names(b),'\s');
    duplicates(b,1) = regexp(names_joined,reg,'match'); 
    
    %locate where there are duplicates
    bool(b,1) = (length(duplicates{b,1}) > 1);
    num_rep(b,1) = length(duplicates{b,1});
    
end

index1 = 1;
index2 = 0;
values_averaged = zeros(length(unique_names),1);

for c = 1:length(unique_names)
    
    index2 = index2 + num_rep(c,1);
    value_average = sum(values_sorted(index1:index2,1))./num_rep(c,1); 
    values_averaged(c,1) = value_average;
    
    index1 = index2 + 1;
end

averaged_values = table(unique_names, values_averaged);