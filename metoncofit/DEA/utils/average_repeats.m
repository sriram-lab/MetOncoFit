function [ genes_and_fcs ] = average_repeats( gene, fold_changes )

%alphebetize matches
[~,i] = sort(gene); 

%sort fc in same order
fc_sorted = fold_changes(i); 


matches_joined = join(gene); 
matches_joined = strcat({' '}, matches_joined,{' '}); 
unique_matches = unique(gene);

duplicates = cell(length(unique_matches),1);
bool = zeros(length(unique_matches),1);
num_rep = zeros(length(unique_matches),1);

for b = 1:length(unique_matches)
   
    reg = strcat('\s',unique_matches(b),'\s');
    duplicates(b,1) = regexp(matches_joined,reg,'match'); 
    
    %locate where there are duplicates
    bool(b,1) = (length(duplicates{b,1}) > 1);
    num_rep(b,1) = length(duplicates{b,1});
    
end

index1 = 1;
index2 = 0;
fc_averaged = zeros(length(unique_matches),1);
for c = 1:length(unique_matches)
    
    index2 = index2 + num_rep(c,1);
    fc_average = sum(fc_sorted(index1:index2,1))./num_rep(c,1); 
    fc_averaged(c,1) = fc_average;
    
    index1 = index2 + 1;
end

genes_and_fcs = table(unique_matches, fc_averaged);
