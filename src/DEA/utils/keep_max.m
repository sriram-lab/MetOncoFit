function [ max_values ] = keep_max( names, values )

%{ 
This function takes in mapped names and values that may contain duplicate
names where you want to keep the max value out of all the values the same 
name is mapped to (rather than averaging, etc) and outputs a table with 
unique names mapped to its respective max value

ex. 
input:
a 1
a 2
b 3
c 2
c 1

output:
a 2
b 3 
c 2

Inputs: 
names: a cell array containing names that are mapped to values
values: a vector containing values that are mapped to names

Outputs: 
max_values: a table that contains unique names and their respective
max values
%}

    %alphebetize matches
    [~,i] = sort(names); 

    %sort fc in same order
    values_sorted = values(i); 

    %combine all names into one string for searching
    names_joined = join(names); 
    names_joined = strcat({' '}, names_joined,{' '}); 
    unique_names = unique(names);

    %create empty vectors to store data during for loop
    duplicates = cell(length(unique_names),1);
    bool = zeros(length(unique_names),1);
    num_rep = zeros(length(unique_names),1);

    for b = 1:length(unique_names)
        
        %create a regular expression to find if there is a match within the
        %names string
        %the spaces added to the beginning and end ensure a substring will
        %not be chosen as a match (ex. 'earth' as a match for 'art')
        reg = strcat('\s',unique_names(b),'\s');
        duplicates(b,1) = regexp(names_joined,reg,'match'); 

        %locate where there are duplicates
        bool(b,1) = (length(duplicates{b,1}) > 1);
        num_rep(b,1) = length(duplicates{b,1});

    end

    %set initial indices 
    %index1 = start of duplicated entries
    %index2 = end of duplicated entries
    %if index1=index2 there are no duplicated entries
    index1 = 1;
    index2 = 0;
    values_max = zeros(length(unique_names),1);

    for c = 1:length(unique_names)

        %set index2 based on number of repeated entries and find max value
        %from entries
        index2 = index2 + num_rep(c,1);
        value_max = max(values_sorted(index1:index2,1)); 
        values_max(c,1) = value_max;

        %reset index1 based on index2
        index1 = index2 + 1;
    end
    
    %convert to table
    max_values = table(unique_names, values_max);
    
end 