function[unique_sub, subsystem_fluxes] = subsystem_mapping(subsystems,flux_matrix)

%{ 
This function maps a knockOut flux matrix to the appropriate subsystem and
averages the rows within each subsystem
Inputs: 
subsystems: cell array containing the subsystems from the desired model
flux_matrix: the flux matrix obtained from performing gene knockouts 

%}

    % Get a list of unique subsystems
    unique_sub = unique(subsystems);

    % cycle through the unique subsystems and all subsystems to match 
    % each subsystem to the unique subsystem list
    locations = zeros(length(subsystems),1);
    for a = 1:length(unique_sub)
        for b = 1:length(subsystems)
            if strcmp(unique_sub(a), subsystems(b))
                locations(b) = a;
            end
        end
    end

    % cycle through all of the unique subsystems and perform a row wise 
    % average to get the subsystem average flux for each gene
    subsystem_fluxes = zeros(size(flux_matrix,1), length(unique_sub));
    for a = 1:length(unique_sub)

        temp = flux_matrix(:, (locations == a) );
        % take the row wise average
        avg = mean(temp,2);
        subsystem_fluxes(:,a) = avg;
        
    end

end

