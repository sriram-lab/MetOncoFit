function[unique_sub, subsystem_fluxes] = subsystem_mapping(subsystems,flux_matrix)

%{ 
This function maps a knockOut flux matrix to the appropriate subsystem and
averages the rows within each subsystem
Inputs: 
subsystems: cell array containing the subsystems from the desired model
flux_matrix: the flux matrix obtained from performing gene knockouts 

%}

    unique_sub = unique(subsystems);

    locations = zeros(length(subsystems),1);
    for a = 1:length(unique_sub)
        for b = 1:length(subsystems)
            if strcmp(unique_sub(a), subsystems(b))
                locations(b) = a;
            end
        end
    end

    subsystem_fluxes = zeros(size(flux_matrix,1), length(unique_sub));
    for a = 1:length(unique_sub)

        temp = flux_matrix(:, (locations == a) );
        avg = mean(temp,2);
        subsystem_fluxes(:,a) = avg;
    end

end

