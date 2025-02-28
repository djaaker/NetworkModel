function [indx,actual_distance] = find_closest(neuron_idx,target_distance,X,Y)

x_coord_neuron = X(neuron_idx);
y_coord_neuron = Y(neuron_idx);

D = sqrt((X-x_coord_neuron).^2 + (Y-y_coord_neuron).^2);

if max(D)<target_distance
    indx = -1;
    actual_distance = 0;
    return
end

% Subtracting target distance
D = D-target_distance;

% Making the distance to self infinite
D(neuron_idx) = Inf;

D = abs(D);

% find the minimum distance 
min_d =  min(D);

% finding all the elements with minimum distance as there can be multiple
% neurons with the same distance 
min_indxs = find(D==min_d);

% Get the random element using the generated index
indx = min_indxs(randi(length(min_indxs)));

actual_distance = D(indx);

end

