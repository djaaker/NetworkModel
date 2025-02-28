function [indx] = find_closest(neuron_idx,target_distance,X,Y)

x_coord_neuron = X(neuron_idx);
y_coord_neuron = Y(neuron_idx);

D = sqrt((X-x_coord_neuron).^2 + (Y-y_coord_neuron).^2);



end

