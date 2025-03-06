function [fig] = plotNetworkProperties(select_neuron,SNN)
    
%% Plotting connections

x_coords = cell2mat(arrayfun(@(s) s.x_coord,SNN.neurons,'UniformOutput',false));
x_coords_i = x_coords(SNN.N_E+1:end);
x_coords = x_coords(1:SNN.N_E); % Selecting just the excitatory neurons

y_coords = cell2mat(arrayfun(@(s) s.y_coord,SNN.neurons,'UniformOutput',false));
y_coords_i = y_coords(SNN.N_E+1:end);
y_coords = y_coords(1:SNN.N_E); % Selecting just the excitatory neurons

connections = SNN.neurons(select_neuron).connections_E;

connections_i = SNN.neurons(select_neuron).connections_I - SNN.N_E;

X = reshape(x_coords,SNN.grid_size_e,[]);
Y = reshape(y_coords,SNN.grid_size_e,[]);
% Define Gaussian function parameters
mu_x = x_coords(select_neuron);     % Mean in X direction
mu_y = y_coords(select_neuron);     % Mean in Y direction
sigma_x = 400e-6;  % Standard deviation in X direction
sigma_y = 400e-6;  % Standard deviation in Y direction

% 2D Gaussian function
Z = exp(-(((X - mu_x).^2 / (2 * sigma_x^2)) + ((Y - mu_y).^2 / (2 * sigma_y^2))));


% Plot all neurons
fig = figure; hold on;
% Plot the Gaussian function using imagesc
imagesc(x_coords, y_coords, Z); 
colormap(sky);
colorbar;
axis xy; % Ensures the Y-axis is not flipped
xlabel('X Coordinate');ylabel('Y Coordinate');title('Neural Network Visualization');

highlight_color = 'r';

% Draw rectangle to represent network borders
rectangle('Position', [min(x_coords), min(y_coords), max(x_coords)-min(x_coords), max(y_coords)-min(y_coords)], 'EdgeColor', 'k', 'LineWidth', 2);


% Plot connections of the selected neuron
for j = 1:length(connections_i)
    conn_neuron = connections_i(j);
    plot([x_coords(select_neuron), x_coords_i(conn_neuron)], ...
         [y_coords(select_neuron), y_coords_i(conn_neuron)], ...
         'Color', [255/256 132/256 0 1],'LineWidth', 1.5);
end

% Plot connections of the selected neuron
for j = 1:length(connections)
    conn_neuron = connections(j);
    plot([x_coords(select_neuron), x_coords(conn_neuron)], ...
         [y_coords(select_neuron), y_coords(conn_neuron)], ...
         'Color', [0 0 0 0.4],'LineWidth', 1.5);
end

% Highlight target neuron
scatter(x_coords(select_neuron), y_coords(select_neuron), 50, highlight_color, 'filled');
box off;


end