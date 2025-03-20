function [fig] = plotNetworkProperties(select_neuron,SNN,spikes)
    
%% Plotting connections

x_coords = cell2mat(cellfun(@(s) s.x_coord,SNN.neurons,'UniformOutput',false));
x_coords_i = x_coords(SNN.N_E+1:end);
x_coords = x_coords(1:SNN.N_E); % Selecting just the excitatory neurons

y_coords = cell2mat(cellfun(@(s) s.y_coord,SNN.neurons,'UniformOutput',false));
y_coords_i = y_coords(SNN.N_E+1:end);
y_coords = y_coords(1:SNN.N_E); % Selecting just the excitatory neurons

connections = SNN.neurons{select_neuron}.connections_E;

connections_i = SNN.neurons{select_neuron}.connections_I - SNN.N_E;

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
fig = figure;
subplot(2,2,[1 3]);
hold on;
% Plot the Gaussian function using imagesc
imagesc(x_coords, y_coords, Z); 
colormap(bone);
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
xlim([-SNN.grid_length/2 SNN.grid_length/2]);
ylim([-SNN.grid_length/2 SNN.grid_length/2]);


subplot(2,2,2);
total_time = size(spikes,2)*SNN.dt;
post_stim = 0.3; % in seconds
mean_FR_e = sum(spikes(1:SNN.N_E,0.3/SNN.dt:end),2)/(total_time-post_stim); % spikes/sec
mean_FR_i = sum(spikes(SNN.N_E+1:end,0.3/SNN.dt:end),2)/(total_time-post_stim); % spikes/sec

[counts, edges] = histcounts(mean_FR_e,20);
binCenters = edges(1:end-1) + diff(edges)/2;
plot(binCenters, counts/SNN.N_E,'Color',[0 0 0 0.4], 'LineWidth', 2); hold on;
xlabel('Mean firing rate (spks/sec)');ylabel('Normalized Density');

[counts, edges] = histcounts(mean_FR_i, 10);
binCenters = edges(1:end-1) + diff(edges)/2;
plot(binCenters, counts/SNN.N_I,'Color',[255/256 132/256 0 1], 'LineWidth', 2); hold on;

subplot(2,2,4);
% Assume spikeMatrix is an N x T matrix (N = neurons, T = time steps)
[N, ~] = size(spikes);  % Get matrix size
% Find spike times for all neurons
[neuronIdx, spikeTimes] = find(spikes');  % Transpose to get time first
% Convert to cell array: each neuron has its own spike time list
spikeTimesCell = accumarray(neuronIdx, spikeTimes, [N, 1], @(x){sort(x)});
% Compute interspike intervals (ISI)
ISI = cellfun(@diff, spikeTimesCell, 'UniformOutput', false);
% Filter neurons with at least 3 spikes
validNeurons = cellfun(@(x) numel(x) >= 2, ISI);
% Compute CV only for valid neurons
cv_values = nan(N, 1);  % Preallocate
cv_values(validNeurons) = cellfun(@(x) std(x) / mean(x), ISI(validNeurons));
% Remove NaN values (neurons with <3 spikes)
valid_cv = cv_values(~isnan(cv_values));

[counts, edges] = histcounts(valid_cv, 50);
binCenters = edges(1:end-1) + diff(edges)/2;
plot(binCenters, counts/numel(valid_cv),'Color',[0 0 0 0.4], 'LineWidth', 2); hold on;
xlabel('Mean firing rate (spks/sec)');ylabel('Normalized Density');xlim([0 2]);


end