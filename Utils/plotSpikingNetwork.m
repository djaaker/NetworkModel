function [fig] = plotSpikingNetwork(SNN,spikes,ge_neurons,gi_neurons)
    % Extract spike data
    spike_matrix = spikes;  % Binary matrix of spikes (neurons x time)
    EI_tag = SNN.EI_tag;        % 0 for excitatory, 1 for inhibitory
    dt = SNN.dt;                % Time step of simulation
    time_vector = (1:size(spike_matrix, 2)) * dt; % Time in seconds

    % Initialize raster image with a white background
    [num_neurons, num_time_steps] = size(spike_matrix);
    raster_image = ones(num_neurons, num_time_steps, 3); % RGB image (white background)

    % Assign colors for spikes
    excitatory_color = [0, 0, 1]; % Blue for excitatory spikes
    inhibitory_color = [1, 0, 0]; % Red for inhibitory spikes

    % Overlay spikes with appropriate colors
    for ch = 1:3  % Apply to each color channel
        raster_image(:,:,ch) = raster_image(:,:,ch) - (spike_matrix .* (EI_tag' == 0) * (1 - excitatory_color(ch))) ...
                                                   - (spike_matrix .* (EI_tag' == 1) * (1 - inhibitory_color(ch)));
    end

    % Compute firing rates
    window_size = round(0.05 / dt); % 50 ms window size
    firing_rates = movsum(spike_matrix, window_size, 2) / (window_size * dt); % Firing rate in Hz
    mean_firing_rate = mean(firing_rates, 1); % Average across neurons


    mean_ge = mean(ge_neurons, 1);
    mean_gi = mean(gi_neurons, 1);

    % Create figure with subplots
    fig = figure;
    
    % Raster plot (top)
    subplot(2,1,1);
    imagesc(time_vector, 1:num_neurons, raster_image);
    xlabel('');
    ylabel('Neuron Index');
    title('Spike Rastergram');
    set(gca, 'YDir', 'normal'); % Ensure neuron index is displayed correctly

    % Firing rate & conductance plot (bottom)
    subplot(2,1,2);
    yyaxis left;
    plot(time_vector, mean_firing_rate, 'k', 'LineWidth', 2);
    ylabel('Avg Firing Rate (Hz)');
    title('Population Firing Rate & Conductance');

    yyaxis right;
    plot(time_vector, mean_ge, 'b', 'LineWidth', 2); hold on;
    plot(time_vector, mean_gi, 'r', 'LineWidth', 2);
    ylabel('Conductance (S)');
    legend({'Firing Rate', 'ge (exc)', 'gi (inh)'}, 'Location', 'northeast');
    hold off;
    
    xlabel('Time (s)');

    % Link x-axes
    linkaxes(findobj(gcf, 'Type', 'axes'), 'x');
    
end
