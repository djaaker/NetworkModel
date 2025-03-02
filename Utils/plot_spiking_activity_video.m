function plot_spiking_activity_video(SNN, saveflag , filename)
    % PLOT_SPIKING_ACTIVITY_VIDEO Generates a video of neuron spiking activity
   
    % Video settings
    if saveflag == 1
        frame_rate = 1 / 10*dt; % Frame rate based on simulation time step
        writerObj = VideoWriter(filename, 'MPEG-4'); 
        writerObj.FrameRate = frame_rate;
        open(writerObj);
    end

    % Number of neurons and time steps
    spike_matrix = SNN.spikes(1:SNN.N_E,:);
    num_time_steps = size(spike_matrix, 2);
    
    x_coords = cell2mat(arrayfun(@(s) s.x_coord, SNN.neurons,'UniformOutput',false));
    y_coords = cell2mat(arrayfun(@(s) s.y_coord, SNN.neurons,'UniformOutput',false));

    % Initialize figure and scatter plot
    figure;
    hold on;
    axis equal;
    xlim([min(x_coords)-1, max(x_coords)+1]);
    ylim([min(y_coords)-1, max(y_coords)+1]);
    title('Neuron Spiking Activity');

    % Initial scatter plot: hollow circles (unfilled)
    scatterPlot = scatter(x_coords, y_coords, 50, 'k'); % Hollow circles (default)
    
    % Loop through time steps and update colors dynamically
    for t = 1:num_time_steps
        spiking_neurons = find(spike_matrix(:, t)); % Neurons that spike
        
        % Set default marker appearance
        scatterPlot.MarkerEdgeColor = 'k'; % Black outline
        scatterPlot.MarkerFaceColor = 'none'; % Unfilled by default
        
        % Fill spiking neurons in red
        if ~isempty(spiking_neurons)
            scatterPlot.CData(spiking_neurons, :) = repmat([1 0 0], length(spiking_neurons), 1); % Red
            scatterPlot.MarkerFaceColor = 'flat'; % Apply color fill for spiking neurons
        end
        
        title(sprintf('Time: %.3f s', t * dt));
        drawnow;
        if saveflag == 1
            % Capture frame and write to video
            frame = getframe(gcf);
            writeVideo(writerObj, frame);
        end
    end
    if saveflag == 1
        close(writerObj);
        disp(['Video saved as ', filename]);
    end
end