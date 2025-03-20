function plot_spiking_activity_video(SNN, spikes,dt, saveflag , filename)
    % PLOT_SPIKING_ACTIVITY_VIDEO Generates a video of neuron spiking activity
   
    % Video settings
    if saveflag == 1
        frame_rate = 1 / (100*dt); % Frame rate based on simulation time step
        writerObj = VideoWriter(filename); 
        writerObj.FrameRate = frame_rate;
        open(writerObj);
    end

    % Number of neurons and time steps
    spike_matrix = spikes(1:SNN.N_E,:);
    num_time_steps = size(spike_matrix, 2);

    % Initialize figure and scatter plot
    figure;
    colormap([1 1 1; 1 0 0]); % White (no spike), Red (spike)
    colorbar;
    caxis([0 1]); % Ensure consistent scaling
    axis equal;
    axis off;
    box off;

    % Create neuron activity grid
    activity_grid = zeros(SNN.grid_size_e, SNN.grid_size_e);

    % Initialize imagesc plot handle
    hImage = imagesc(activity_grid);
    titleHandle = title('Neuron Spiking Activity - Time: 0.000 s');

    % Loop through time steps and update colors dynamically
    for t = 1:num_time_steps
        % Reset activity grid
        activity_grid = reshape(spike_matrix(:,t),SNN.grid_size_e,SNN.grid_size_e,[]); 

        % Update imagesc data instead of re-plotting
        set(hImage, 'CData', activity_grid);
        set(titleHandle, 'String', sprintf('Neuron Spiking Activity - Time: %.3f s', t * dt));
%         drawnow;
        pause(0.00001);
        if saveflag == 1
            % Capture frame for video
            frame = getframe(gcf);
            writeVideo(writerObj, frame);
        end
    end
    if saveflag == 1
        close(writerObj);
        disp(['Video saved as ', filename]);
    end
end