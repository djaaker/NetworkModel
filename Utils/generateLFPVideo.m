function generateLFPVideo(lfp_data, dt, filename)
    % GENERATE_LFP_VIDEO Creates a video of LFP activity using imagesc.
    %
    % Inputs:
    %   lfp_data - N x N x T matrix (LFP data over time)
    %   dt       - Time step between frames (seconds)
    %   filename - Output video file name (e.g., 'lfp_video.mp4')
    
    % Get matrix size
    [N1, N2, T] = size(lfp_data);
    
    % Create video writer
    writerObj = VideoWriter(filename); 
    writerObj.FrameRate = 1 / dt; % Frame rate based on simulation step
    open(writerObj);

    % Initialize figure
    figure;
    hImage = imagesc(lfp_data(:, :, 1)); % Initial plot
    colormap(jet); % Color map for visualization
    colorbar;
    caxis([min(lfp_data(:)), max(lfp_data(:))]); % Fixed color scale
    axis equal;
    axis off;
    titleHandle = title(sprintf('LFP Activity - Time: %.3f s', 0));

    % Loop through time and update CData
    for t = 1:T
        set(hImage, 'CData', lfp_data(:, :, t)); % Update data
        set(titleHandle, 'String', sprintf('LFP Activity - Time: %.3f s', t * dt)); % Update title
        drawnow;

        % Capture frame for video
        frame = getframe(gcf);
        writeVideo(writerObj, frame);
    end

    % Close video writer
    close(writerObj);
    disp(['Video saved as ', filename]);
end
