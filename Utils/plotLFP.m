function [fig] = plotLFP(LFP,SNN,spikes,spike_train_ext)
    
    LPF_stack = reshape(LFP.LFP_downsampled,LFP.N_blocks_combined^2,[]);
    
    time = 0:1/LFP.target_sampling_freq:(size(LPF_stack,2)-1)*1/LFP.target_sampling_freq;
    
    %% Ploting a stack plot of the LFP 
    n_plots = 5;
    fig = figure;
    ax1 = subplot(n_plots,1,1);
    stack_plot(LPF_stack,0,20,LFP.target_sampling_freq);
    ylabel('LFP Amplitude (a.u.)');xlabel(['Time ' ...
        '(s)']);
    title('Local Field Potential from all Electrodes');
    hold on;
    % Marking start up 
    y_limits = ylim; % Automatically get Y-axis limits
    x_startup = [0 0.2]; % Define the x-range to shade
    x_fill = [x_startup(1) x_startup(2) x_startup(2) x_startup(1)]; % Rectangle x-coordinates
    y_fill = [y_limits(1) y_limits(1) y_limits(2) y_limits(2)]; % Use current y-limits
    % Fill the shaded region
    fill(x_fill, y_fill, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Semi-transparent shading
    x_label = mean(x_startup); % Middle of shaded region
    y_label = y_limits(1) + 0.2*(y_limits(2)-y_limits(1)); % Middle of Y-range
    text(x_label, y_label, 'Start-up', 'Color', 'r', 'FontSize', 12, ...
        'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    
    % Marking stim 
    y_limits = ylim; % Automatically get Y-axis limits
    x_stim = [0.6 0.625]; % Define the x-range to shade
    x_fill = [x_stim(1) x_stim(2) x_stim(2) x_stim(1)]; % Rectangle x-coordinates
    y_fill = [y_limits(1) y_limits(1) y_limits(2) y_limits(2)]; % Use current y-limits
    % Fill the shaded region
    fill(x_fill, y_fill, 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Semi-transparent shading
    x_label = mean(x_stim); % Middle of shaded region
    y_label = y_limits(1) + 0.2*(y_limits(2)-y_limits(1)); % Middle of Y-range
    text(x_label, y_label, 'Stim', 'Color', 'b', 'FontSize', 12, ...
        'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    
    ax2 = subplot(n_plots,1,2);
    plot(time,mean(LPF_stack,1),'LineWidth',2,'Color',[0.1 0.2 0.9]);
    ylabel('LFP Amplitude (a.u.)');xlabel(['Time ' ...
        '(s)']);
    hold on;title('Mean Local Field Potential from all Electrodes');
    % Marking start up 
    y_limits = ylim; % Automatically get Y-axis limits
    x_startup = [0 0.2]; % Define the x-range to shade
    x_fill = [x_startup(1) x_startup(2) x_startup(2) x_startup(1)]; % Rectangle x-coordinates
    y_fill = [y_limits(1) y_limits(1) y_limits(2) y_limits(2)]; % Use current y-limits
    % Fill the shaded region
    fill(x_fill, y_fill, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Semi-transparent shading
    x_label = mean(x_startup); % Middle of shaded region
    y_label = y_limits(1) + 0.2*(y_limits(2)-y_limits(1)); % Middle of Y-range
    text(x_label, y_label, 'Start-up', 'Color', 'r', 'FontSize', 12, ...
        'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    
    % Marking stim 
    y_limits = ylim; % Automatically get Y-axis limits
    x_stim = [0.6 0.625]; % Define the x-range to shade
    x_fill = [x_stim(1) x_stim(2) x_stim(2) x_stim(1)]; % Rectangle x-coordinates
    y_fill = [y_limits(1) y_limits(1) y_limits(2) y_limits(2)]; % Use current y-limits
    % Fill the shaded region
    fill(x_fill, y_fill, 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Semi-transparent shading
    x_label = mean(x_stim); % Middle of shaded region
    y_label = y_limits(1) + 0.2*(y_limits(2)-y_limits(1)); % Middle of Y-range
    text(x_label, y_label, 'Stim', 'Color', 'b', 'FontSize', 12, ...
        'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    ylim([y_limits]);
    
    ax3 = subplot(n_plots,1,3);
    time_spikes = 0:SNN.dt:(size(spike_train_ext,2)-1)*SNN.dt;
    imagesc(time_spikes,1:size(spike_train_ext,1),spike_train_ext);
    colormap(flipud(gray));xlabel('Time (s)');ylabel('Neurons');
    title('External spiking input to the excitatory network');
    yyaxis right;
    FR_binsize = 10e-3;
    firingRates_input = estimateFiringRate(spike_train_ext,FR_binsize,1/SNN.dt);
    time_firing_rate = FR_binsize/2:FR_binsize:FR_binsize/2+(size(firingRates_input,2)-1)*FR_binsize;
    plot(time_firing_rate,mean(firingRates_input,1),'LineWidth',2,'Color',[1 0 0]);
    ylabel('Mean firing rate (Hz)');
    
    ax4 = subplot(n_plots,1,4);
    time_spikes = 0:SNN.dt:(size(spikes,2)-1)*SNN.dt;
    imagesc(time_spikes,1:SNN.N_E,spikes(1:SNN.N_E,:));
    colormap(flipud(gray));xlabel('Time (s)');ylabel('Neurons');
    title('Spiking rater of the excitatory network');
    yyaxis right;
    % FR_binsize = 10e-3;
    firingRates_input = estimateFiringRate(spikes(1:SNN.N_E,:),FR_binsize,1/SNN.dt);
    time_firing_rate = FR_binsize/2:FR_binsize:FR_binsize/2+(size(firingRates_input,2)-1)*FR_binsize;
    plot(time_firing_rate,mean(firingRates_input,1),'LineWidth',2,'Color',[1 0 0]);
    ylabel('Mean firing rate (Hz)'); ylim([0 50]);
    
    ax5 = subplot(n_plots,1,5);
    time_spikes = 0:SNN.dt:(size(spikes,2)-1)*SNN.dt;
    imagesc(time_spikes,1:SNN.N_I,spikes(SNN.N_E+1:end,:));
    colormap(flipud(gray));xlabel('Time (s)');ylabel('Neurons');
    title('Spiking rater of the inhibitory network');
    yyaxis right;
    % FR_binsize = 10e-3;
    firingRates_input = estimateFiringRate(spikes(SNN.N_E+1:end,:),FR_binsize,1/SNN.dt);
    time_firing_rate = FR_binsize/2:FR_binsize:FR_binsize/2+(size(firingRates_input,2)-1)*FR_binsize;
    plot(time_firing_rate,mean(firingRates_input,1),'LineWidth',2,'Color',[1 0 0]);
    ylabel('Mean firing rate (Hz)');ylim([0 50]);
    
    linkaxes([ax1,ax2,ax3,ax4,ax5],'x');
    xlim([0 1.05]);
end