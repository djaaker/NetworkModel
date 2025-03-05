function [SNN] =  NetworkSimulation2Layer_CPP(grid_size,excitatory_ratio,grid_length,sigma,vAP,tau_syn,kden,GE,GI,time_step,total_time,spike_train_ext)
    %#codegen
    % n_time = round(total_time/time_step);
    %% Create spiking neural network SNN 
    SNN = Network2Layer_CPP(grid_size,grid_length,excitatory_ratio, sigma,vAP,tau_syn,time_step,total_time,GE,GI,kden);
    
    %% Simulating network
    % h = waitbar(0, 'Initializing...'); % Initialize waitbar
    total_iterations = total_time / time_step;
    start_time = tic; % Start timer
    %% Simulating network
    
    for i = 1:total_iterations
        I_ext = zeros(SNN.num_neurons, 1);
    
        spike_ext = (squeeze(spike_train_ext(:,i))).*(~SNN.EI_tag)';
    
        % Update network
        SNN = SNN.update(i, i * time_step, I_ext,spike_ext);
        
        % Calculate elapsed time and estimate time remaining
        elapsed_time = toc(start_time);avg_time_per_iter = elapsed_time / i;estimated_time_left = (avg_time_per_iter * (total_iterations - i))/3600;
        disp(estimated_time_left);
        % Update waitbar with progress and estimated time left
        % waitbar(i / total_iterations, h,sprintf('SNN Simulation Progress: %.2f%% | Time left: %.2f hrs',(i / total_iterations) * 100, estimated_time_left));
        sprintf('SNN Simulation Progress: %.2f%% | Time left: %.2f hrs',(i / total_iterations) * 100, estimated_time_left);
    end
    
    % close(h); % Close waitbar when done

end