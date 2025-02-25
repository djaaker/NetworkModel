%% Network Paramters 
%Number of neurons in one dimension 
parameters.grid_size = 50;
parameters.excitatory_ratio = 0.8;
% Size of the patch of cortex to simulate 
parameters.grid_length = 1e-3; % in m


% Decay constant for connectivity 
parameters.sigma = 0.4e-3; %% in m 
% Speed of AP 
parameters.vAP = 0.2; % in m/s
% Synaptic delay 
parameters.tau_syn = 300e-6; % in seconds

%% Simulation Parameters
% time step of simulation in s
parameters.dt = 0.1e-3;
% Total time of simulation in s 
parameters.total_time = 2;

%% Create spiling neural network SNN 
[SNN, spikes,v_neurons,parameters] = generateNetwork(parameters);


%% Saving session
savepath = uigetdir(path);
sessionName = [savepath,'/','SNNStimulation50x50PoissonInput.mat'];

%% Simulating network
h = waitbar(0, 'Initializing...'); % Initialize waitbar
total_iterations = parameters.total_time / parameters.dt;
start_time = tic; % Start timer
save_interval = total_iterations/20;

rate = 20;                       % Firing rate in Hz
duration = 1;    
spike_train_ext = gpuArray(generate_poisson_spikes_N(parameters.num_neurons, rate, duration, parameters.dt));

for i = 1:total_iterations

    I_ext = gpuArray(zeros(parameters.num_neurons, 1));
   
    if i <= size(spike_train_ext,2)
        spike_ext = squeeze(spike_train_ext(:,i));
    else
        spike_ext = gpuArray(false(parameters.num_neurons,1));
    end

    % Update network
    [SNN,spikes,v_neurons] = updateNetworkGPU(i,i*parameters.dt,I_ext,spike_ext,SNN,spikes,v_neurons,parameters.num_neurons,parameters.dt);
    
    % Calculate elapsed time and estimate time remaining
    elapsed_time = toc(start_time);
    avg_time_per_iter = elapsed_time / i;
    estimated_time_left = avg_time_per_iter * (total_iterations - i);
    
    % Update waitbar with progress and estimated time left
    waitbar(i / total_iterations, h, sprintf('SNN Simulation Progress: %.2f%% | Time left: %.2f sec',(i / total_iterations) * 100, estimated_time_left));
    if mod(i, save_interval) == 0 || i == total_iterations
        save(sessionName,"-v7.3");
        fprintf('Workspace saved at %d%% completion.\n', round((i / total_iterations) * 100));
    end
end

close(h); % Close waitbar when done

%% Plotting

figure;
plot(SNN.v_neurons(2,:));

%% 
figure();
imagesc(SNN.spikes);
colormap(flipud(gray))