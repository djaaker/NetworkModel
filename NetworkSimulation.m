%% Network Paramters 
%Number of neurons in one dimension 
grid_size = 10;
excitatory_ratio = 0.8;
% Size of the patch of cortex to simulate 
grid_length = 0.25e-3; % in m


% Decay constant for connectivity 
sigma = 0.4e-3; %% in m 
% Speed of AP 
vAP = 0.2; % in m/s
% Synaptic delay 
tau_syn = 300e-6; % in seconds

%% Simulation Parameters
% time step of simulation in s
time_step = 0.1e-3;
% Total time of simulation in s 
total_time = 0.4;

%% Create spiling neural network SNN 
SNN = Network(grid_size,grid_length,excitatory_ratio, sigma,vAP,tau_syn,time_step,total_time);
% Checking structure 
% SNN.plot_network();
% SNN.plot_neuron_connection(45);

%% Saving session
saveFlag = 0;
if saveFlag == 1
    savepath = uigetdir(path);
    sessionName = [savepath,'/','SNNStimulation20x20PoissonInputModulated.mat'];
end

%% Simulating network
h = waitbar(0, 'Initializing...'); % Initialize waitbar
total_iterations = total_time / time_step;
start_time = tic; % Start timer
save_interval = total_iterations/20;

N = grid_size*grid_size;         % Number of neurons
rate_baseline = 20;              % baseline firing rate
rate_stimulation = 50;           % stimulation modulated firing rate
baseline_duration = 0.2;    
dt = time_step; 
spike_train_ext = [generate_poisson_spikes_N(N, rate_baseline, baseline_duration, dt) generate_poisson_spikes_N(N, rate_stimulation, total_time-baseline_duration, dt)];

for i = 1:total_iterations
    % % Apply external current condition
    % if i * time_step > 0.001 && i * time_step < 1
    %     I_ext = 0.5e-9 * ones(SNN.num_neurons, 1);
    % else
    I_ext = zeros(SNN.num_neurons, 1);
    % end
   
%     if i <= size(spike_train_ext,2)
%         spike_ext = (squeeze(spike_train_ext(:,i))).*(~SNN.EI_tag)';
%     else
%         spike_ext = false(N,1);
%     end

    spike_ext = (squeeze(spike_train_ext(:,i))).*(~SNN.EI_tag)';

    % Update network
    SNN = SNN.update_par(i, i * time_step, I_ext,spike_ext);
    
    % Calculate elapsed time and estimate time remaining
    elapsed_time = toc(start_time);
    avg_time_per_iter = elapsed_time / i;
    estimated_time_left = (avg_time_per_iter * (total_iterations - i))/3600;
    
    % Update waitbar with progress and estimated time left
    waitbar(i / total_iterations, h, ...
        sprintf('SNN Simulation Progress: %.2f%% | Time left: %.2f hrs', ...
        (i / total_iterations) * 100, estimated_time_left));
    if (mod(i, save_interval) == 0 || i == total_iterations) && saveFlag == 1
        save(sessionName,"-v7.3");
        fprintf('Workspace saved at %d%% completion.\n', round((i / total_iterations) * 100));
    end
end

close(h); % Close waitbar when done

%% Plotting

figure;
plot(SNN.v_neurons(3,:));hold on;
plot(spike_train_ext(3,:));
ylim([-1 1.5])

figure();
subplot(4,1,1)
imagesc(SNN.v_neurons);
subplot(4,1,2)
imagesc(SNN.ge_neurons);
subplot(4,1,3)
imagesc(SNN.gi_neurons);
subplot(4,1,4)
imagesc(SNN.spikes);
colormap(flipud(gray))

%% 
figure();
imagesc(SNN.spikes);
colormap(flipud(gray))


figure();
imagesc(spike_train_ext);
colormap(flipud(gray))
