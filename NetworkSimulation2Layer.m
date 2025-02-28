%% Network Paramters 
%Number of neurons in one dimension 
grid_size = 150;
excitatory_ratio = 0.8;
% Size of the patch of cortex to simulate 
grid_length = 0.9e-3; % in m


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
total_time = 0.8;

%% Create spiling neural network SNN 
SNN = Network2Layer(grid_size,grid_length,excitatory_ratio, sigma,vAP,tau_syn,time_step,total_time);
% Checking structure 
% SNN.plot_network();
% SNN.plot_neuron_connection(45);

%% Saving session
saveFlag = 1;
if saveFlag == 1
    savepath = uigetdir(path);
    sessionName = [savepath,'/','SNN2LayerStimulation150x150PoissonInputModulated.mat'];
end

%% Simulating network
h = waitbar(0, 'Initializing...'); % Initialize waitbar
total_iterations = total_time / time_step;
start_time = tic; % Start timer
save_interval = total_iterations/10;

%% Generating external input to the network
N = SNN.num_neurons;             % Number of neurons
rate_baseline = 20;              % baseline firing rate
rate_stimulation = 40;           % stimulation modulated firing rate
rate_baseline_post = 0; 
baseline_duration = 0.2; 
stim_duration = 0.2;
dt = time_step; 
spike_train_ext = [generate_poisson_spikes_N(N, rate_baseline, baseline_duration, dt) generate_poisson_spikes_N(N, rate_stimulation, stim_duration, dt) generate_poisson_spikes_N(N, rate_baseline_post, total_time-baseline_duration-stim_duration, dt)];

%% Simulating network

for i = 1:total_iterations
    I_ext = zeros(SNN.num_neurons, 1);

    spike_ext = (squeeze(spike_train_ext(:,i))).*(~SNN.EI_tag)';

    % Update network
    SNN = SNN.update(i, i * time_step, I_ext,spike_ext);
    
    % Calculate elapsed time and estimate time remaining
    elapsed_time = toc(start_time);avg_time_per_iter = elapsed_time / i;estimated_time_left = (avg_time_per_iter * (total_iterations - i))/3600;
    
    % Update waitbar with progress and estimated time left
    waitbar(i / total_iterations, h,sprintf('SNN Simulation Progress: %.2f%% | Time left: %.2f hrs',(i / total_iterations) * 100, estimated_time_left));
    if (mod(i, save_interval) == 0 || i == total_iterations) && saveFlag == 1
        save(sessionName,"-v7.3");fprintf('Workspace saved at %d%% completion.\n', round((i / total_iterations) * 100));
    end
end

close(h); % Close waitbar when done

%% Plotting

plotSpikingNetwork(SNN)



% figure;
% plot(SNN.v_neurons(4,:));hold on;
% plot(SNN.spikes(3,:));
% ylim([-1 1.5])
% 
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

% %% 
% figure();
% imagesc(SNN.spikes);
% colormap(flipud(gray))
% 
% 
figure();
imagesc(spike_train_ext);
colormap(flipud(gray))

abc = 9;
figure();
plot(SNN.ge_neurons(abc,:));
hold on;
plot(SNN.gi_neurons(abc,:));
yyaxis right;
plot(spike_train_ext(abc,:));

figure();
plot(SNN.v_neurons(abc,:));
hold on;
yyaxis right;
plot(SNN.spikes(abc,:));