
%% Network Paramters 
%Number of neurons in one dimension 
grid_size = 250;
excitatory_ratio = 0.8;
% Size of the patch of cortex to simulate 
grid_length = 1.5e-3; % in m


% Decay constant for connectivity 
sigma = 0.4e-3; %% in m 
% Speed of AP 
vAP = 0.2; % in m/s
% Synaptic delay 
tau_syn = 400e-6; % in seconds

%% Simulation Parameters
% time step of simulation in s
time_step = 0.1e-3;
% Total time of simulation in s 
total_time = 0.6;

%% Saving session
saveFlag = 1;
if saveFlag == 1
    savepath = uigetdir(path);
    sessionName = [savepath,'/','SNN2LayerStimulation250x250startup_20hz_200ms.mat'];
    filename = [savepath,'/','SNN2LayerStimulation250x250startup_20hz_200ms.avi'];
    lfpfilename = [savepath,'/','SNN2LayerStimulation250x250startup_20hz_200ms_LFP.avi'];
    lfpfilename2 = [savepath,'/','SNN2LayerStimulation250x250startup_20hz_200ms_LFPgauss.avi'];
end

%% Create spiking neural network SNN 
SNN = Network2Layer(grid_size,grid_length,excitatory_ratio, sigma,vAP,tau_syn,time_step,total_time);
% Checking structure 
% SNN.plot_network();
% SNN.plot_neuron_connection(45);

%% Simulating network
h = waitbar(0, 'Initializing...'); % Initialize waitbar
total_iterations = total_time / time_step;
start_time = tic; % Start timer
save_interval = total_iterations/10;

%% Generating external input to the network
N = SNN.num_neurons;             % Number of neurons
rate_startup = 20;              % baseline firing rate
rate_stimulation = 200;           % stimulation modulated firing rate
rate_baseline = 0; 

startup_duration = 0.2; 
baseline_prestim = 0.4;
stim_duration = 0.0;
dt = time_step; 
spike_train_ext = [generate_poisson_spikes_N(N, rate_startup, startup_duration, dt) generate_poisson_spikes_N(N, rate_baseline, baseline_prestim, dt) generate_poisson_spikes_N(N, rate_stimulation, stim_duration, dt) generate_poisson_spikes_N(N, rate_baseline, total_time-startup_duration-baseline_prestim-stim_duration, dt)];

frac_neurons_ext_input = 1;
rand_indx_expt_input = (rand(N,1)<frac_neurons_ext_input).*(~SNN.EI_tag)';


%% Simulating network

for i = 1:total_iterations
    I_ext = zeros(SNN.num_neurons, 1);

    spike_ext = (squeeze(spike_train_ext(:,i))).*rand_indx_expt_input;

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