%% Network Paramters 
%Number of neurons in one dimension 
grid_size = 100;
excitatory_ratio = 0.8;
% Size of the patch of cortex to simulate 
grid_length = 0.6e-3; % in m


%% Connectivity Parameters
sigma = 0.4e-3; %% in m 
% Speed of AP 
vAP = 0.2; % in m/s
% Synaptic delay 
tau_syn = 400e-6; % in seconds
kden = 0.0075;

%% Synapse parameters
GE = 3e-9;
GI = 30e-9;


%% Simulation Parameters
% time step of simulation in s
time_step = 0.1e-3;
% Total time of simulation in s 
total_time = 0.4;

%% Saving session
saveFlag = 0;
if saveFlag == 1
    savepath = uigetdir(path);
    sessionName = [savepath,'/','SNN2LayerStimulation300x200startup_20hz_200ms_stim_cpp.mat'];
    % filename = [savepath,'/','SNN2LayerStimulation300x300startup_20hz_200ms.avi'];
end

%% Calculating some network parameters for generating external spikes
target_neuron_density = 2.8125e10; % per m2
% Number of E and I neurons
num_neurons = grid_size^2;
N_E = (ceil(sqrt(excitatory_ratio* num_neurons)))^2;
N_I = (ceil(sqrt((1-excitatory_ratio)* num_neurons)))^2;
grid_size_e = sqrt(N_E);
grid_size_i = sqrt(N_I);
num_neurons =  N_E+N_I;

%% Generating external input to the network
N = num_neurons;             % Number of neurons
rate_startup = 20;              % baseline firing rate
rate_stimulation = 200;           % stimulation modulated firing rate
rate_baseline = 0; 

startup_duration = 0.2; 
baseline_prestim = 0.2;
stim_duration = 0;
dt = time_step; 
spike_train_ext = [generate_poisson_spikes_N(N, rate_startup, startup_duration, dt) generate_poisson_spikes_N(N, rate_baseline, baseline_prestim, dt) generate_poisson_spikes_N(N, rate_stimulation, stim_duration, dt) generate_poisson_spikes_N(N, rate_baseline, total_time-startup_duration-baseline_prestim-stim_duration, dt)];

%% Simulation
[SNN] =  NetworkSimulation2Layer_CPP(grid_size,excitatory_ratio,grid_length,sigma,vAP,tau_syn,kden,GE,GI,time_step,total_time,spike_train_ext);

%% Saving workspace
save(sessionName,"-v7.3");