
%% Network Paramters 
%Number of neurons in one dimension 
grid_size = 10;
excitatory_ratio = 0.8;
% Size of the patch of cortex to simulate 
grid_length = 2e-3; % in mm


% Decay constant for connectivity 
sigma = 0.4e-3; %% in mm 
% Speed of AP 
vAP = 0.2*1000*1e-3; % in mm/s
% Synaptic delay 
tau_syn = 300e-6; % in seconds

%% Simulation Parameters
% time step of simulation in s
time_step = 0.2e-3;
% Total time of simulation in s 
total_time = 0.5;

%% Create spiling neural network SNN 
SNN = Network(grid_size,grid_length,excitatory_ratio, sigma,vAP,tau_syn,time_step,total_time);
% Checking structure 
% SNN.plot_network();
% SNN.plot_neuron_connection(45);

%% Simulating network

for i=1:total_time/time_step
    if i*time_step > 0.01 && i*time_step<0.2 
        I_ext = 10e-9*(ones(SNN.num_neurons,1)); 
    else
        I_ext = zeros(SNN.num_neurons,1);
    end
    SNN = SNN.update(time_step,i*time_step,I_ext);
    i
end