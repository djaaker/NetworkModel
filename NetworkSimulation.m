%% Network Paramters 
<<<<<<< HEAD
target_neuron_density = 2.8125e10; % per m2
=======
%Number of neurons in one dimension 
grid_size = 150;
>>>>>>> d59761fc73fbfc3c0243bda96adde8fdfe8ef1f6
excitatory_ratio = 0.8;
grid_size = 50;
grid_length = 0.3e-3;
% Decay constant for connectivity 
sigma = 0.3e-3; %% in m 
% Speed of AP 
vAP = 0.2; % in m/s
% Synaptic delay 
tau_syn = 300e-6; % in seconds


%% Simulation Parameters
% time step of simulation in s
dt = 0.1e-3;
% Total time of simulation in s 
total_time = 0.2;

%% Create spiling neural network SNN 
SNN = Network(grid_size,grid_length,excitatory_ratio, sigma,vAP,tau_syn, dt,total_time);
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
total_iterations = total_time /  dt;
start_time = tic; % Start timer
save_interval = total_iterations/20;

%% Generating external input to the network
N = grid_size*grid_size;         % Number of neurons
rate_baseline = 20;              % baseline firing rate
rate_stimulation = 50;           % stimulation modulated firing rate
rate_baseline_post = 0; 
baseline_duration = 0.1; 
stim_duration = 0;
spike_train_ext = [generate_poisson_spikes_N(N, rate_baseline, baseline_duration, dt) generate_poisson_spikes_N(N, rate_stimulation, stim_duration, dt) generate_poisson_spikes_N(N, rate_baseline_post, total_time-baseline_duration-stim_duration, dt)];

%% Simulating network

for i = 1:total_iterations
    I_ext = zeros(SNN.num_neurons, 1);

    spike_ext = (squeeze(spike_train_ext(:,i))).*(~SNN.EI_tag)';

    % Update network
    SNN = SNN.update_par(i, i * dt, I_ext,spike_ext);
    
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