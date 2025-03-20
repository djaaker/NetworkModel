format compact;
% set(0,'DefaultFigureWindowStyle','normal')
addpath(genpath('Utils'));
addpath(genpath('Wave'));
addpath(genpath('Dependancies'));
addpath(genpath('CPP_Implementation'));

%% Get directory to save the plots and videos
% Prompt user to select a directory
save_dir = uigetdir('', 'Select Directory to save plots and videos');
folderpath1 = [save_dir, '/' , 'Plots'];
mkdir(folderpath1);

%% Loading saved files
spikes_binfile = [save_dir,'/','spikes'];
v_neurons_binfile = [save_dir,'/','v_neurons'];
gi_neurons_binfile = [save_dir,'/','gi_neurons'];
ge_neurons_binfile = [save_dir,'/','ge_neurons'];
snn_matfile = [save_dir,'/','SNN.mat'];

spikes = read_logical_matrix_bin(spikes_binfile);
v_neurons = read_matrix_bin(v_neurons_binfile);
ge_neurons = read_matrix_bin(ge_neurons_binfile);
gi_neurons = read_matrix_bin(gi_neurons_binfile);

v_neurons = single(v_neurons);
ge_neurons = single(ge_neurons);
gi_neurons = single(gi_neurons);
load (snn_matfile);

%% 
spikes = spikes(:,1:round(total_time/dt));
v_neurons = v_neurons(:,1:round(total_time/dt));
ge_neurons = ge_neurons(:,1:round(total_time/dt));
gi_neurons = gi_neurons(:,1:round(total_time/dt));

%% Plot network properties
fig(1) = plotNetworkProperties(35000,SNN,spikes);
savefig(fig(1),[folderpath1,'/','networkConnections.fig']);

%% Plotting and generating video for spiking activity
filename = [folderpath1, '/', 'spiking.avi'];
plot_spiking_activity_video(SNN,spikes,dt,1,filename);

%% PLotting single neuron dynamics
fig(2) = plotSingleNeuronDynamics(35000,SNN,spikes,v_neurons,ge_neurons,gi_neurons);
savefig(fig(2),[folderpath1,'/','singleNeuronDynamic.fig']);

%% Generating LFP
LFP = calculateLFP(SNN,v_neurons,ge_neurons,gi_neurons);

% Plotting LFP, spikes and input 
fig(3) = plotLFP(LFP,SNN,spikes,spike_train_ext);
savefig(fig(3),[folderpath1,'/','LFP_Spikes.fig']);

% Genreating a video of the LFP
lfpfilename = [folderpath1,'/','LFP_block.avi'];
generateLFPVideo(LFP.LFP_downsampled, 1/LFP.target_sampling_freq, lfpfilename);
lfpfilename = [folderpath1,'/','LFP_gauss.avi'];
generateLFPVideo(LFP.LFP_gaussian_downsampled, 1/LFP.target_sampling_freq, lfpfilename);

%% Wave detection
% Filtering and GP 
filterOrder = 4;
filterLP1 = 5;
filterLP2 = 40;
LFP.xf = bandpass_filter(LFP.LFP_gaussian_downsampled,filterLP1,filterLP2,filterOrder,LFP.target_sampling_freq);
[LFP.xgp, LFP.wt] = generalized_phase(LFP.xf,LFP.target_sampling_freq,0);

% Initializing parameters for wave detection
parameters.rows = LFP.N_blocks_combined;parameters.cols = LFP.N_blocks_combined;
[parameters.X,parameters.Y] = meshgrid( 1:parameters.cols, 1:parameters.rows );
SNN.grid_length = grid_length;
parameters.xspacing = LFP.lfp_patch_size*(SNN.grid_length/SNN.grid_size_e)*1e3 ;parameters.yspacing = LFP.lfp_patch_size*(SNN.grid_length/SNN.grid_size_e)*1e3; % in mm need to change 
parameters.rhoThres= 0.5;

% Wave detection
allwaves.LFPIndex = (1:1:size(LFP.LFP_combined,3))';
xf{1,1} = LFP.xf;
xgp{1,1} = LFP.xgp;
wt{1,1} = LFP.wt;
Wavesall = detectWaves(xf,xgp,wt,allwaves,parameters,parameters.rhoThres);

%% Plotting
fig(4) = plotSpikingNetwork(SNN,spikes,ge_neurons,gi_neurons);
savefig(fig(4),[folderpath1,'/','spikingActivity.fig']);

%%

% figure;
% plot(SNN.v_neurons(4,:));hold on;
% plot(SNN.spikes(3,:));
% ylim([-1 1.5])
% 
figure();
subplot(4,1,1)
imagesc(v_neurons);
subplot(4,1,2)
imagesc(ge_neurons);
subplot(4,1,3)
imagesc(gi_neurons);
subplot(4,1,4)
imagesc(spikes);
colormap(flipud(gray))

% 
figure();
imagesc(spikes);
colormap(flipud(gray))

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