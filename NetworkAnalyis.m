%% Plotting and generating video for spiking activity
plot_spiking_activity_video(SNN,dt,saveFlag,filename);

%% Generating LFP
LFP = calculateLFP(SNN);

% Plotting LFP, spikes and input 
plotLFP(LFP,SNN,spike_train_ext);
% Genreating a video of the LFP                 
generateLFPVideo(LFP.LFP_downsampled, 1/LFP.target_sampling_freq, lfpfilename);
generateLFPVideo(LFP.LFP_gaussian_downsampled, 1/LFP.target_sampling_freq, lfpfilename2);

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
parameters.xspacing = LFP.lfp_patch_size*(SNN.grid_length/SNN.grid_size_e) ;parameters.yspacing = LFP.lfp_patch_size*(SNN.grid_length/SNN.grid_size_e);
parameters.rhoThres= 0.5;
% Wave detection

allwaves.LFPIndex = (1:1:size(LFP.LFP_combined,3))';
xf{1,1} = LFP.xf;
xgp{1,1} = LFP.xgp;
wt{1,1} = LFP.wt;
Wavesall = detectWaves(xf,xgp,wt,allwaves,parameters,parameters.rhoThres);

%% Plotting
plotSpikingNetwork(SNN)

firingRates = estimateFiringRate(SNN.spikes,20e-3,1/dt);
%%

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