function [fig] = plotSingleNeuronDynamics(select_neuron,SNN,spikes,v_neurons,ge_neurons,gi_neurons)

spike_time = find(spikes(select_neuron,:) == 1);
v_neuron = v_neurons(select_neuron,:);
v_neuron(spike_time) = 0;
time_vector = (1:size(v_neuron, 2)) * SNN.dt; % Time in seconds
ge_neuron = ge_neurons(select_neuron,:);
gi_neuron = gi_neurons(select_neuron,:);

input_spikes =spikes([SNN.neurons{select_neuron}.connections_E;SNN.neurons{select_neuron}.connections_I],:);

z_g_neuron = zscore(ge_neuron-gi_neuron);

fig = figure();
subplot(3,1,1)
plt1 = plot(time_vector,v_neuron,'Color',[0.1 0.1 1],'LineWidth',1.5);hold on;
yline(SNN.neurons{select_neuron}.V_thresh);
ylabel('Membrane Voltage (V)');xlabel('Time (s)');
xline(0.2);
ylim([-80e-3 10e-3]);

subplot(3,1,2)
% plot(time_vector,ge_neuron,'Color',[1 0.2 0],'LineWidth',1.5,'LineStyle','-');hold on;
% plot(time_vector,gi_neuron,'Color',[1 0 0.2],'LineWidth',1.5,'LineStyle','-');
plt2 = plot(time_vector,z_g_neuron,'Color',[1 0 0],'LineWidth',1.5,'LineStyle','-');
ylabel('g_e - g_i (z scored)');xlabel('Time (s)');
yline(0);
ylim([-4 4]);

subplot(3,1,3)
plt3 = imagesc(time_vector,1:1:numel(SNN.neurons{select_neuron}.connections_E),input_spikes);
colormap(flipud(gray));
axis xy
yline(numel(SNN.neurons{select_neuron}.connections_E));
ylabel('Input Spikes');xlabel('Time (s)');
[FR,timeFR] = estimateFiringRate(input_spikes,25e-3,10e3);
yyaxis right
plot(timeFR,smoothdata(mean(FR,1)));
ylabel('Mean Input firing rate (Hz)');

linkaxes([plt1 plt2 plt3],'x');

%% Plotting



end