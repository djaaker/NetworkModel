function [fig] = plotSingleNeuronDynamics(select_neuron,SNN)

spike_time = find(SNN.spikes(select_neuron,:) == 1);
v_neuron = SNN.v_neurons(select_neuron,:);
v_neuron(spike_time) = 0;
time_vector = (1:size(v_neuron, 2)) * dt; % Time in seconds
ge_neuron = SNN.ge_neurons(select_neuron,:);
gi_neuron = SNN.gi_neurons(select_neuron,:);

input_spikes = SNN.spikes([SNN.neurons(select_neuron).connections_E;SNN.neurons(select_neuron).connections_I],:);

z_g_neuron = zscore(ge_neuron-gi_neuron);

fig = figure();
subplot(3,1,1)
plot(time_vector,v_neuron,'Color',[0.1 0.1 1],'LineWidth',1.5);hold on;
yline(SNN.neurons(select_neuron).V_thresh);
ylabel('Membrane Voltage (V)');xlabel('Time (s)');
xline(0.2);
ylim([-80e-3 10e-3]);

subplot(3,1,2)
% plot(time_vector,ge_neuron,'Color',[1 0.2 0],'LineWidth',1.5,'LineStyle','-');hold on;
% plot(time_vector,gi_neuron,'Color',[1 0 0.2],'LineWidth',1.5,'LineStyle','-');
plot(time_vector,z_g_neuron,'Color',[1 0 0],'LineWidth',1.5,'LineStyle','-');
ylabel('g_e - g_i (z scored)');xlabel('Time (s)');
yline(0);

subplot(3,1,3)
imagesc(time_vector,1:1:numel(SNN.neurons(select_neuron).connections_E),input_spikes);
colormap(flipud(gray));
axis xy
yline(numel(SNN.neurons(select_neuron).connections_E));
ylabel('Input Spikes');xlabel('Time (s)');
[FR,timeFR] = estimateFiringRate(input_spikes,25e-3,10e3);
yyaxis right
plot(timeFR,smoothdata(mean(FR,1)));
ylabel('Mean Input firing rate (Hz)');

end