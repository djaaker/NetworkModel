function [SNN,spikes,v_neurons] = updateNetworkGPU(time_indx,time,I_ext,spike_ext,SNN,spikes,v_neurons,num_neurons,dt)
            
% Parallelized neuron updates
parfor i = 1:num_neurons
    [SNN(i), spikes(i, time_indx)] = SNN(i).update(I_ext(i),spike_ext(i), dt, time);
    v_neurons(i, time_indx) = SNN(i).V; % Store new voltage
end

% Serial update of spikes in all neurons
spikes_current = squeeze(spikes(:, time_indx)); 
parfor i = 1:num_neurons
    SNN(i).update_spikes_buff(spikes_current);
end

end

