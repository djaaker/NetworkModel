function [SNN] =  NetworkSimulation2Layer_CPP_block(grid_size,excitatory_ratio,grid_length,sigma,vAP,tau_syn,kden,GE,GI,time_step,total_time,spike_train_ext,folderpath)
    %#codegen
    % n_time = round(total_time/time_step);
    % Create spiking neural network SNN 

    size_stim = size(spike_train_ext,2);


    spikes_binfile = [folderpath,'/','spikes'];
    v_neurons_binfile = [folderpath,'/','v_neurons'];
    gi_neurons_binfile = [folderpath,'/','gi_neurons'];
    ge_neurons_binfile = [folderpath,'/','ge_neurons'];
    snn_matfile = [folderpath,'/','SNN.mat'];


    
    time_block = 0.1/time_step; % the number of elements in each time block;
    
    % Initializing the network
    disp('Creating a SNN ...');
    SNN = Network2Layer_CPP_block(grid_size,grid_length,excitatory_ratio, sigma,vAP,tau_syn,time_step,total_time,GE,GI,kden);
    
    disp('Saving the SNN ...');
    save(snn_matfile,"SNN","-v7.3");
    disp('Complete.')
    
    num_neurons = SNN.num_neurons;

    % Initializing the spikes,v_neurons,g_e and g_i
    spikes = FIFOSpikes(num_neurons,time_block);
    v_neurons = FIFO(num_neurons,time_block);
    ge_neurons = FIFO(num_neurons,time_block);
    gi_neurons = FIFO(num_neurons,time_block);

    % Preallocating
    % spikes_current = false(num_neurons,1);
    % v_current = zeros(num_neurons,1);
    % ge_current = zeros(num_neurons,1);
    % gi_current = zeros(num_neurons,1);

    
    % Simulating network
    total_iterations = total_time / time_step;
    start_time = tic; % Start timer
    % Simulating network
    
    for i = 1:total_iterations
        I_ext = zeros(num_neurons, 1);
        if i <= size_stim
            spike_ext = (squeeze(spike_train_ext(:,i)))&(~SNN.EI_tag)';
        else
            spike_ext = false(num_neurons,1);
        end
    
        % Update network
        [SNN,spikes_current,v_current,ge_current,gi_current] = SNN.update(i*time_step, I_ext,spike_ext);
        
        spikes = spikes.push(spikes_current);
        v_neurons = v_neurons.push(v_current);
        ge_neurons = ge_neurons.push(ge_current);
        gi_neurons = gi_neurons.push(gi_current);

        if spikes.isFull == 1
            disp('Buffers full. Saving to bin files ...')
            % Saving the buffers 
            if i<=time_block
                save_logical_matrix_bin(spikes_binfile, spikes.buffer, 'w');
                save_matrix_bin(v_neurons_binfile,v_neurons.buffer,'w');
                save_matrix_bin(ge_neurons_binfile,ge_neurons.buffer,'w');
                save_matrix_bin(gi_neurons_binfile,gi_neurons.buffer,'w');
            else
                save_logical_matrix_bin(spikes_binfile, spikes.buffer, 'a');
                save_matrix_bin(v_neurons_binfile,v_neurons.buffer,'a');
                save_matrix_bin(ge_neurons_binfile,ge_neurons.buffer,'a');
                save_matrix_bin(gi_neurons_binfile,gi_neurons.buffer,'a');
            end
            % clearing the buffers 
            spikes = spikes.clearBuffer();
            v_neurons = v_neurons.clearBuffer();
            ge_neurons = ge_neurons.clearBuffer();
            gi_neurons = gi_neurons.clearBuffer();
            disp('Done saving. Cleared buffers.');
        end

        if (i==total_iterations && spikes.head ~= spikes.tail) % Last iterration and the buffer has some elements 
            save_logical_matrix_bin(spikes_binfile, spikes.buffer, 'a');
            save_matrix_bin(v_neurons_binfile,v_neurons.buffer,'a');
            save_matrix_bin(ge_neurons_binfile,ge_neurons.buffer,'a');
            save_matrix_bin(gi_neurons_binfile,gi_neurons.buffer,'a');
        end
       
        % Calculate elapsed time and estimate time remaining
        elapsed_time = toc(start_time);avg_time_per_iter = elapsed_time / i;estimated_time_left = (avg_time_per_iter * (total_iterations - i))/3600;
        disp(estimated_time_left);
        sprintf('SNN Simulation Progress: %.2f%% | Time left: %.2f hrs',(i / total_iterations) * 100, estimated_time_left);

    end

end