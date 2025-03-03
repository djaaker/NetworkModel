classdef Network2Layer
    properties
        grid_size_e             % the size of the grid - excitatory
        grid_size_i             % the size of the grid - inhibitory
        grid_length             % length of one side of the grid
        num_neurons             % Total number of neurons in the network
        N_E                     % Total number of excitatory  neurons
        N_I                     % Total number of inhibitory neurons 
        K_e                     % Number of inhibitory synapses per neuron
        K_i                     % Number of excitatory synapses per neuron
        neurons                 % Array of all neurons in the network
        dt                      % Time step of the simulation
        spikes                  % matrix storing all the spikes
        v_neurons               % matrix storing all the v of neurons
        ge_neurons              % matrix storing all the ge of neurons
        gi_neurons              % matrix storing all the gi of neurons
        total_time              % total time of the simulation
        EI_tag                  % EI tag of all neurons
    end
    
    methods
        function obj = Network2Layer(grid_size,grid_length,excitatory_ratio, sigma,vAP,tau_syn,dt,total_time)
            
            target_neuron_density = 2.8125e10; % per m2
            % Number of E and I neurons
            num_neurons = grid_size^2;
            N_E = (ceil(sqrt(excitatory_ratio* num_neurons)))^2;
            N_I = (ceil(sqrt((1-excitatory_ratio)* num_neurons)))^2;
            grid_size_e = sqrt(N_E);
            grid_size_i = sqrt(N_I);
            num_neurons =  N_E+N_I;
            
            neuron_density = num_neurons/(grid_length^2);
            
            if ( neuron_density/target_neuron_density > 1.2 )
                warning('Network is too dense')
                disp('Sugested length :')
                suggested_grid_length = sqrt(num_neurons/target_neuron_density)
                
            elseif (neuron_density/target_neuron_density <0.8)
                warning('Network is too sparse')
                disp('Sugested length :')
                suggested_grid_length = sqrt(num_neurons/target_neuron_density)
            else
                disp('Network has target neuron density');
            end
            
            if (grid_length > 5*sigma)
                K = 1200;
                K_e = round(0.8*K);
                K_i = round(0.2*K);
                fprintf('The number of synapses from each neuron : %d\n', K);
            else
                K = round(0.0075*num_neurons);
                K_e = round(0.8*K);
                K_i = round(0.2*K);
                K = K_e + K_i;
                fprintf('The number of synapses from each neuron : %d\n', K);
            end
            
            spacing_e =  grid_length/ grid_size_e;
            spacing_i = grid_length/ grid_size_i; 


            % Distribute E and I neurons 
            EI_tag = [zeros(1,N_E), ones(1,N_I)];
            
            % Creating an array of Neurons 
            neurons = Neuron.empty( num_neurons, 0);
            
            % Assigning x-y coordinates to the neurons
            x_origin = (grid_size_e*spacing_e)/2;
            y_origin = (grid_size_e*spacing_e)/2;
            
            [X_e,Y_e] = meshgrid( (1: grid_size_e), (1: grid_size_e));
            X_e = (reshape(X_e,[],1)*spacing_e)-x_origin;
            Y_e = (reshape(Y_e,[],1)*spacing_e)-y_origin;
            
            [X_i,Y_i] = meshgrid( (1: grid_size_i), (1: grid_size_i));
            X_i = (reshape(X_i,[],1)*spacing_i)-x_origin;
            Y_i = (reshape(Y_i,[],1)*spacing_i)-y_origin;
            

            %% For excitatory neurons 
            for i=1:grid_size_e
                for j=1:grid_size_e
                    neuron_idx = (i-1)*grid_size_e+j;        
                    x_coord_neuron = X_e(neuron_idx);
                    y_coord_neuron = Y_e(neuron_idx);
                    
            
                    % Connections to E neurons in network
                    % getting distance from the neuron
                    D = sqrt((X_e-x_coord_neuron).^2 + (Y_e-y_coord_neuron).^2);
                    % getting connection probability distribution from the neuron
                    conn_prob = exp((-D.^2)/(2*(sigma^2)));
                    % removing self connection
                    conn_prob(neuron_idx) = 0;
                    random_matrix = rand(size(conn_prob));
                    conn = random_matrix<conn_prob;
                    delays = tau_syn + ((conn.*D)/vAP);
            
                    connections_E = find(conn == 1);
                    connection_delays_E = round(delays(connections_E)/(dt));
                    [connections_E,randomindx] = randomSampleFromArray(connections_E,K_e);
                    connection_delays_E = connection_delays_E(randomindx);
                    
                    % Connections to I neurons in network
                    % getting distance from the neuron
                    D = sqrt((X_i-x_coord_neuron).^2 + (Y_i-y_coord_neuron).^2);
                    % getting connection probability distribution from the neuron
                    conn_prob = exp((-D.^2)/(2*(sigma^2)));
                    random_matrix = rand(size(conn_prob));
                    conn = random_matrix<conn_prob;
                    delays = tau_syn + ((conn.*D)/vAP);
            
                    connections_I = find(conn == 1);
                    connection_delays_I = round(delays(connections_I)/(dt));
                    [connections_I,randomindx] = randomSampleFromArray(connections_I,K_i);
                    connection_delays_I = connection_delays_I(randomindx);
                    % adjusting for e neurons before i neurons
                    connections_I = connections_I + N_E;
            
                    neurons(neuron_idx) = Neuron(neuron_idx,~EI_tag(neuron_idx),connections_E,connections_I,X_e(neuron_idx),Y_e(neuron_idx),connection_delays_E,connection_delays_I);
                end
            end

            %% For inhibitory neurons 
            for i=1:grid_size_i
                for j=1:grid_size_i
                    neuron_idx = (i-1)*grid_size_i+j;        
                    x_coord_neuron = X_i(neuron_idx);
                    y_coord_neuron = Y_i(neuron_idx);
            
                    % Connections to E neurons in network
                    % getting distance from the neuron
                    D = sqrt((X_e-x_coord_neuron).^2 + (Y_e-y_coord_neuron).^2);
                    % getting connection probability distribution from the neuron
                    conn_prob = exp((-D.^2)/(2*(sigma^2)));
                    random_matrix = rand(size(conn_prob));
                    conn = random_matrix<conn_prob;
                    delays = tau_syn + ((conn.*D)/vAP);
            
                    connections_E = find(conn == 1);
                    connection_delays_E = round(delays(connections_E)/(dt));
                    [connections_E,randomindx] = randomSampleFromArray(connections_E,K_e);
                    connection_delays_E = connection_delays_E(randomindx);
                    
                    % Connections to I neurons in network
                    % getting distance from the neuron
                    D = sqrt((X_i-x_coord_neuron).^2 + (Y_i-y_coord_neuron).^2);
                    % getting connection probability distribution from the neuron
                    conn_prob = exp((-D.^2)/(2*(sigma^2)));
                    % removing self connection
                    conn_prob(neuron_idx) = 0;
                    random_matrix = rand(size(conn_prob));
                    conn = random_matrix<conn_prob;
                    delays = tau_syn + ((conn.*D)/vAP);
            
                    connections_I = find(conn == 1);
                    connection_delays_I = round(delays(connections_I)/(dt));
                    [connections_I,randomindx] = randomSampleFromArray(connections_I,K_i);
                    connection_delays_I = connection_delays_I(randomindx);
                    % adjusting for e neurons before i neurons
                    connections_I = connections_I + N_E;
            
                    neurons(neuron_idx+N_E) = Neuron(neuron_idx+N_E,~EI_tag(neuron_idx+N_E),connections_E,connections_I,X_i(neuron_idx),Y_i(neuron_idx),connection_delays_E,connection_delays_I);
                end
            end


            obj.grid_size_e = grid_size_e;            
            obj.grid_size_i = grid_size_i;             
            obj.grid_length             
            obj.num_neurons = num_neurons;
            obj.neurons = neurons;
            obj.N_E = N_E;
            obj.N_I = N_I;
            obj.K_e = K_e;
            obj.K_i = K_i;
            obj.dt = dt;
            n_time = round(obj.total_time/obj.dt);
            obj.spikes = false(obj.num_neurons,n_time);
            obj.v_neurons = zeros(obj.num_neurons,n_time);
            obj.ge_neurons = zeros(obj.num_neurons,n_time);
            obj.gi_neurons = zeros(obj.num_neurons,n_time);      
            obj.total_time = total_time;     
            obj.EI_tag = EI_tag;             

        end

        function obj = plot_network(obj)
            figure();hold on;
            for i=1:obj.num_neurons
                if obj.neurons(i).is_excitatory == 1
                    color = 'b';
                else
                    color = 'r';
                end
                scatter(obj.neurons(i).x_coord,obj.neurons(i).y_coord,10,color);
            end 
        end

        function obj = plot_neuron_connection(obj,neuron_idx)
            A = zeros(obj.num_neurons,obj.num_neurons);
            A(neuron_idx,:) = obj.neurons(neuron_idx).connections_E + -1*obj.neurons(neuron_idx).connections_I;
            g = digraph(A);
            [X,Y] = meshgrid( (1:obj.grid_size), (1:obj.grid_size));
            X = reshape(X,[],1);
            Y = reshape(Y,[],1);
            figure();
            plot(g,'XData',X,'YData',Y);
        end

        function obj = update(obj, time_indx, time, I_ext,spike_ext)
            num_neurons_tmp = obj.num_neurons;  % Avoid repeated property access
        
            % Preallocate arrays for neuron properties
            spikes_tmp = false(num_neurons_tmp, 1); 
            v_tmp = zeros(num_neurons_tmp, 1);
            ge_tmp = zeros(num_neurons_tmp, 1);
            gi_tmp = zeros(num_neurons_tmp, 1);
            
            % Extract necessary properties to reduce handle object overhead
            dt = obj.dt;
            neuron_data = obj.neurons; % Store handle references (not used in `parfor` directly)
        
            % Parallelized neuron updates
            for i = 1:num_neurons_tmp
                neuron = neuron_data(i); % Copy handle object
                [neuron, spikes_tmp(i)] = neuron.update(I_ext(i),spike_ext(i), dt, time);
                v_tmp(i) = neuron.V; % Store new voltage
                ge_tmp(i) = neuron.ge; % Store new ge
                gi_tmp(i) = neuron.gi; % Store new gi
                neuron_data(i) = neuron; % Store modified neuron
            end
        
            % Assign back computed results
            obj.spikes(:, time_indx) = spikes_tmp;
            obj.v_neurons(:, time_indx) = v_tmp;
            obj.ge_neurons(:, time_indx) = ge_tmp;
            obj.gi_neurons(:, time_indx) = gi_tmp;
            
        
            % Serial update of spikes in all neurons
            spikes_current = squeeze(obj.spikes(:, time_indx)); 
            for i = 1:num_neurons_tmp
                neuron_data(i) = neuron_data(i).update_spikes_buff(spikes_current);
            end

            obj.neurons = neuron_data; % Assign updated neurons back
        end

        function obj = update_par(obj, time_indx, time, I_ext,spike_ext)
            num_neurons_tmp = obj.num_neurons;  % Avoid repeated property access
        
            % Preallocate arrays for neuron properties
            spikes_tmp = false(num_neurons_tmp, 1); 
            v_tmp = zeros(num_neurons_tmp, 1);
            ge_tmp = zeros(num_neurons_tmp, 1);
            gi_tmp = zeros(num_neurons_tmp, 1);
            
            % Extract necessary properties to reduce handle object overhead
            dt = obj.dt;
            neuron_data = obj.neurons; % Store handle references (not used in `parfor` directly)
        
            % Parallelized neuron updates
            parfor i = 1:num_neurons_tmp
                neuron = neuron_data(i); % Copy handle object
                [neuron, spikes_tmp(i)] = neuron.update(I_ext(i),spike_ext(i), dt, time);
                v_tmp(i) = neuron.V; % Store new voltage
                ge_tmp(i) = neuron.ge; % Store new ge
                gi_tmp(i) = neuron.gi; % Store new gi
                neuron_data(i) = neuron; % Store modified neuron
            end
        
            % Assign back computed results
            obj.spikes(:, time_indx) = spikes_tmp;
            obj.v_neurons(:, time_indx) = v_tmp;
            obj.ge_neurons(:, time_indx) = ge_tmp;
            obj.gi_neurons(:, time_indx) = gi_tmp;
            
        
            % Serial update of spikes in all neurons
            spikes_current = squeeze(obj.spikes(:, time_indx)); 
            parfor i = 1:num_neurons_tmp
                neuron_data(i) = neuron_data(i).update_spikes_buff(spikes_current);
            end

            obj.neurons = neuron_data; % Assign updated neurons back
        end

    end
end






