classdef Network
    properties
        grid_size               % the size of the grid to be simulated in number of neurons
        grid_length             % length of one side of the grid
        num_neurons             % Total number of neurons in the network
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
        function obj = Network(grid_size,grid_length,excitatory_ratio, sigma,vAP,tau_syn,time_step,total_time)
            
            obj.dt = time_step;
            obj.total_time = total_time;
            
            obj.grid_size = grid_size;
            obj.grid_length = grid_length;
            obj.num_neurons =  grid_size*grid_size;

            n_time = round(obj.total_time/obj.dt);

            obj.spikes = false(obj.num_neurons,n_time);
            obj.v_neurons = zeros(obj.num_neurons,n_time);
            obj.ge_neurons = zeros(obj.num_neurons,n_time);
            obj.gi_neurons = zeros(obj.num_neurons,n_time);

            % Number of E and I neurons
            N_E = round(excitatory_ratio*obj.num_neurons);
            N_I = obj.num_neurons-N_E;
            spacing = obj.grid_length/obj.grid_size; % in m

            % Distribute E and I neurons 
            EI_tag = [ones(1,N_I), zeros(1,N_E)];
            EI_tag = EI_tag(randperm(length(EI_tag)));
            obj.EI_tag = EI_tag;

            % Creating an array of Neurons 
            obj.neurons = Neuron.empty(obj.num_neurons, 0);
            
            f = waitbar(0,'Generating network connections');
            for i=1:obj.grid_size
                for j=1:obj.grid_size
                    % getting distance from the neuron
                    [X,Y] = meshgrid( (1:obj.grid_size)-i, (1:obj.grid_size)-j);
                    D = sqrt( (X*spacing).^2 + (Y*spacing).^2 );
                    % getting connection probability distribution from the neuron
                    conn_prob = exp((-D.^2)/(2*(sigma^2)));
                    % removing self connection
                    conn_prob(i,j) = 0;conn_prob(j,i) = 0;
                    random_matrix = rand(size(conn_prob));
                    conn = random_matrix<conn_prob;
                    delays = (conn.*D)/vAP;
                    delays = reshape(delays,[],1);
                    conn = reshape(conn,[],1);
                    connections_E = find(conn.*(~EI_tag)' == 1);
                    connections_I = find(conn.*(EI_tag)' == 1);
                    connection_delays_E = round((tau_syn*ones(size(connections_E)) + delays(connections_E))/(obj.dt));
                    connection_delays_I = round((tau_syn*ones(size(connections_I)) + delays(connections_I))/(obj.dt));
                    id = (i-1)*obj.grid_size+j;
                    obj.neurons(id) = Neuron(id,~EI_tag(id),connections_E,connections_I,i,j,connection_delays_E,connection_delays_I);
                end
                waitbar(i/obj.grid_size,f,sprintf('Generating network connections: %d %%', floor(i/obj.grid_size*100)));
                pause(0.01)
            end
            clear X Y D conn conn_prob delays random_matrix
            close (f);
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

        function obj = update(obj,time_indx,time,I_ext,spike_ext)
            % Getting spikes from all neurons
            for i=1:obj.num_neurons
                [obj.neurons(i), obj.spikes(i,time_indx)] = obj.neurons(i).update(I_ext(i),spike_ext(i),obj.dt,time);
                obj.v_neurons(i,time_indx) = obj.neurons(i).V;
            end
            % updating spikes in all neurons 
            for i=1:obj.num_neurons
                obj.neurons(i).update_spikes_buff(squeeze(obj.spikes(:,time_indx)));
            end
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
            obj.neurons = neuron_data; % Assign updated neurons back
        
            % Serial update of spikes in all neurons
            spikes_current = squeeze(obj.spikes(:, time_indx)); 
            for i = 1:num_neurons_tmp
                obj.neurons(i) = obj.neurons(i).update_spikes_buff(spikes_current);
            end
        end

    end
end






