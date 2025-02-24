classdef Network
    properties
        grid_size               % the size of the grid to be simulated in number of neurons
        grid_length             % length of one side of the grid
        num_neurons             % Total number of neurons in the network
        neurons                 % Array of all neurons in the network
        dt                      % Time step of the simulation
        spikes                  % matrix storing all the spikes
        v_neurons               % matrix 
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

            n_time = obj.total_time/obj.dt;

            obj.spikes = false(obj.num_neurons,n_time);
            obj.v_neurons = zeros(obj.num_neurons,n_time);

            % Number of E and I neurons
            N_E = round(excitatory_ratio*obj.num_neurons);
            N_I = obj.num_neurons-N_E;
            spacing = obj.grid_length/obj.grid_size; % in mm

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
                    conn_prob = exp(-D.^2/(2*sigma^2));
                    % removing self connection
                    conn_prob(i,j) = 0;conn_prob(j,i) = 0;
                    random_matrix = rand(size(conn_prob));
                    conn = random_matrix<conn_prob;
                    delays = (conn.*D)/vAP;
                    delays = reshape(delays,[],1);
                    conn = reshape(conn,[],1);
                    connections_E = conn.*(~EI_tag)';
                    connections_I = conn.*(EI_tag)';
                    connection_delays_E = round(((delays+tau_syn).*connections_E)/(obj.dt));
                    connection_delays_I = round(((delays+tau_syn).*connections_I)/(obj.dt));
                    % synaptic_weights = rand(obj.num_neurons,1)*0.05;
                    % synaptic_weights_E = synaptic_weights.*connections_E;
                    % synaptic_weights_I = synaptic_weights.*connections_I;
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

        function obj = update(obj,dt,time,I_ext)
            for i=1:obj.num_neurons
                time_indx = round(time/dt);
                [obj.neurons(i), obj.spikes(i,time_indx)] = obj.neurons(i).update(obj.spikes,I_ext(i),obj.dt,time);
                obj.v_neurons(i,time_indx) = obj.neurons(i).V;
            end
        end
    end
end






