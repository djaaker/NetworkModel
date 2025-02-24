classdef Network_old
    properties
        neurons
        num_neurons
        grid_size
        connectivity_matrix
    end
    
    methods
        function obj = Network_old(grid_size, excitatory_ratio, connection_prob)
            obj.grid_size = grid_size;
            obj.num_neurons = grid_size^2;
            obj.neurons = Neuron.empty(obj.num_neurons, 0);
            
            for i = 1:obj.num_neurons
                is_excitatory = rand() < excitatory_ratio;
                obj.neurons(i) = Neuron(is_excitatory, obj.num_neurons);
            end
            obj = obj.initialize_connections(connection_prob);
        end
        
        function obj = initialize_connections(obj, connection_prob)
            obj.connectivity_matrix = zeros(obj.num_neurons);
            
            for i = 1:obj.num_neurons
                for j = 1:obj.num_neurons
                    if i ~= j
                        distance = obj.compute_distance(i, j);
                        prob = exp(-distance^2 / connection_prob^2);
                        if rand() < prob
                            obj.connectivity_matrix(i, j) = 1;
                        end
                    end
                end
            end
        end
        
        function distance = compute_distance(obj, i, j)
            [xi, yi] = ind2sub([obj.grid_size, obj.grid_size], i);
            [xj, yj] = ind2sub([obj.grid_size, obj.grid_size], j);
            distance = sqrt((xi - xj)^2 + (yi - yj)^2);
        end
        
        function obj = update(obj, dt, time)
            spikes = false(1, obj.num_neurons);
            
            for i = 1:obj.num_neurons
                I = sum(obj.connectivity_matrix(:, i) .* [obj.neurons.synaptic_weights]);
                [obj.neurons(i), spikes(i)] = obj.neurons(i).update(I, dt, time);
            end
            obj = obj.apply_stdp(spikes);
        end
        
        function obj = apply_stdp(obj, spikes)
            for i = 1:obj.num_neurons
                for j = 1:obj.num_neurons
                    if obj.connectivity_matrix(i, j) && spikes(i)
                        obj.neurons(j).synaptic_weights(i) = obj.neurons(j).synaptic_weights(i) + 0.01;
                    elseif obj.connectivity_matrix(i, j) && spikes(j)
                        obj.neurons(i).synaptic_weights(j) = obj.neurons(i).synaptic_weights(j) - 0.01;
                    end
                end
            end
        end
    end
end