function [SNN,spikes,v_neurons,parameters] = generateNetwork(parameters)

parameters.num_neurons = parameters.grid_size^2;
parameters.n_time = round(parameters.total_time/parameters.dt);

spikes = gpuArray(false(parameters.num_neurons,parameters.n_time));
v_neurons = gpuArray(zeros(parameters.num_neurons,parameters.n_time));

% Number of E and I neurons
parameters.N_E = round(parameters.excitatory_ratio*parameters.num_neurons);
parameters.N_I = parameters.num_neurons-parameters.N_E;
parameters.spacing = parameters.grid_length/parameters.grid_size; % in m

% Distribute E and I neurons 
parameters.EI_tag = [ones(1,parameters.N_I), zeros(1,parameters.N_E)];
parameters.EI_tag = parameters.EI_tag(randperm(length(parameters.EI_tag)));
parameters.EI_tag = parameters.EI_tag;


% Creating an array of Neurons 
SNN = NeuronGPU.empty(parameters.num_neurons, 0);

f = waitbar(0,'Generating network connections');
for i=1:parameters.grid_size
    for j=1:parameters.grid_size
        % getting distance from the neuron
        [X,Y] = meshgrid( (1:parameters.grid_size)-i, (1:parameters.grid_size)-j);
        parameters.D = sqrt( (X*parameters.spacing).^2 + (Y*parameters.spacing).^2 );
        % getting connection probability distribution from the neuron
        conn_prob = exp(-parameters.D.^2/(2*parameters.sigma^2));
        % removing self connection
        conn_prob(i,j) = 0;conn_prob(j,i) = 0;
        random_matrix = rand(size(conn_prob));
        conn = random_matrix<conn_prob;
        delays = (conn.*parameters.D)/parameters.vAP;
        delays = reshape(delays,[],1);
        conn = reshape(conn,[],1);
        connections_E = find(conn.*(~parameters.EI_tag)' == 1);
        connections_I = find(conn.*(parameters.EI_tag)' == 1);
        connection_delays_E = round((parameters.tau_syn*ones(size(connections_E)) + delays(connections_E))/(parameters.dt));
        connection_delays_I = round((parameters.tau_syn*ones(size(connections_I)) + delays(connections_I))/(parameters.dt));
        id = (i-1)*parameters.grid_size+j;
        SNN(id) = NeuronGPU(id,~parameters.EI_tag(id),connections_E,connections_I,i,j,connection_delays_E,connection_delays_I);
    end
    waitbar(i/parameters.grid_size,f,sprintf('Generating network connections: %d %%', floor(i/parameters.grid_size*100)));
    pause(0.01)
end
clear X Y D conn conn_prob delays random_matrix
close (f);

end

