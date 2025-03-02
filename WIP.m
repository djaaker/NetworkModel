
target_neuron_density = 2.8125e10; % per m2
excitatory_ratio = 0.8;
grid_size = 50;
grid_length = 0.3e-3;
% Decay constant for connectivity 
sigma = 0.3e-3; %% in m 
% Speed of AP 
vAP = 0.2; % in m/s
% Synaptic delay 
tau_syn = 300e-6; % in seconds

dt = 0.1e-3;


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

if (grid_length > 2.6*sigma)
    K = 1200;
    K_e = round(0.8*K);
    K_i = round(0.2*K);
    fprintf('The number of synapses from each neuron : %d\n', K);
else
    K = round(0.005*num_neurons);
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
Y_e = (reshape(Y_e,[],1)*spacing_e)-x_origin;

[X_i,Y_i] = meshgrid( (1: grid_size_i), (1: grid_size_i));
X_i = (reshape(X_i,[],1)*spacing_i)-x_origin;
Y_i = (reshape(Y_i,[],1)*spacing_i)-x_origin;

% X = [X_e;X_i];
% Y = [Y_e;Y_i];

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

% 
% while num_synapses > 0
%     distance_draw = randn * sigma;
%     [indx,actual_distance] = find_closest(id,distance_draw,X_e,Y_e);
%     if indx == -1
%         continue;
%     else
%         connection_matrix_e(num_synapses) = indx;
%         delay_matrix_e(num_synapses) = round((tau_syn + (actual_distance/vAP))/dt);
%     end
%     num_synapses = num_synapses - 1;
% end
% 
% num_synapses = K_i;
% 
% while num_synapses > 0
%     distance_draw = randn * sigma;
%     [indx,actual_distance] = find_closest(id,distance_draw,X_i,Y_i);
%     if indx == -1
%         continue;
%     else
%         connection_matrix_i(num_synapses) = N_E + indx;
%         delay_matrix_i(num_synapses) = round((tau_syn + (actual_distance/vAP))/dt);
%     end
%     num_synapses = num_synapses - 1;
% end