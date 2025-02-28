
target_neuron_density = 2.8125e10; % per m2
excitatory_ratio = 0.8;
grid_size = 100;
grid_length = 0.6e-3;
% Decay constant for connectivity 
sigma = 0.3e-3; %% in m 
% Speed of AP 
vAP = 0.2; % in m/s
% Synaptic delay 
tau_syn = 300e-6; % in seconds


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

X = [X_e;X_i];
Y = [Y_e;Y_i];




dummy_connection_matrix_e = zeros(K_e,1);
dummy_connection_matrix_i = zeros(K_i,1);

for id = 1: num_neurons
    neurons(id) = Neuron(id,~EI_tag(id),dummy_connection_matrix_e,dummy_connection_matrix_e,X(id),Y(id),dummy_connection_matrix_e,dummy_connection_matrix_e);
end

