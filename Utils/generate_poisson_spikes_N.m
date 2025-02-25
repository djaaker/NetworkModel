function spike_trains = generate_poisson_spikes_N(N, rate, duration, dt)
    % GENERATE_POISSON_SPIKES_N Generates Poisson spike trains for N neurons
    %
    % Inputs:
    %   N        - Number of neurons
    %   rate     - Firing rate (Hz) [scalar or N-length vector]
    %   duration - Total duration of spike train (seconds)
    %   dt       - Time step (seconds)
    %
    % Output:
    %   spike_trains - Binary matrix (N x num_steps), where each row is a neuron

    % Time vector
    time_steps = 0:dt:duration-dt; 
    num_steps = length(time_steps);

    % Ensure rate is a vector of size N
    if isscalar(rate)
        rate = rate * ones(N, 1); 
    end

    % Compute spike probability for each neuron
    spike_prob = rate * dt; % Probability of spiking at each time step

    % Generate random spikes for each neuron
    spike_trains = rand(N, num_steps) < spike_prob;
end
