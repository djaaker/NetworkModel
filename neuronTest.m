% Simulate a single neuron
single_neuron = Neuron(true, 1);
dt = 1; % Time step (ms)
total_time = 1000; % 1 second (1000 ms)
time_values = 0:dt:total_time;
voltage_trace = zeros(size(time_values));

for t = 1:length(time_values)
    time = time_values(t);
    I = 100; % Constant input current
    [single_neuron, ~] = single_neuron.update(I, dt, time);
    voltage_trace(t) = single_neuron.V;
end

figure;
plot(time_values, voltage_trace);
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
title('Single Neuron Simulation');
grid on;
