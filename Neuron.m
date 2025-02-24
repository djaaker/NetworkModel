classdef Neuron
    properties
        id                          % the id of the neuron 1 <= id <= N where N is the total number of neurons
        is_excitatory               % True if excitatory, false if inhibitory
        connections_E               % Connections to other E neurons 1 - if connected, 2 - not connected 
        connections_I               % Connections to other I neurons 1 - if connected, 2 - not connected 
        x_coord                     % x coordinate 
        y_coord                     % y coordinate 
        connection_delays_E         % Delays based on connection length for E neurons
        connection_delays_I         % Delays based on connection length for I neurons

        V                           % Membrane potential
        E_L = -65e-3;               % Resting potential (mV)
        V_thresh = -50e-3;          % Threshold potential (mV)
        V_reset = -70e-3;           % Reset potential (mV)
        refractory_period = 5e-3;   % Refractory period (ms)
        last_spike_time = -inf;     % Last spike time (ms)
        E_e = 0;                    % excitatory reversal potential
        E_i = -80e-3;               % inhibitory reversal potential
        Cm = 200e-12                % membrance capacitance

        ge                          % excitatory conductance
        gi                          % inhibitory conductance

        GL = 10e-9;                 % leak conductance
        GE = 1e-9;                  % excitatory synaptic weight
        GI = 10e-9;                 % inhibitory synaptic weight

        tau_syn_e = 5e-3;           % excitatory synpatic time constant
        tau_syn_i = 5e-3;           % inhibitory synpatic time constant
    end
    
    methods
        function obj = Neuron(id,is_excitatory,connections_E,connections_I,x_coord,y_coord,connection_delays_E,connection_delays_I)

            obj.id = id;
            obj.is_excitatory = is_excitatory;
            obj.connections_E = connections_E;
            obj.connections_I = connections_I;
            obj.x_coord = x_coord;
            obj.y_coord = y_coord;
            obj.connection_delays_E = connection_delays_E;
            obj.connection_delays_I = connection_delays_I;

            obj.V = obj.E_L;
            obj.ge = 0;
            obj.gi = 0;
        end
        
        function [obj, spike] = update(obj,spikes,I_ext,dt,time)
            time_indx = round(time/dt);
            % For excitatory conductance
            spikes_delay = circshiftSpikeDelay(spikes,obj.connection_delays_E');
            try
                delta_ge = sum(squeeze(spikes_delay(:,time_indx)))*obj.GE;
            catch
                warning('breakpoint');
                end
            obj.ge = obj.ge - obj.ge/obj.tau_syn_e + delta_ge;

            % For inhibitory conductance
            spikes_delay = circshiftSpikeDelay(spikes,obj.connection_delays_I');
            delta_gi = sum(squeeze(spikes_delay(:,time_indx)))*obj.GI;
            obj.gi = obj.gi - obj.gi/obj.tau_syn_i + delta_gi;

            if (time - obj.last_spike_time) < obj.refractory_period
                obj.V = obj.V_reset; % Refractory period
                spike = false;
                return;
            end
            obj.V = obj.V + dt*(1/obj.Cm)*(obj.GL*(obj.E_L - obj.V) ... 
                + obj.ge*(obj.E_e - obj.V) ...
                + obj.gi*(obj.E_i - obj.V) ...
                + I_ext);
            spike = obj.V >= obj.V_thresh;
            if spike
                obj.V = obj.V_reset;
                obj.last_spike_time = time;
            end
        end
    end
end

