classdef Neuron_CPP_block
    %#codegen
    properties
        id = -1;                         % the id of the neuron 1 <= id <= N where N is the total number of neurons
        is_excitatory = true;              % True if excitatory, false if inhibitory
        connections_E  = [];              % Connections to other E neurons  - id of the connected  
        connections_I  = [];             % Connections to other I neurons - id of the connected  
        x_coord = 0;                    % x coordinate 
        y_coord = 0;                   % y coordinate 
        connection_delays_E = [];         % Delays based on connection length for E neurons
        connection_delays_I = [];       % Delays based on connection length for I neurons

        V = -65e-3;                          % Membrane potential
        E_L = -65e-3;               % Resting potential (mV)
        V_thresh = -50e-3;          % Threshold potential (mV)
        V_reset = -70e-3;           % Reset potential (mV)
        refractory_period = 5e-3;   % Refractory period (ms)
        last_spike_time = -inf;     % Last spike time (ms)
        E_e = 0;                    % excitatory reversal potential
        E_i = -80e-3;               % inhibitory reversal potential
        Cm = 200e-12                % membrance capacitance

        ge  = 0;                        % excitatory conductance
        gi  = 0;                        % inhibitory conductance

        GL = 10e-9;                 % leak conductance
        GE = 3e-9;                  % excitatory synaptic weight
        GI = 30e-9;                 % inhibitory synaptic weight

        tau_syn_e = 5e-3;           % excitatory synpatic time constant
        tau_syn_i = 5e-3;           % inhibitory synpatic time constant

        spikes_e_buff = FIFOSpikes(0, 1);             % Buffer for containing all past spikes from excitatory connections 
        spikes_i_buff = FIFOSpikes(0, 1);             % Buffer for containing all past spikes from inhibitory connections

    end
    
    methods

        function obj = Neuron_CPP_block(id,is_excitatory,connections_E,connections_I,x_coord,y_coord,connection_delays_E,connection_delays_I,GE,GI)
            if nargin > 0
                obj.id = id;
                obj.is_excitatory = is_excitatory;
                obj.connections_E = connections_E;
                obj.connections_I = connections_I;
                obj.x_coord = x_coord;
                obj.y_coord = y_coord;
                obj.connection_delays_E = connection_delays_E;
                obj.connection_delays_I = connection_delays_I;
                obj.GE = GE;
                obj.GI = GI;
                
                % Define voltage range in volts
                minVoltage = -70e-3; % -70 mV in V
                maxVoltage = 0;      % 0 V
                
                % Generate a random voltage within the range
                obj.V = minVoltage + (maxVoltage - minVoltage) * rand;
                obj.ge = 0;
                obj.gi = 0;
    
                max_e_delay = max(obj.connection_delays_E);
                max_i_delay = max(obj.connection_delays_I);
    
                obj.spikes_e_buff = FIFOSpikes(numel(connection_delays_E),max_e_delay);
                obj.spikes_i_buff = FIFOSpikes(numel(connection_delays_I),max_i_delay);
            else
                obj.id = -1;
                obj.is_excitatory = true;
                obj.connections_E = [];
                obj.connections_I = [];
                obj.x_coord = 0;
                obj.y_coord = 0;
                obj.connection_delays_E = [];
                obj.connection_delays_I = [];
                obj.GE = 0;
                obj.GI = 0;
               
                obj.V = -65e-3;
                obj.ge = 0;
                obj.gi = 0;
    
                obj.spikes_e_buff = FIFOSpikes(0, 1);
                obj.spikes_i_buff = FIFOSpikes(0, 1);
            end
        end
        
        function [obj, spike] = update(obj,I_ext,spike_ext,dt,time)
            % For excitatory conductance
            spikes_delay = obj.spikes_e_buff.get_behind(obj.connection_delays_E+1);
            delta_ge = (sum(spikes_delay,'all')+spike_ext)*obj.GE;
            obj.ge = obj.ge - (dt*obj.ge)/obj.tau_syn_e + delta_ge;

            % For inhibitory conductance
            spikes_delay = obj.spikes_i_buff.get_behind(obj.connection_delays_I+1);
            delta_gi = sum(spikes_delay,'all')*obj.GI;
            obj.gi = obj.gi - (dt*obj.gi)/obj.tau_syn_i + delta_gi;

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
        
        function obj = update_spikes_buff(obj,last_spikes)
            obj.spikes_e_buff = obj.spikes_e_buff.push(last_spikes(obj.connections_E));
            obj.spikes_i_buff = obj.spikes_i_buff.push(last_spikes(obj.connections_I));
        end

    end
end

