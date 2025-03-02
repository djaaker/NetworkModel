function firingRates = estimateFiringRate(spikingData, binDuration, samplingRate)
% estimateFiringRate - Estimate firing rates by binning spiking data.
%
% Syntax:  firingRates = estimateFiringRate(spikingData, binDuration, samplingRate)
%
% Inputs:
%    spikingData  - Binary matrix (neurons x samples) where 1 indicates a spike.
%    binDuration  - Duration of each bin in seconds (default: 0.05, i.e., 50 ms).
%    samplingRate - Sampling rate in Hz (default: 10000 Hz).
%
% Outputs:
%    firingRates - Matrix (neurons x numBins) containing the firing rate (Hz)
%                  for each neuron in each bin.

    % Set default parameters if not provided
    if nargin < 2
        binDuration = 0.05; % 50 ms default bin duration
    end
    if nargin < 3
        samplingRate = 10000; % default sampling rate 10 kHz
    end

    % Compute number of samples per bin
    samplesPerBin = floor(binDuration * samplingRate);

    % Get the dimensions of spikingData
    [numNeurons, numSamples] = size(spikingData);

    % Determine number of complete bins
    numBins = floor(numSamples / samplesPerBin);

    % Preallocate firingRates matrix
    firingRates = zeros(numNeurons, numBins);

    % Loop through each bin, sum the spikes, and convert count to Hz
    for bin = 1:numBins
        idxStart = (bin - 1) * samplesPerBin + 1;
        idxEnd = bin * samplesPerBin;
        % Sum spikes in each bin for each neuron
        binSpikeCount = sum(spikingData(:, idxStart:idxEnd), 2);
        % Convert spike count to firing rate (spikes per second)
        firingRates(:, bin) = binSpikeCount / binDuration;
    end
end
