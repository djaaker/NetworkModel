function [randomSample,randomIndices]= randomSampleFromArray(sample, n)
% randomSampleFromArray - Randomly samples n elements from a given array.
%
% Syntax: randomSample = randomSampleFromArray(sample, n)
%
% Inputs:
%    sample - An array of elements from which to sample.
%    n      - Number of elements to sample.
%
% Outputs:
%    randomSample - An array containing n randomly selected elements.

% Validate that n is not greater than the number of available elements
if n > numel(sample)
    error('n must be less than or equal to the number of elements in the sample.');
end

% Generate n unique random indices from 1 to the number of elements in sample
randomIndices = randperm(numel(sample), n);

% Retrieve the randomly sampled elements
randomSample = sample(randomIndices);

end
