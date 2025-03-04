function [PGD] = phase_gradient_directionality(pm,dx,dy)
% *WAVE*
%
%PHASE_GRADIENT_DIRECTIONALITY - measure of how well gradients align across
%                                the array as a function of time 
% Reference - Rubino, D., Robbins, K. & Hatsopoulos, N. 
%                   Propagating waves mediate information transfer in the motor cortex.
%                   Nat Neurosci 9, 1549–1557 (2006). https://doi.org/10.1038/nn1802
%
% INPUT
% pm - phase gradient magnitude
% dx - x component of the phase gradient 
% dy - y component of the phase gradient 
%
% OUTPUT
% pdg - phase gradient directionality time series 
%
num = sqrt(mean(dx,[1 2],'omitnan').^2 + mean(dy,[1 2],'omitnan').^2);
den = (2*pi).*mean(pm,[1 2],'omitnan');
PGD = squeeze(num./den)';
end

