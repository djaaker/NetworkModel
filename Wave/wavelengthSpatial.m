function [l] = wavelengthSpatial( pm )

% *WAVE*
%
% WAVELENGTH     makes the instantaneous wavelength measurement
%                           via the analytic signal framework
%
%   Reference: https://github.com/mullerlab/davis2021ncomms/blob/main/external/gradient_measure.m
%
% INPUT
% pm - phase gradient magnitude datacube (r,c,t)
%
% OUTPUT
% l - instantaneous wavelength estimate at t spatially average
%

assert(( ndims(pm) == 3 ), 'datacube inputs required' );
l = mean((0.1./pm),[1 2],'omitnan'); % in cm/s