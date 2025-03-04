function [s] = speedSpatial( wt, pm )

% *WAVE*
%
% SPEED     makes the instantaneous speed measurement
%                           via the analytic signal framework
%
% INPUT
% wt - instantaneous frequency datacube (r,c,t)
% pm - phase gradient magnitude datacube (r,c,t)
%
% OUTPUT
% s - instantaneous speed estimate at t spatiallly average
%

assert( ( ndims(wt) == 3 ) & ( ndims(pm) == 3 ), 'datacube inputs required' );
assert( isequal( size(wt), size(pm) ), 'datacube sizes must be equal' );

s = 0.1*squeeze(mean(wt,[1 2],'omitnan') ./ mean(abs(pm),[1 2],'omitnan'))'; % in cm/s