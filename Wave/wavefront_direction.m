function [vx, vy] = wavefront_direction(pd,s)
% *WAVE*
%
% WAVEFRONT DIRECTION    calculates the vecocity vector 
%
% INPUT
% pd - phase gradient direction
% s -  instantaneous speed 
% 
%  direction of velocity = - spaitalmean(phase gradient)
%
% OUTPUT
% vx - the x component of the vector pointing in velocity direction
% vy - the y component of the vector pointing in velocity direction
%

assert( ( ndims(pd) >=2  ) , 'datacube inputs required' ); %& ( ndims(s) == 3 )
% assert( isequal( size(pd), size(s) ), 'datacube sizes must be equal' );

mag = -1;
% mag = -abs(s)/10; %in cm/s
vx = squeeze(mean(cos(pd),[1 2],'omitnan').*mag)' ;
vy = squeeze(mean(sin(pd),[1 2],'omitnan').*mag)';
end

