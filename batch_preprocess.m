close all; clear all;

[g, s] = WIN_RS_params;

% Step one: Motion Corection
groovy_realign(g,s);

% Step two: Slice time correction
groovy_slice(g,s);

% Step three: Normalization
groovy_norm_calc(g,s);
groovy_norm_write(g,s);

% Step four: Smooth
%groovy_smooth(g,s);

% Step five: Band pass filter
%melodic_bpfilter(g,s);