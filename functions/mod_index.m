function [mi] = mod_index(binned_amp,uni_dist)
% Calculates the modulation index from Tort et al.
% this is really just a normalized KLDistance
%   Inputs are 
% binned_amp:   binned amplitudes, size nChannels x nBins
% uni:          a uniforn distribution of the same size to compare it to
% Outputs are
% mi:           modulation index, bounded 0-1

d = KLDiv(binned_amp, uni_dist);
mi = d./log10(size(binned_amp,2));

end

