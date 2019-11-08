function [surrogate_dist] = pac_surr_cut(nSim, phase, amp, nBin, edges)
% Get surrogate data for PAC analysis using the single sut method (Aru et
% al). This method workds across epoched data
%   Inputs
% nSim:     The number of simulations for the surrogate distribution
% min_cut:  The minimum index for a cut. Intuitive, you probably dont want
%           to test a cut of 0, this is the same
% phase:    The vectorized phase - a vector of size nTimePt x 1
% amp:      The vectorized amplitude - a vector of size nTimePt x 1
% nBin:     the number of bins to but the phase into
% edges:    edges of bins, obtained by descretizing the range -pi:pi

% Outputs
% surrogate_dist:   the normalized distribution of amplitudes for each simulation

nTime = numel(phase);
surrogate_dist = zeros(nBin,nSim);

for n = 1:nSim
    % repeat process with single random cut through the trials (keep at
    % least min_cut together before cutting)
    % vectorize
    cut_ind = 1;
    while cut_ind == 1 % you'll get an error if this equals 1
        cut_ind = randi(nTime,1);
    end
    
    phase_vect_surr = phase([cut_ind:nTime, 1:(cut_ind-1)]);
    
    % get indices for bins
    ind_surr = discretize(phase_vect_surr,edges);
    
    % get avg amp
    surrogate_dist(:, n) = accumarray(ind_surr', amp, [nBin, 1], @mean);
    surrogate_dist(:, n) = surrogate_dist(:, n)./sum(surrogate_dist(:, n));
end
end

