function [pac, mod_idx] = get_pac(phase,amp)
% get phase-amplitude coupling for oscillation in theta and broadband gamma
% get MI via KL distance method
% compare to surrogate (single cut), return z-value

%   Detailed explanation goes here

%% Initialize

nElec = numel(amp.label);
nPair = (nElec^2-nElec)/2;
nTrial = numel(amp.trial);
nBin = 18; % this parameter has been shown to not make a huge difference in MI
nSim = 100;
pac =  zeros(nTrial, nPair);

%% Get phase and amp

% from complex hilbert data, get phase
for t = 1:nTrial
    phase.trial{t} = atan2(imag(phase.trial{t}),real(phase.trial{t}));
end

% from complext hilbert data, get amp
for t = 1:nTrial
    amp.trial{t} = abs(amp.trial{t});
end

%% Get PAC

% Get modulation index
uni_dist = unifpdf(linspace(-pi, pi, nBin),-pi,pi);
mod_idx = zeros(nTrial,nPair);

% -pi to pi
angles = linspace(-pi, pi, 100);
[~,edges] = discretize(angles, nBin);

for i = 1:nTrial
    % vectorize
    curr_phase = phase.trial{i};
    curr_amp = amp.trial{i};
    
    % get indices for bins
    ind = discretize(curr_phase,edges);
    
    % get MI
    cnt = 0;
    for j = 1:nElec
        for k = (j+1):nElec
            cnt = cnt + 1;
            binned_amp = accumarray(ind(j,:)', curr_amp(k,:), [nBin, 1], @mean);
            % normalize over all bins
            binned_amp = binned_amp./sum(binned_amp,1);
            mod_idx(i,cnt) = mod_index(binned_amp', uni_dist);
            
            % get surrogate
            surrogate_dist = pac_surr_cut(nSim, curr_phase(j,:), curr_amp(k,:), nBin, edges);
            mi_surr = zeros(nBin, nSim);
            for m = 1:nSim
                mi_surr(:,m) = mod_index(surrogate_dist(:,m)', uni_dist);
            end
            
            % z-score
            pac(i, cnt) = (mod_idx(i,cnt) - mean(mi_surr))/std(mi_surr);
        end
    end
end


end

