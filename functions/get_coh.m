function [C_band] = get_coh(data, bands)
% Gets coherence from fieldtrip mtmfft freqanalysis call.
% Inputs:
%   data        output of fieldtrip call to ft_freqnaalysis. Data should 
%               already be averaged over repetitions (tapers, trials, 
%               timepoints, etc). It should have fields crsspctrm,
%               powspctrm, freq, label, and labelcmb
%   bands       frequency boundaries for each band. Size Nx2 where N is the
%               number of bands. Band i spans from (bands(i,1), bands(i,2)
%               inclusive
% Outputs
%   C_bands     coherence, where Cxy = (abs(Pxy).^2)./(Pxx.*Pyy)

% constants
nPair = numel(data.labelcmb(:,1));
nFreq = numel(data.freq);
nTrial = numel(data.cumsumcnt);
nBand = size(bands, 1);

% initialize
C = zeros(nFreq, nPair, nTrial);
C_band = zeros(nBand, nPair, nTrial);

for i = 1:nPair
    % get csd
    Pxy = squeeze(data.crsspctrm(:,i,:));
    % get psd
    x = strcmp(data.labelcmb{i,1}, data.label);
    y = strcmp(data.labelcmb{i,2}, data.label);
    Pxx = squeeze(data.powspctrm(:,x,:));
    Pyy = squeeze(data.powspctrm(:,y,:));
    % get coherence
    C(:,i,:) = ((abs(Pxy).^2)./(Pxx.*Pyy))';
end

% average over freq bands
for i = 1:nBand
    idx = data.freq >= bands(i,1) & data.freq <= bands(i,2);
    C_band(i,:,:) = mean(C(idx,:,:));
end

end

