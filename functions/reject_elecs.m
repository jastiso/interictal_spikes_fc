function [rmv] = reject_elecs(data,thr, srate)
%identify bad elecs by kurtosis, weird PSD, or line length
%   this is good for elecs with large artifacts, but not for persistently
%   noisy elecs
% Inputs:
%       thr         the threshold for rejection in z-scores
%       data        data, nElec x nTime
%       srate       sampling rate
% Outputs
%       rmv         logical index of elecs to remove
% @author JStiso 07/2019


k = kurtosis(data');                % kurtosis of ts -- are there weird peaks? if so, maybe artifact?
zk = zscore(k);                     % zscore kurtosis

rm = zk >= thr;                % remove channels that show more kurtosis

[pxx,~] = pwelch(data',hanning(500),0, [], srate);    % calculate power spectral density w/ welchs method
d = squareform(pdist(log10(pxx)','spearman'));    % find average similarity of power spectral density plots
zd = zscore(mean(d));

rm = rm | (zd >= thr);

% now check line length (Esteller et al., 2001)
ll = sum(abs(diff(data,1,2)),2);
rm_ll = ll > 3*mean(ll);

rmv = rm' | rm_ll;
end

