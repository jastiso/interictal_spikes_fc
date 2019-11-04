function [wpli] = get_wpli(data)
% get wavelet weighted phaselocking value from fieldtrip freq data
%   Detailed explanation goes here



%% get dWPLI

cfg = [];
cfg.method = 'wpli_debiased';
wpli = ft_connectivityanalysis(cfg, wave);

end

