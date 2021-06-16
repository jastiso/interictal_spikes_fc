function [spec] = get_IRASA_spec(data, exp_st, exp_en, srate, win_length_ms, step_ms, filter)
% This helper function makes it easier to run IRASA for different
% parameters

% data should be a VECTOR for a specific electrode
% times should be in MILISECONDS

win_length = floor((win_length_ms)*(srate/1000)) - 1; % in samples
step = floor(step_ms*(srate/1000));
nWin = floor(((exp_en-exp_st)*(srate/1000) - win_length)/step); % (length - win)/step

sig = zeros(win_length + 1, nWin);
for i = 1:nWin
    st = ceil((exp_st*(srate/1000)) + (i-1)*step + 1);
    en = ceil(st + win_length);
    
    sig(:,i) = data(1,st:en);
end

% spec = amri_sig_fractal(sig,srate,...); optional arguments frange,
% detrend, filter
% make frange short to avoid knee in PSD
spec = amri_sig_fractal(sig, srate, 'frange', [2, 70], 'detrend', 1, 'filter', filter);


end

