function [f,pxx] = fft_plot(x,win_len, srate)
% makes an FFT plot
% function treats the columns as independent channels
% win_length is the length of the sample for each estimate

% Note that this is mostly for plotting, I would not use these parameters
% to obtain power estimates for analysis

[pxx,f] = pwelch(x,hanning(win_len),0, [], srate);

plot(f,10*log10(pxx), 'linewidth', 2)
title('FFT')
xlabel('Frequency (Hz)')
ylabel('Power (dB)')

end

