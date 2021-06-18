%Original time series
Fs = 1000;
t = 0:1/Fs:1-1/Fs;
x = cos(2*pi*100*t) + randn(size(t));
T = 1;
% psd
N = length(x);
xdft = fft(x);
phase = angle(xdft);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx = 2*psdx;
freq = 0:Fs/length(x):Fs/2;
% back to timeseries
amp = sqrt((1/2).*psdx)*(Fs);
z = amp.*exp(1i*phase);
x2 = ifft(z);
plot(t,x); hold on; plot(t, real(x2), '--k')

