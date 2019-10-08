function [] = plot_lfp(data,srate,i)
% Plot non-interactive trace of LFP, for use in parfor loops
% Inputs:
% data      channel x time
% srate     sample rate

%demean 
data = data - mean(data, 2);

% get offset for each channel
maxes = zeros(size(data,1),1) + mean(max(data, [], 2));
offset = zeros(size(maxes));
for i = 1:numel(offset)
    offset(i) = sum(maxes(i:end));
end
% add to data
data = data + offset;

% get x_vector
x = 1:size(data,2);
x = x./srate;

% plot
figure(i); clf
set(gcf, 'Position',  [100, 100, 4000, 500]);
plot(x,data, 'linewidth', 1.5)
xlabel('Time (s)')
ylim([min(min(data)), max(max(data))])

end

