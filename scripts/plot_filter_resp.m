%% Plot filter responses
% make fake data
data = zeros(1,1000);
data(499) = 1;
srate = 500;
labels = {'1'};

ft_data = [];
ft_data.trial{1,1} = data;
ft_data.time{1,1} = (1/srate):(1/srate):(size(data,2)/srate);
ft_data.label = labels;
ft_data.fsample = srate;
ft_data.nSamples = size(data,2);
ft_data.sampleinfo = [1, size(data,2)];

cfg = [];
ft_data = ft_preprocessing(cfg,ft_data);


%% Notch filter

[b,a] = butter(4, [59/(ft_data.fsample/2), 61/(ft_data.fsample/2)], 'stop');
ft_data_notch = filtfilt(b,a,ft_data.trial{1}')';

figure(1); clf
subplot(1,3,1)
plot(ft_data.trial{1},'k', 'linewidth', 2);hold on
xlabel('Time (samples)')
ylabel('Signal')
plot(ft_data_notch,'b','linewidth',2); hold off
    
%% Each bandpass

subplot(1,3,2)
plot(ft_data.trial{1},'k', 'linewidth', 2);hold on
xlabel('Time (samples)')
ylabel('Signal')

freqs = unique(round(logspace(log10(4),log10(150),30)));
bands = [4, 8; 9, 15; 16 25; 36, 70; 71, 150];
for i = 1:size(bands,1)
    curr_range = bands(i,:);
    cfg = [];
    cfg.bpfilter = 'yes';
    cfg.bpfreq = curr_range;
    cfg.bpfiltdf = 0.5; % bandpass transition width
    cfg.bsfiltdf = 0.5; % bandstop transition width
    cfg.bpfiltdev = 0.01; % bandpass max passband deviation
    cfg.bsfiltdev = 0.05; % bandstp max passband deviation
    cfg.bpfilttype = 'firws'; % or 'firls' (slower), but avoid the default 'but' (= not well-suited for hilbert phase estimate)
    cfg.hilbert = 'complex';
    %cfg.plotfiltresp = 'yes';
    bp_data = ft_preprocessing(cfg, ft_data);


    plot(abs(bp_data.trial{1}),'linewidth',2)
end


%% Low pass

cfg = [];
cfg.lpfilt = 'yes';
cfg.lpfilttype = 'but';
cfg.lpfreq = 200;
cfg.lpfiltord = 4;

lfp = ft_preprocessing(cfg, ft_data);

subplot(1,3,3)
plot(ft_data.trial{1},'k', 'linewidth', 2);hold on
xlabel('Time (samples)')
ylabel('Signal')
plot(lfp.trial{1},'linewidth',2)
saveas(gca, ['/Users/stiso/Documents/Code/interictal_spikes_fc/img/filters.png'], 'png')

