clear
clc
close

% global variables and packages
addpath(genpath('/Users/stiso/Documents/Code/interictal_spikes_fc/'))
addpath(genpath('/Users/stiso/Documents/MATLAB/eeglab_current/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')

release = '1';
protocol = 'r1';

top_dir = '/Volumes/bassett-data/Jeni/RAM/';
release_dir = [top_dir, 'release', release '/'];
eval(['cd ', top_dir])

win = 0.05; % size of the window to look for the minimum number of channels, in seconds
discharge_tol=0.005; % taken from spike function
min_chan = 3;

%% Pick which subject

subj = 'R1045E';
exper = 'YC1';
sess = '1';

subj = 'R1031M';
exper = 'PAL2';
sess = '1';

subj = 'R1015J';
exper = 'FR1';
sess = '0';

subj = 'R1021D';
exper = 'YC1';
sess = '0';

subj = 'R1001P';
exper = 'FR2';
sess = '0';

subj = 'R1060M';
exper = 'FR2';
sess = '0';

%% Get spikes

%eval(['info = info.protocols.', protocol,';']);
out_clean = [];
marker = [];

if ~exist([top_dir, 'processed/release',release, '/', protocol, '/', subj, '/'], 'dir')
    mkdir([top_dir, 'processed/release',release, '/', protocol, '/', subj, '/']);
end

fprintf('******************************************\nStarting preprocessing for subject %s...\n', subj)



% get the path names for this session, loaded from a json file
%eval(['curr_info = info.subjects.' subj, '.experiments.' exper, '.sessions.x', sess, ';'])

% folders
raw_dir = [release_dir, 'protocols/', protocol, '/subjects/', subj, '/experiments/', exper, '/sessions/', sess, '/ephys/current_processed/'];
save_dir = [top_dir, 'processed/release',release, '/', protocol, '/', subj, '/', exper, '/', sess, '/'];
img_dir = [top_dir, 'img/diagnostics/release',release, '/', protocol, '/', subj, '/', exper, '/', sess, '/'];
if ~exist(img_dir, 'dir')
    mkdir(img_dir);
end
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

fprintf('\nExperiment %s session %s\n\n', exper, sess)

if exist([save_dir, 'data_clean.mat'], 'file')
    load([save_dir, 'data_clean.mat'])
    load([save_dir, 'header.mat'])
    load([save_dir, 'channel_info.mat'])
    
    % combine soz and interictal
    all_soz = unique([cell2mat(cellfun(@(x) find(strcmpi(x, ft_data.label)), soz, 'UniformOutput', false)); cell2mat(cellfun(@(x) find(strcmpi(x, ft_data.label)), interictal_cont, 'UniformOutput', false))]);
    nChan = numel(ft_data.label);
    
    if isempty(all_soz)
        warning('The interictal and SOZ channels were not marked for this subject')
    end
    
    nTrial = numel(ft_data.trial);
    for j = 1:nTrial
        % run detection alg
        
        [out,MARKER] = ...
            spike_detector_hilbert_v16_byISARG(ft_data.trial{j}', header.sample_rate);
        
        % eliminate some spikes
        nSamp = size(MARKER.d,1);
        nSpike = numel(out.pos);
        kept_spike = false(size(out.pos));
        
        for i = 1:nSpike
            curr_pos = out.pos(i);
            curr_chan = out.chan(i);
            
            if kept_spike(i) == 0
                win_spike = (out.pos > curr_pos & out.pos < (curr_pos + win));
                win_chan = out.chan(win_spike);
                
                if numel(unique(win_chan)) >= min_chan
                    kept_spike(win_spike) = true;
                end
                
            end
        end
        % select only good spikes
        out_clean(j).pos = out.pos(kept_spike);
        out_clean(j).dur = out.dur(kept_spike);
        out_clean(j).chan = out.chan(kept_spike);
        out_clean(j).weight = out.weight(kept_spike);
        out_clean(j).con = out.con(kept_spike);
        
        % get new M
        m_clean = zeros(size(MARKER.M));
        for i=1:size(out_clean(j).pos,1)
            m_clean(round(out_clean(j).pos(i)*MARKER.fs:out_clean(j).pos(i)*MARKER.fs+discharge_tol*MARKER.fs),...
                out_clean(j).chan(i))=out_clean(j).con(i);
        end
        marker(j).m_clean = m_clean;
        marker(j).d = MARKER.d;
        marker(j).fs = MARKER.fs;
        
    end
end




%% Plot everything

eegplot(marker(1).d', 'srate', marker(1).fs)
eegplot(marker(1).M', 'srate', marker(1).fs)
eegplot(marker(1).m_clean', 'srate', marker(1).fs)

st = 24*marker(1).fs;
en = 44*marker(1).fs;
figure(1); clf
plot_lfp(marker(1).d(st:en,:)', 24*marker(1).fs); hold on
plot_lfp(marker(1).d(logical(marker(1).m_clean(st:en,:)))', 24*marker(1).fs, 'r*');



%% Plot spikes

nElec = size(marker(1).d,2);

% plot length in samples
l = 4000;

% randomly pick a spike
st = size(marker(1).d,1);
while st > (size(marker(1).d,1) - l - 101)
    curr = randi(numel(out_clean(1).pos));
    st = round(out_clean(1).pos(curr)*marker(1).fs) - 100;
end

% get data ready for plotting
data = (marker(1).d - mean(marker(1).d))';

% get offset for each channel
maxes = zeros(size(data,1),1) + mean(max(data, [], 2));
offset = zeros(size(maxes));
for i = 1:numel(offset)
    offset(i) = sum(maxes(i:end));
end
% add to data
data = data + offset;

% get x_vector
x = 1:l;
x = x./marker(1).fs;

% offset markers
markers = marker(1).m_clean';
for i = 1:nElec
    markers(i,markers(i,:) ~= 0) = markers(i,markers(i,:) ~= 0) + offset(i);
end

% plot
figure(1); clf
set(gcf, 'Position',  [100, 100, 800, 2000]); 
plot(x,data(:,st:(st + l - 1)), 'k', 'linewidth', 1.5); hold on
xlabel('Time (s)')
ylim([min(min(data(:,st:(st + l - 1)))), max(max(data(:,st:(st + l - 1))))])
plot(x, markers(:,st:(st + l - 1)), 'rx', 'linewidth', 2)
hold off
saveas(gca, [top_dir, 'img/spike_examples/', subj, '_', num2str(curr),'.png'], 'png')





