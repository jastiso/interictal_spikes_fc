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



%% Plot single spike

l = 2000;
for i = 1:numel(chann_idx)
    curr = i;
    
    st = find(MARKER.M(:,curr) ~= 0);
    if ~isempty(st)
        st = st(1);
        st = st - 100;
    else
        st = 1;
    end
    
    figure(4); clf
    %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    subplot(2,1,1)
    plot(1:numel(st:(st+l)), MARKER.d(st:(st+l),curr), 'linewidth', 1.5); hold on
    plot(find(MARKER.M(st:(st+l),curr) ~= 0), MARKER.d(MARKER.M(st:(st+l),curr) ~= 0,curr), 'rx', 'linewidth', 2)
    title('Voltage')
    xlim([1,l])
    subplot(2,1,2)
    plot(1:numel(st:(st+l)), envelope(st:(st+l),curr), 'linewidth', 1.5); hold on
    plot(find(MARKER.M(st:(st+l),curr) ~= 0), envelope(MARKER.M(st:(st+l),curr) ~= 0,curr), 'rx', 'linewidth', 2)
    title('Envelope')
    xlim([1,l])
    hold off
    pause(0.5)
    %saveas(gca, [top_dir, 'img/spike_examples/',num2str(i),'.png'], 'png')
end




