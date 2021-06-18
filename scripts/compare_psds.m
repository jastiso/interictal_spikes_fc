%% Get PSDs for IED and normal windows

clear
clc
close all
warning ON

restoredefaultpath
addpath(genpath('/Users/stiso/Documents/Code/interictal_spikes_fc/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')
addpath(genpath('/Users/stiso/Documents/MATLAB/IRASA/'))

% global variables
top_dir = '/Volumes/bassett-data/Jeni/RAM/';

releases = ['1', '2', '3'];

spike_win = 0.05; %for loading spike data
win_length = 1; % in seconds
detector = '';
win_length_ir = 400; % in ms
step = 50; %in ms
filter = 0; % without lowpass alaising filter

all_psd = [];
all_ied = [];
all_aper = [];
all_aperied = [];
cnt = 1;

% global variables
top_dir = '/Volumes/bassett-data/Jeni/RAM/';

releases = ['1', '2', '3'];

spike_win = 0.05; %for loading spike data
win_length = 1; % in seconds
detector = '_delphos_auto';

%bands
freqs = unique(round(logspace(log10(4),log10(150),30)));
bands = [4, 8; 9, 15; 16 25; 36, 70];
band_names = [{'theta'}, {'alpha'}, {'beta'}, {'gamma'}];
% constants
nBand = size(bands,1);
all_osc = {};
cnt = 1;
nFreq = 34; % number of frequencies in the IRASA

% subjects not to use
load([top_dir, 'bad_datasets.mat'])
errors = struct('files', [], 'message', []);

for r = 1:numel(releases)
    release = releases(r);
    
    release_dir = [top_dir, 'release', release '/'];
    
    % remove parent and hidden directories, then get protocols
    folders = dir([release_dir '/protocols']);
    folders = {folders([folders.isdir]).name};
    protocols = folders(cellfun(@(x) ~contains(x, '.'), folders));
    
    % just for testing that qsub fuction will work
    for p = 1:numel(protocols)
        protocol = protocols{p};
        
        % get global info struct
        fname = [release_dir 'protocols/', protocol, '.json'];
        fid = fopen(fname);
        raw = fread(fid);
        str = char(raw');
        fclose(fid);
        info = jsondecode(str);
        info = info.protocols.r1;
        
        % get subjects
        subjects = fields(info.subjects);
        for s = 1:numel(subjects)
            subj = subjects{s};
            
            % constants
            nBand = size(bands,1);
            
            release_dir = [top_dir, 'release', release '/'];
            
            % get global info struct
            fname = [release_dir 'protocols/', protocol, '.json'];
            fid = fopen(fname);
            raw = fread(fid);
            str = char(raw');
            fclose(fid);
            info = jsondecode(str);
            eval(['info = info.protocols.', protocol,';']);
            
            % subjects not to use
            load([top_dir, 'bad_datasets.mat'])
            errors = struct('files', [], 'message', []);
            
            % make subject directory
            subj_dir = [top_dir, 'FC/release',release, '/', protocol, '/', subj, '/'];
            if ~exist(subj_dir, 'dir')
                mkdir(subj_dir);
            end
            
            if ~exist([top_dir, 'processed/release',release, '/', protocol, '/', subj, '/'], 'dir')
                mkdir([top_dir, 'processed/release',release, '/', protocol, '/', subj, '/']);
            end
            
            fprintf('\n******************************************\nStarting psd comparison for subject %s...\n', subj)
            
            % get experiements
            eval(['experiments = fields(info.subjects.' subj, '.experiments);'])
            for e = 1:numel(experiments)
                exper = experiments{e};
                
                % get seesions
                eval(['sessions = fields(info.subjects.' subj, '.experiments.', exper, '.sessions);'])
                for n = 1:numel(sessions)
                    
                    sess = sessions{n};
                    sess = strsplit(sess, 'x');
                    sess = sess{end};
                    % get the path names for this session, loaded from a json file
                    eval(['curr_info = info.subjects.' subj, '.experiments.' exper, '.sessions.x', sess, ';'])
                    
                    % folders
                    data_dir = [top_dir, 'processed/release',release, '/', protocol, '/', subj, '/', exper, '/', sess, '/'];
                    save_dir = [top_dir, 'FC/release',release, '/', protocol, '/', subj, '/', exper, '/', sess, '/'];
                    img_dir = [top_dir, 'img/FC/release',release, '/', protocol, '/', subj, '/', exper, '/', sess, '/'];
                    if ~exist(img_dir, 'dir')
                        mkdir(img_dir);
                    end
                    if ~exist(save_dir, 'dir')
                        mkdir(save_dir);
                    end
                    if ~exist([subj_dir, 'win_', num2str(win_length), '/'], 'dir')
                        mkdir([subj_dir, 'win_', num2str(win_length), '/']);
                    end
                    
                    if exist([data_dir, 'data_clean.mat'], 'file') && (exist([data_dir, 'spike_info', detector, '.mat'], 'file') || exist([data_dir, 'spike_info', num2str(win), '.mat'], 'file'))
                        load([data_dir, 'data_clean.mat'])
                        load([data_dir, 'header.mat'])
                        load([data_dir, 'channel_info.mat'])
                        if strcmp(detector, '')
                            load([data_dir, 'spike_info_', num2str(spike_win), '.mat'])
                        else
                            load([data_dir, 'spike_info', detector, '.mat'])
                        end
                        load([data_dir, 'artifact.mat'])
                        load([data_dir, 'demographics.mat'])
                        %try
                        % check if this subect has clean data
                        reject = zeros(numel(ft_data.trial),1);
                        for i = 1:numel(ft_data.trial)
                            curr_ext = [subj, '_' exper, '_', sess, '_', num2str(i)];
                            reject(i) = any(strcmp(curr_ext, bad_datasets));
                        end
                        fprintf('\nRejected %d datasets\n', sum(reject))
                        
                        ft_data.trial = ft_data.trial(~reject);
                        ft_data.time = ft_data.time(~reject);
                        out_clean = out_clean(~reject);
                        artifact_all = artifact_all(~reject);
                        ft_data.sampleinfo = ft_data.sampleinfo(~reject,:);
                        
                        if ~isempty(ft_data.trial)
                            
                            nElec = numel(ft_data.label);
                            nPair = (nElec^2-nElec)/2;
                            
                            % constants
                            lower_tri = reshape(tril(true(nElec),-1),[],1);
                            
                            % get window start times
                            trl = [];
                            spike_index = [];
                            spike_spread = [];
                            spike_num = [];
                            spike_chan = {};
                            time_vec = [];
                            cnt = 1;
                            for i = 1:numel(ft_data.trial)
                                idx = 1;
                                curr_data = ft_data.trial{i};
                                curr_spike = out_clean(i);
                                curr_artifact = artifact_all(i);
                                dur = size(curr_data,2);
                                trl_offset = ft_data.sampleinfo(i,1);
                                while (idx + (win_length*header.sample_rate)) <= dur
                                    st = round(idx); % gets rid of scientific notation
                                    en = round(st + (win_length*header.sample_rate));
                                    st_ms = st/header.sample_rate;
                                    en_ms = en/header.sample_rate;
                                    spike_flag = 0;
                                    % check if there is an artifact
                                    if ~any(curr_artifact.idx(st:en))
                                        % record if there is a spike, if so
                                        % move idx up
                                        if any((curr_spike.pos >= st_ms) & (curr_spike.pos <= en_ms))
                                            spike_flag = 1;
                                            % if there are spikes in this
                                            % window, move the start to the
                                            % first spike in the window
                                            st = min(curr_spike.pos((curr_spike.pos >= st_ms) & (curr_spike.pos <= en_ms)))*header.sample_rate - 1;
                                            en = st + (win_length*header.sample_rate - 1);
                                        end
                                        % check that we havent gone
                                        % past the end of the data
                                        if en <= dur
                                            % get all the spikes in the window
                                            curr_idx = (curr_spike.pos >= st_ms) & (curr_spike.pos <= en_ms);
                                            seqs = unique(curr_spike.seq(curr_idx));
                                            % update spike idx
                                            spike_index(cnt) = spike_flag;
                                            spike_num(cnt) = numel(unique(curr_spike.seq(curr_idx)));
                                            if spike_flag
                                                spread = zeros(spike_num(cnt), 1);
                                                for m = 1:numel(spread)
                                                    spread(m) = numel(unique(curr_spike.chan(curr_spike.seq == seqs(m) & curr_idx)));
                                                end
                                                spike_spread(cnt) = mean(spread);
                                            else
                                                spike_spread(cnt) = 0;
                                            end
                                            spike_chan(cnt) = {curr_spike.chan(curr_idx)};
                                            % update trl
                                            trl(cnt,:) = [st + (trl_offset - 1), en + (trl_offset - 1), 0];
                                            % add time
                                            time_vec(cnt) = (st + trl_offset - 1)/header.sample_rate;
                                            % update cnt
                                            cnt = cnt + 1;
                                        end
                                    end
                                    idx = en + 1;
                                end
                            end
                        end
                        
                        % check for oscillations
                        % check that we found at least one window with
                        % a spike
                        if sum(spike_num) > 0
                            %try
                                % redefine trial
                                cfg = [];
                                cfg.trl = round(trl);
                                win_data = ft_redefinetrial(cfg,ft_data);
                                nTrial = numel(win_data.trial);
                                clear ft_data artifact_all
                                
                                % prewhiten
                                cfg = [];
                                cfg.derivative = 'yes';
                                ft_preprocessing(cfg, win_data);
                                
                                
                                psds = zeros(nTrial,nElec,nFreq);
                                aper = zeros(nTrial,nElec);
                                for i = 1:nTrial
                                    for j = 1:nElec
                                        spec = get_IRASA_spec(win_data.trial{i}(j,:), 1, numel(win_data.trial{i}(j,:)), win_data.fsample, win_length_ir, step, filter);
                                        psds(i,j,:) = mean(spec.mixd,2);
                                        x = spec.freq;
                                        y = mean(spec.frac,2);
                                        P = polyfit(log10(x),log10(y),1);
                                        aper(i,j) = P(1);
                                    end
                                end
                                curr_clean = squeeze(mean(mean(psds(~spike_index,:,:))));
                                curr_ied = squeeze(mean(mean(psds(logical(spike_index),:,:))));
                                all_psd(cnt,:) = curr_clean;
                                all_ied(cnt,:) = curr_ied;
                                all_aper(cnt) = mean(mean(aper(~spike_index,:)));
                                all_aperied(cnt) = mean(mean(aper(logical(spike_index),:)));
                                cnt = cnt + 1;
%                             catch
%                                 fprintf("\nError for subj %s", subj)
%                             end
                        end
                    end
                end
            end
        end
    end
end