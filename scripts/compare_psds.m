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
                            
                            psds = zeros(nElec,nFreq);
                            aper = zeros(nElec);
                            for j = 1:nElec
                                spec = get_IRASA_spec(ft_data.trial{1}(j,:), 1, numel(ft_data.trial{1}(j,:)), ft_data.fsample, win_length_ir, step, filter);
                                psds(j,:) = mean(spec.mixd,2);
                                x = spec.freq;
                                y = mean(spec.frac,2);
                                P = polyfit(log10(x),log10(y),1);
                                aper(j) = P(1);
                            end
                            curr_clean = squeeze(mean(mean(psds)));
                            all_psd(cnt,:) = curr_clean;
                            all_aper(cnt) = mean(mean(aper));
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
