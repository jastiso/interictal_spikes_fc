%% Count how much data each subject has
clear
clc
close all
warning ON

addpath(genpath('/Users/stiso/Documents/Code/interictal_spikes_fc/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')
addpath('/Users/stiso/Documents/MATLAB/arfit/')

%% Run this on a remote machine (main FC calculation)


% global variables
top_dir = '/Volumes/bassett-data/Jeni/RAM/';

releases = ['1', '2', '3'];

spike_win = 0.05; %for loading spike data
win_length = 1; % in seconds
detector = '';
nWin = 1;
cnt = 1;
n_wins = [];
n_wins_post = [];

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
        
        % subjects not to use
        load([top_dir, 'bad_datasets.mat'])
        
        % get subjects
        subjects = fields(info.subjects);
        for s = 1:numel(subjects)
            curr_w = 0;
            curr_w_post = 0;
            subj = subjects{s};
            
            % main function for functional connectivity - helps with paralelizing
            % do you want to save all your FC mtrices? faster if no
            save_flag = false;
            release_dir = [top_dir, 'release', release '/'];
            
            % get global info struct
            fname = [release_dir 'protocols/', protocol, '.json'];
            fid = fopen(fname);
            raw = fread(fid);
            str = char(raw');
            fclose(fid);
            info = jsondecode(str);
            eval(['info = info.protocols.', protocol,';']);
            
            % make subject directory
            subj_dir = [top_dir, 'FC/release',release, '/', protocol, '/', subj, '/'];
            if ~exist(subj_dir, 'dir')
                mkdir(subj_dir);
            end
            
            
            if ~exist([top_dir, 'processed/release',release, '/', protocol, '/', subj, '/'], 'dir')
                mkdir([top_dir, 'processed/release',release, '/', protocol, '/', subj, '/']);
            end
            
            % check that we need data for this subj
            if ~exist([top_dir, 'FC/release',release, '/', protocol, '/', subj, '/', 'win_', num2str(win_length), '/alt_spike', detector, '.csv'], 'file')
                fprintf('\n******************************************\nCounting for subject %s...\n', subj)
                
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
                        
                        if exist([data_dir, 'data_clean.mat'], 'file') && exist([data_dir, 'spike_info_', num2str(spike_win), '.mat'], 'file')
                            load([data_dir, 'data_clean.mat'])
                            load([data_dir, 'header.mat'])
                            
                            if ~isempty(ft_data.trial)
                                w = size(ft_data.trial{1},2)/500;
                                curr_w = curr_w + w;
                            end
                            
                            %check if this subect has clean data
                            reject = zeros(numel(ft_data.trial),1);
                            for i = 1:numel(ft_data.trial)
                                curr_ext = [subj, '_' exper, '_', sess, '_', num2str(i)];
                                reject(i) = any(strcmp(curr_ext, bad_datasets));
                            end
                            fprintf('\nRejected %d datasets\n', sum(reject))
                            
                            ft_data.trial = ft_data.trial(~reject);
                            ft_data.time = ft_data.time(~reject);
                            ft_data.sampleinfo = ft_data.sampleinfo(~reject,:);
                            
                            if ~isempty(ft_data.trial)
                                w = size(ft_data.trial{1},2)/500;
                                curr_w_post = curr_w_post + w;
                            end
                            
                        end
                    end
                end
            end
            n_wins = [curr_w, n_wins];
            n_wins_post = [curr_w_post, n_wins_post];
        end
    end
end

% print stuff
quantile(nonzeros(n_wins),[0 0.25 0.50 0.75 1])
quantile(n_wins_post(n_wins ~= 0),[0 0.25 0.50 0.75 1])