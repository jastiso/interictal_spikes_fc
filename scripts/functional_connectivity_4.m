%% Get functional connectivity

clear
clc
close all
warning ON

addpath(genpath('/Users/stiso/Documents/Code/interictal_spikes_fc/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')

%%

% global variables
top_dir = '/Volumes/bassett-data/Jeni/RAM/';

% constants
win_length = 1; % in seconds
releases = ['1', '2', '3'];


for r = 1:numel(releases)
    release = releases(r);
    
    release_dir = [top_dir, 'release', release '/'];
    
    % for catching errors
    errors = struct('files', [], 'message', []);
    warnings = struct('files', [], 'message', []);
    
    % remove parent and hidden directories, then get protocols
    folders = dir([release_dir '/protocols']);
    folders = {folders([folders.isdir]).name};
    protocols = folders(cellfun(@(x) ~contains(x, '.'), folders));
    
    for p = 1:numel(protocols)
        protocol = protocols{p};
        
        % get global info struct
        fname = [release_dir 'protocols/', protocol, '.json'];
        fid = fopen(fname);
        raw = fread(fid);
        str = char(raw');
        fclose(fid);
        info = jsondecode(str);
        eval(['info = info.protocols.', protocol,';']);
        
        % get subjects
        subjects = fields(info.subjects);
        for s = 1:numel(subjects)
            subj = subjects{s};
            
            % save command window
            %clc
            if ~exist([top_dir, 'processed/release',release, '/', protocol, '/', subj, '/'], 'dir')
                mkdir([top_dir, 'processed/release',release, '/', protocol, '/', subj, '/']);
            end
            %eval(['diary ', [top_dir, 'processed/release',release, '/', protocol, '/', subj, '/log.txt']]);
            
            fprintf('\n******************************************\nStarting functional connectivity for subject %s...\n', subj)
            
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
                    save_dir = [top_dir, 'processed/release',release, '/', protocol, '/', subj, '/', exper, '/', sess, '/'];
                    img_dir = [top_dir, 'img/artifact/release',release, '/', protocol, '/', subj, '/', exper, '/', sess, '/'];
                    if ~exist(img_dir, 'dir')
                        mkdir(img_dir);
                    end
                    
                    % check if this subect has clean data
                    if exist([save_dir, 'data_clean.mat'], 'file')
                        load([save_dir, 'data_clean.mat'])
                        load([save_dir, 'header.mat'])
                        load([save_dir, 'channel_info.mat'])
                        
                        nElec = numel(ft_data.label);
                        nTrial = numel(ft_data.trial);
                        dur = cellfun(@(x) size(x,2), ft_data.trial);
                        
                        % prewhiten
                        
                        
                        % redefine trial
                        cfg = [];
                        cfg.length = win_length;
                        cfg.overlap = 0;
                        tmp = ft_redefinetrial(cfg,ft_data)
                        
                        % fc
                        % band limited
                        cfg = [];
                        cfg.method     = 'wavelet';
                        cfg.width      = 6;
                        cfg.output     = 'powandcsd';
                        cfg.foi        = freqs;
                        cfg.toi        = st:1/srate:en;
                        cfg.pad        = 'nextpow2';
                        cfg.trials     = t;
                        wave = ft_freqanalysis(cfg, data);
                        
                        % broadband
                        
                        
                        
                        % save things
                        
                        
                        
                        
                    end
                end
            end
        end
    end
end
